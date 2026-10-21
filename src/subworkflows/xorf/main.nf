#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xorf

    Main subworkflow: chunks input, gets candidates, predicts ORFs, concatenates results
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller

    GitHub:  https://github.com/hillerlab/xorf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHUNKER }      from '../../modules/chunker/main.nf'
include { CHUNKER as UNMASKED_CHUNKER } from '../../modules/chunker/main.nf'

include { CONCAT }       from '../../modules/concat/main.nf'
include { CONCAT as CONCAT_RAW }   from '../../modules/concat/main.nf'
include { CONCAT as CONCAT_RENAMED }   from '../../modules/concat/main.nf'
include { CONCAT as CONCAT_RENAMED_RAW }   from '../../modules/concat/main.nf'

include { GENOMEMASK_SELENO } from '../../modules/genomemask/seleno/main.nf'
include { GENEPRED_LINT } from '../../modules/genepred/lint/main.nf'
include { DETACH_DUPLICATES } from '../../modules/detach/main.nf'
include { ISOTOOLS_TRUNCATION_DETECTOR } from '../../modules/isotools/utr/main.nf'
include { STRIP_OCCURRENCES as STRIP_TRUNCATIONS } from '../../modules/strip/main.nf'

include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_UNMASKED } from '../../modules/bedtools/intersect/main.nf'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_MASKED } from '../../modules/bedtools/intersect/main.nf'
include { BEDTOOLS_EXCLUDE as BEDTOOLS_EXCLUDE_UNMASKED } from '../../modules/bedtools/exclude/main.nf'
include { BEDTOOLS_EXCLUDE as BEDTOOLS_EXCLUDE_MASKED } from '../../modules/bedtools/exclude/main.nf'

include { RENAME_PREDICTIONS as RENAME_PREDICTIONS_RAW } from '../../modules/rename/main.nf'
include { RENAME_PREDICTIONS } from '../../modules/rename/main.nf'

include { PREDICT_ORFS } from '../predict_orfs/main.nf'
include { GET_CANDIDATES } from '../candidates/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XORF {
    take:
      regions        // [ [id:id, chr:chr] , file ]
      sequence       // [ file ]
      database       // [ file ]
      output_dir     // path
      chunk_size     // int
      samba_weights  // channel: [ meta, path ]
      esm_weights          // channel: path (HF hub cache dir, or dummy when !params.esm)
      predict_keep_raw     // boolean
      selenocysteine_sites // path
      skip_netstart        // boolean
      skip_transaid        // boolean
      rename_deactivate    // boolean
      do_polishing         // boolean
      strip_truncations    // boolean
      skip_joined_concat   // boolean
      run_only_on          // boolean
      run_only_mode        // string [ options: mask, unmask ]
      run_only_target      // string [ options: intersect, exclude ]
      database_versions    // channel: versions.yml from the protein database preparation steps

    main:
      def ch_regions  = regions
      def ch_sequence = sequence
      def ch_database = database
      def chunkSize   = chunk_size ?: 20

      // Use deterministic hash for reproducibility; previously randomHash() caused
      // non-deterministic output file names and downstream ordering.
      def localHash = params.prefix ? params.prefix.hashCode().toString().replaceAll("-", "").take(6) : "xorf"

      def ch_versions = Channel.empty()

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          LINTING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      GENEPRED_LINT(
        ch_regions
      )

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHUNKING OWNED CHANNELS
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_masked_chunked_regions = Channel.empty()
      ch_masked_chunked_sequences = Channel.empty()

      ch_unmasked_chunked_regions = Channel.empty()
      ch_unmasked_chunked_sequences = Channel.empty()

      ch_rescue_chunked_regions = Channel.empty()
      ch_rescue_chunked_sequences = Channel.empty()

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHUNKING + ROUTING RUNNING MODE [ --run_only_on true/false ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: run_only_on is able to neglect either masked/unmasked channels
      // in order to run only on masked or unmasked sequences
      if (run_only_on) {
          ch_chunker_target = Channel.empty()

          switch (run_only_mode) {
              case 'mask':
                if (selenocysteine_sites) {
                  Channel.fromPath(selenocysteine_sites, checkIfExists: true)
                    .map { it -> [ [id: it.baseName], it ] }
                    .set { ch_selenocysteine_sites }

                  GENOMEMASK_SELENO(
                      ch_sequence.map { it -> [ [id: it.baseName], it ] },
                      ch_selenocysteine_sites
                  )

                  switch (run_only_target) {
                      case 'intersect':
                        BEDTOOLS_INTERSECT_MASKED(
                            ch_regions,
                            ch_selenocysteine_sites
                        )

                        ch_chunker_target = BEDTOOLS_INTERSECT_MASKED.out.bed
                        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_MASKED.out.versions)

                      break

                      case 'exclude':
                        BEDTOOLS_EXCLUDE_MASKED(
                            ch_regions,
                            ch_selenocysteine_sites
                        )

                        ch_chunker_target = BEDTOOLS_EXCLUDE_MASKED.out.bed
                        ch_versions = ch_versions.mix(BEDTOOLS_EXCLUDE_MASKED.out.versions)
                      break

                      default:
                        error """
                        ERROR: run_only_target is not recognized.
                        Please provide a valid run_only_target: 'intersect' or 'exclude'.
                        """.stripIndent()
                  }

                  CHUNKER(
                      ch_chunker_target
                        .map { meta, file -> [ meta + [ masked: true, hash: localHash ], file ] },
                      GENOMEMASK_SELENO.out.fasta.first(),
                      chunkSize,
                  )

                  ch_masked_chunked_regions = CHUNKER.out.chunked_regions
                  ch_masked_chunked_sequences = CHUNKER.out.chunked_sequences

                  ch_rescue_chunked_regions = CHUNKER.out.rescue_regions
                    .map { meta, file -> [ meta + [ upstream: 0, downstream: 0 ], file ] }
                  ch_rescue_chunked_sequences = CHUNKER.out.rescue_sequences
                    .map { meta, file -> [ meta + [ upstream: 0, downstream: 0 ], file ] }

                  ch_versions = ch_versions.mix(GENOMEMASK_SELENO.out.versions)
                } else {
                  error """
                  ERROR: run_only_on is true but no selenocysteine_sites were provided.
                  Please provide a selenocysteine_sites file or set run_only_on to false.
                  """.stripIndent()
                }
              break

              case 'unmask':
                if (selenocysteine_sites) {
                  switch (run_only_target) {
                      case 'intersect':
                        BEDTOOLS_INTERSECT_UNMASKED(
                            ch_regions,
                            ch_selenocysteine_sites
                        )

                        ch_chunker_target = BEDTOOLS_INTERSECT_UNMASKED.out.bed
                        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_UNMASKED.out.versions)

                      break

                      case 'exclude':
                        BEDTOOLS_EXCLUDE_UNMASKED(
                            ch_regions,
                            ch_selenocysteine_sites
                        )

                        ch_chunker_target = BEDTOOLS_EXCLUDE_UNMASKED.out.bed
                        ch_versions = ch_versions.mix(BEDTOOLS_EXCLUDE_UNMASKED.out.versions)
                      break

                      default:
                        error """
                        ERROR: run_only_target is not recognized.
                        Please provide a valid run_only_target: 'intersect' or 'exclude'.
                        """.stripIndent()
                  }

                  UNMASKED_CHUNKER(
                      ch_chunker_target
                        .map { meta, file -> [ meta + [ chr:'UNMSK', masked: false, hash: localHash ], file ] },
                      ch_sequence.map { it -> [ [ id: it.baseName ], it ] },
                      chunkSize,
                  )

                  ch_unmasked_chunked_regions = UNMASKED_CHUNKER.out.chunked_regions
                  ch_unmasked_chunked_sequences = UNMASKED_CHUNKER.out.chunked_sequences

                  ch_rescue_chunked_regions = UNMASKED_CHUNKER.out.rescue_regions
                    .map { meta, file -> [ meta + [ upstream: 0, downstream: 0 ], file ] }
                  ch_rescue_chunked_sequences = UNMASKED_CHUNKER.out.rescue_sequences
                    .map { meta, file -> [ meta + [ upstream: 0, downstream: 0 ], file ] }

                } else {
                  error """
                  ERROR: run_only_on is true but no selenocysteine_sites were provided.
                  Please provide a selenocysteine_sites file or set run_only_on to false.
                  """.stripIndent()
                }

              break

              default:
                error """
                ERROR: run_only_mode is not recognized.
                Please provide a valid run_only_mode: 'mask' or 'unmask'.
                """.stripIndent()
              break
          }
      } else {
          if (selenocysteine_sites) {
              Channel.fromPath(selenocysteine_sites, checkIfExists: true)
                .map { it -> [ [id: it.baseName], it ] }
                .set { ch_selenocysteine_sites }

              GENOMEMASK_SELENO(
                  ch_sequence.map { it -> [ [id: it.baseName], it ] },
                  ch_selenocysteine_sites
              )

              CHUNKER(
                  ch_regions.map { meta, file -> [ meta + [ masked: true, hash: localHash ], file ] },
                  GENOMEMASK_SELENO.out.fasta.first(),
                  chunkSize,
              )

              ch_masked_chunked_regions = CHUNKER.out.chunked_regions
              ch_masked_chunked_sequences = CHUNKER.out.chunked_sequences

              BEDTOOLS_INTERSECT_UNMASKED(
                  ch_regions,
                  ch_selenocysteine_sites
              )
              UNMASKED_CHUNKER(
                  BEDTOOLS_INTERSECT_UNMASKED.out.bed
                    .map { meta, file -> [ meta + [ chr:'UNMSK', masked: false, hash: localHash ], file ] },
                  ch_sequence.map { it -> [ [ id: it.baseName ], it ] },
                  chunkSize,
              )
              ch_unmasked_chunked_regions = UNMASKED_CHUNKER.out.chunked_regions
              ch_unmasked_chunked_sequences = UNMASKED_CHUNKER.out.chunked_sequences

              ch_rescue_chunked_regions = CHUNKER.out.rescue_regions
                .concat(UNMASKED_CHUNKER.out.rescue_regions)
                .map { meta, file -> [ meta + [ upstream: 0, downstream: 0 ], file ] }
              ch_rescue_chunked_sequences = CHUNKER.out.rescue_sequences
                .concat(UNMASKED_CHUNKER.out.rescue_sequences)
                .map { meta, file -> [ meta + [ upstream: 0, downstream: 0 ], file ] }

              ch_versions = ch_versions.mix(GENOMEMASK_SELENO.out.versions)
              ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_UNMASKED.out.versions)
          } else {
              CHUNKER(
                  ch_regions,
                  ch_sequence.map { it -> [ [ id: it.baseName, hash: localHash ], it ] },
                  chunkSize,
              )

              ch_masked_chunked_regions = CHUNKER.out.chunked_regions
              ch_masked_chunked_sequences = CHUNKER.out.chunked_sequences

              ch_rescue_chunked_regions = CHUNKER.out.rescue_regions
                .map { meta, file -> [ meta + [ upstream: 0, downstream: 0 ], file ] }
              ch_rescue_chunked_sequences = CHUNKER.out.rescue_sequences
                .map { meta, file -> [ meta + [ upstream: 0, downstream: 0 ], file ] }
          }
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          JOINING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_masked_chunked_regions
          .concat(ch_unmasked_chunked_regions)
          .concat(ch_rescue_chunked_regions)
          .flatMap { meta, region -> 
              def regions = region instanceof List ? region : [region]
              regions.collect { it ->
                [ meta + [ name: it.baseName ], it] }
              }
          .join(
              ch_masked_chunked_sequences
                .concat(ch_unmasked_chunked_sequences)
                .concat(ch_rescue_chunked_sequences)
                .flatMap { meta, fasta -> 
                    def fas = fasta instanceof List ? fasta : [fasta]
                    fas.collect { it ->
                      [ meta + [ name: it.baseName ], it] }
                }
          )
          .set { ch_pairs }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PREDICTION
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      GET_CANDIDATES(
          ch_pairs,
          ch_database,
          samba_weights,
          esm_weights,
          skip_netstart,
          skip_transaid
      )

      PREDICT_ORFS(
          GET_CANDIDATES.out.candidates,
          GET_CANDIDATES.out.counts,
          predict_keep_raw
      )

      PREDICT_ORFS.out.orfs
          .filter { meta, bed, tsv -> bed.size() > 0 }
          .map { meta, bed, tsv -> 
              def groupKey = "${meta.id}@${meta.chr}"
              tuple(groupKey, meta, bed, tsv)
          }
          .groupTuple()
          .filter { groupKey, metas, beds, tsvs -> !beds.isEmpty() }
          .map { groupKey, metas, beds, tsvs ->
              return tuple([ id: groupKey, name: metas[0].id, chr: metas[0].chr, masked: metas[0].masked, hash: localHash ], beds, tsvs)
          }
          .set { ch_all }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          RENAMING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_raw_renamed = Channel.empty()
      if (predict_keep_raw) {
          PREDICT_ORFS.out.raw
              .filter { meta, bed, tsv -> bed.size() > 0 }
              .map { meta, bed, tsv -> 
                  def groupKey = "${meta.id}@${meta.chr}"
                  tuple(groupKey, meta, bed, tsv)
              }
              .groupTuple()
              .filter { groupKey, metas, beds, tsvs -> !beds.isEmpty() }
              .map { groupKey, metas, beds, tsvs ->
                  return tuple([ name: metas[0].id, chr: metas[0].chr, id: groupKey, masked: metas[0].masked, hash: localHash ], beds, tsvs)
              }
              .set { ch_raw }

          CONCAT_RAW(
              ch_raw
          )

          if (!rename_deactivate) {
              RENAME_PREDICTIONS_RAW(
                  CONCAT_RAW.out.bed,
                  CONCAT_RAW.out.tsv
              )

             RENAME_PREDICTIONS_RAW.out.files
              .map { metas, beds, tsvs ->
                  [beds, tsvs, metas]
              }
              .collect(flat: false)
              .map { pairs ->
                  def bedFiles = pairs.collectMany { pair ->
                      pair[0] instanceof Collection ? pair[0] : [pair[0]]
                  }

                  def tsvFiles = pairs.collectMany { pair ->
                      pair[1] instanceof Collection ? pair[1] : [pair[1]]
                  }

                  def metadataFiles = pairs.collectMany { pair ->
                      pair[2] instanceof Collection ? pair[2] : [pair[2]]
                  }

                  def newId = metadataFiles[0].id 
                  return [
                      [ id: newId, name: localHash ],
                      bedFiles,
                      tsvFiles
                  ]
              }
              .set { ch_raw_renamed }
          }
      }
      
      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CONCATENATION
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      CONCAT(
          ch_all
      )

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          RENAMING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_full_length_with_duplicates = Channel.empty()
      ch_full_length_predictions = Channel.empty()
      ch_truncations = Channel.empty()
      ch_duplicates  = Channel.empty()
      if (!rename_deactivate) {
          RENAME_PREDICTIONS(
              CONCAT.out.bed,
              CONCAT.out.tsv
          )

          if (!skip_joined_concat) {
            RENAME_PREDICTIONS.out.files
                .map { metas, beds, tsvs ->
                    [beds, tsvs, metas]
                }
                .collect(flat: false)
                .map { pairs ->
                    def bedFiles = pairs.collectMany { pair ->
                        pair[0] instanceof Collection ? pair[0] : [pair[0]]
                    }

                    def tsvFiles = pairs.collectMany { pair ->
                        pair[1] instanceof Collection ? pair[1] : [pair[1]]
                    }

                    def metadataFiles = pairs.collectMany { pair ->
                        pair[2] instanceof Collection ? pair[2] : [pair[2]]
                    }

                    def newId = metadataFiles[0].id
                    return [
                        [ id: newId, name: localHash ],
                        bedFiles,
                        tsvFiles
                    ]
                }
                .set { ch_renamed }

            // INFO: process uses 
            // > [meta.id, meta.name].findAll { it }.join('.') as new prefix
            CONCAT_RENAMED(
                ch_renamed
            )
            ch_full_length_with_duplicates = CONCAT_RENAMED.out.files.map { meta, bed, tsv -> tuple(meta, bed) }
            ch_full_length_predictions = CONCAT_RENAMED.out.files

            CONCAT_RENAMED_RAW(
                ch_raw_renamed
            )
          } else {
            ch_full_length_with_duplicates = RENAME_PREDICTIONS.out.bed
            ch_full_length_predictions = RENAME_PREDICTIONS.out.files
          }
      } else {
          ch_full_length_with_duplicates = CONCAT.out.files.map { meta, bed, tsv -> tuple(meta, bed) }
          ch_full_length_predictions = CONCAT.out.files
      }

      if (do_polishing) {
          /*
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              DUPLICATES
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          */

          DETACH_DUPLICATES( ch_full_length_with_duplicates )
          ch_full_length_transcripts = DETACH_DUPLICATES.out.pass
          ch_duplicates = DETACH_DUPLICATES.out.duplicates

          /*
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              TRUNCATION
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          */

          ISOTOOLS_TRUNCATION_DETECTOR(
              ch_full_length_transcripts
          )
          // INFO: strip_truncations false (annotation-grade input) -> a hint no descriptor
          // INFO: line contains, so every prediction stays in 04_results
          STRIP_TRUNCATIONS(
              ch_full_length_transcripts,
              ISOTOOLS_TRUNCATION_DETECTOR.out.descriptor,
              strip_truncations ? "TRUNCATED" : "KEEP_ALL_PREDICTIONS"
          )

          // `files` used to be CONCAT_RENAMED (pre-detach, pre-strip). Downstream
          // consumers then raced polish: they started as soon as concat closed,
          // while DETACH / iso-utr / STRIP ran as unused side tasks. Re-point
          // the emit at stripped HQ so the DAG waits on 04_results.
          ch_pre_polish_tsv = ch_full_length_predictions.map { meta, _bed, tsv ->
              tuple(meta, tsv)
          }
          ch_full_length_predictions = STRIP_TRUNCATIONS.out.hq
              .join(ch_pre_polish_tsv)
              .map { meta, hq_bed, tsv -> tuple(meta, hq_bed, tsv) }
          ch_truncations = STRIP_TRUNCATIONS.out.discard
          ch_versions = ch_versions
              .mix(DETACH_DUPLICATES.out.versions)
              .mix(ISOTOOLS_TRUNCATION_DETECTOR.out.versions)
              .mix(STRIP_TRUNCATIONS.out.versions)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          COUNTS
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      PREDICT_ORFS.out.counts
      .map { meta, initial, transaid, ns_td, tai, blast, samba, all, unique, kept ->
          def line = "${meta.id}@${meta.chr}\t${initial}\t${transaid}\t${ns_td}\t${tai}\t${blast}\t${samba}\t${all}\t${unique}\t${kept}"
          return line
      }
      .collectFile(
        name: 'counts.tsv',
        newLine: true,
        storeDir: "${output_dir}/XORF_PIPELINE_INFO/XORF_COUNTS"
      )
      .set { ch_counts }

      ch_versions = ch_versions.mix(CONCAT.out.versions)
      ch_versions = ch_versions.mix(database_versions)
      ch_pipeline_versions = ch_versions
          .collectFile(
              name:      "xorf.versions.yml",
              storeDir:  "${output_dir}/XORF_PIPELINE_INFO",
              sort:      true,
              keepHeader: false,
              newLine:   true
          )

      ch_versions = ch_versions.mix(CHUNKER.out.versions)
      ch_versions = ch_versions.mix(GET_CANDIDATES.out.versions)
      ch_versions = ch_versions.mix(PREDICT_ORFS.out.versions)

    emit:
      files        = ch_full_length_predictions // [ meta, bed, tsv ] — stripped HQ when do_polishing
      truncations  = ch_truncations             // [ meta, bed ] STRIP discard; empty if !do_polishing
      duplicates   = ch_duplicates              // [ meta, bed ] DETACH duplicates; empty if !do_polishing
      counts       = ch_counts
      versions     = ch_pipeline_versions
}
