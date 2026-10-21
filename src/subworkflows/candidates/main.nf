// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    candidates

    Runs translation, RNASamba, NetStart, TransAID, BLAST to get ORF candidates
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller

    GitHub:  https://github.com/hillerlab/xorf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TRANSLATION } from '../../modules/translation/main.nf'
include { RNASAMBA }    from '../../modules/rnasamba/main.nf'
include { NETSTART }    from '../../modules/netstart/main.nf'
include { TRANSAID }    from '../../modules/transaid/main.nf'
include { BLAST }       from '../../modules/blast/main.nf'
include { JOIN as JOIN_NETS }   from '../../modules/join/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GET_CANDIDATES {
    take:
      ch_pairs       // channel: [ meta, path ]
      ch_database    // channel: [ meta, path ]
      samba_weights  // path
      esm_weights    // path
      skip_netstart  // boolean
      skip_transaid  // boolean

    main:
      ch_versions = Channel.empty()

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          TRANSLATION
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      TRANSLATION(ch_pairs)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          RNASAMBA & TRANSAID
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      RNASAMBA(
        ch_pairs,
        samba_weights
      )

      // INFO: skip_transaid -> a PREDICT model trained without TransAID features needs no
      // INFO: TransAID run; a header-only CSV makes the blast step write the missing value
      // INFO: (-1) to every TransAID column, as for candidates TransAID did not score
      if (skip_transaid) {
        ch_transaid_csv = RNASAMBA.out.bed.map { meta, bed ->
          [ meta, file("${projectDir}/../assets/transaid/empty.transaid.csv", checkIfExists: true) ]
        }
        ch_transaid_count = RNASAMBA.out.bed.map { meta, bed -> [ meta, 0 ] }
      } else {
        TRANSAID(
          RNASAMBA.out.fasta,
          RNASAMBA.out.bed
        )
        ch_transaid_csv = TRANSAID.out.transaid
        ch_transaid_count = TRANSAID.out.count
      }

      if (!skip_netstart) {
        NETSTART(
          RNASAMBA.out.fasta,
          RNASAMBA.out.bed
        )

        ch_netstart = NETSTART.out.netstart

        ch_transaid_csv
          .join(RNASAMBA.out.bed)
          .join(ch_netstart)
          .set { ch_nets }

        ch_nets.map { meta, transaid, bed, netstart -> [ meta, transaid, bed ] }.set { ch_transaid }
        ch_nets.map { meta, transaid, bed, netstart -> [ meta, netstart ] }.set { ch_netstart }

        JOIN_NETS( ch_transaid, ch_netstart, skip_netstart )
        ch_versions = ch_versions.mix(NETSTART.out.versions)
      } else {
        ch_transaid_csv
          .join(RNASAMBA.out.bed)
          .set { ch_nets }

        JOIN_NETS( ch_nets, Channel.value([[:], []]), skip_netstart )
      }

      TRANSLATION.out.predictions
      .join(JOIN_NETS.out.net)
      .set { ch_pre_candidates }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          BLAST
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: BLAST.out.blast is a channel of [meta, bed, tsv]
      // INFO: RNASAMBA.out.samba is a channel of [meta, tsv]
      BLAST(
          ch_pre_candidates,
          ch_database,
          esm_weights
      )
      .blast
      .join(RNASAMBA.out.samba)
      .set { ch_candidates }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          COUNTS
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      BLAST.out.counts
      .join(ch_transaid_count)
      .join(JOIN_NETS.out.count)
      .set { ch_counts }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VERSIONS
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_versions = ch_versions.mix(TRANSLATION.out.versions)
      ch_versions = ch_versions.mix(RNASAMBA.out.versions)
      ch_versions = ch_versions.mix(BLAST.out.versions)

    emit:
      candidates = ch_candidates
      counts = ch_counts
      versions = ch_versions
}
