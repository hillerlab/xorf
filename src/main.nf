#!/usr/bin/env nextflow

/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xorf

    end-to-end robust and comprehensive ORF prediction pipeline
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller

    GitHub:  https://github.com/hillerlab/xorf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info """
    xorf v${workflow.manifest.version}
    End-to-end robust and comprehensive ORF prediction pipeline

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

    Usage (full run):
        mamba activate nextflow
        nextflow run main.nf \\
            --regions        /path/to/regions.{bed, gtf, gff} \\
            --sequence       /path/to/genome.{fa, 2bit, fa.gz} \\
            --outdir         results/ \\
            -profile         apptainer,slurm

    Pass all parameters from a JSON file (replaces old --params_from_file):
        nextflow run main.nf -params-file my_params.json

    Required parameters (full run + fill/clean):
        --regions            PATH    Path to genomic regions (BED, GTF, or GFF)
        --sequence           PATH    Path to genome sequence (FASTA or 2bit)

    Optional parameters (common):
        --outdir              PATH    Output directory [default: ./results]
        --engine              STRING  BLAST engine: diamond or mmseqs2 [default: diamond]
        --esm                 BOOL    Add CPU ESMFold2-Fast mean pLDDT scores [default: false]
        --esmfold_local_weights PATH  Local Hugging Face hub cache; skips DOWNLOAD_ESMFOLD_WEIGHTS [default: null]
        --custom_database     PATH    Path to custom database [default: null]:
                                      diamond: .dmnd/.dmnd.gz replaces the default database;
                                      fasta (.fa/.fasta/.fa.gz/.fasta.gz) is appended to raw_database
                                      mmseqs2: fasta is appended to raw_database; .dmnd is rejected
        --raw_database        URL     Raw protein FASTA (SwissProt). Used for fasta custom_database
                                      and as the default when --engine mmseqs2 [default: UniProt/SwissProt]
        --chunk_size          INT     Chunk size for parallel processing [default: 20]
        --predict_keep_raw    BOOL      Keep raw predictions [default: false]
        --selenocysteine_sites PATH      Selenocysteine masking [default: null]
        --predict_min_score_max_predictions FLOAT   Minimum score for ORF predictions; negative = the model's cutoff (0.29; 0.33 with skip_transaid) [default: -1]
        --predict_max_predictions INT  Maximum number of ORF predictions [default: 3]
        --predict_threshold   FLOAT     Score cutoff for the predicted_class column; negative = the model's cutoff [default: -1]
        --skip_netstart       BOOL      Skip NetStart2; the PREDICT model does not use its scores [default: true]
        --skip_transaid       BOOL      Skip TransAID; PREDICT scores with its model trained without TransAID features [default: false]
        --transaid_tis_cutoff FLOAT     Min TIS probability for TransAID to pair an ATG [default: 0.1]
        --transaid_tts_cutoff FLOAT     Min TTS probability for TransAID to pair a stop [default: 0.1]
        --strip_truncations   BOOL      Drop predictions flagged as 5' truncations; false for annotation-grade input [default: true]
        --run_only_on         BOOL      Run only on masked or unmasked sequences [default: false]
        --run_only_mode       STRING    Run only on masked or unmasked sequences [default: null; options: mask, unmask]
        --run_only_target     STRING    Run only on masked or unmasked sequences [default: null; options: intersect, exclude]

    Profiles:
        local       Run on local machine (default)
        slurm       Submit jobs to SLURM cluster
        conda       Use conda environments
        apptainer   Use Apptainer containers
        singularity Use Singularity containers
        docker      Use Docker containers
        test        Run with bundled test data

    Use --help to show this message.
    """.stripIndent()
    System.exit(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EMAIL }        from './modules/email/main.nf'
include { XORF as MAIN }         from './subworkflows/xorf/main.nf'
include { WGET as WGET_SAMBA_WEIGHTS } from './modules/wget/main.nf'
include { DOWNLOAD_ESMFOLD_WEIGHTS } from './modules/esmfold/main.nf'
include { UNTAR } from './modules/untar/main.nf'
include { WGET as WGET_PROTEIN_DATABASE } from './modules/wget/main.nf'
include { WGET as WGET_MMSEQS_FASTA } from './modules/wget/main.nf'
include { GUNZIP as GUNZIP_DATABASE } from './modules/gunzip/main.nf'
include { FASTA_MERGE } from './modules/diamond/merge/main.nf'
include { FASTA_MERGE as FASTA_MERGE_MMSEQS } from './modules/diamond/merge/main.nf'
include { DIAMOND_MAKEDB } from './modules/diamond/makedb/main.nf'
include { MMSEQS_CREATEDB } from './modules/mmseqs/createdb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateRun() {
    def errors = []
    if (!params.regions)   errors << "  --regions is required"
    if (!params.sequence)  errors << "  --sequence is required"
    if (!(params.engine in ['diamond', 'mmseqs2'])) {
        errors << "  --engine must be 'diamond' or 'mmseqs2' (got: ${params.engine})"
    }
    if (params.engine == 'mmseqs2' && params.custom_database && isDiamondDb(params.custom_database)) {
        errors << "  --engine mmseqs2 cannot use a diamond .dmnd database; pass FASTA or omit custom_database"
    }
    if (params.esm && params.esmfold_local_weights) {
        def weights = file(params.esmfold_local_weights)
        if (!weights.exists()) {
            errors << "  --esmfold_local_weights does not exist: ${params.esmfold_local_weights}"
        } else if (!weights.isDirectory()) {
            errors << "  --esmfold_local_weights must be a Hugging Face hub cache directory"
        } else if (!file("${weights}/models--biohub--ESMFold2-Fast").exists() || !file("${weights}/models--biohub--ESMC-6B").exists()) {
            errors << "  --esmfold_local_weights is not a Hugging Face hub cache (need models--biohub--ESMFold2-Fast and models--biohub--ESMC-6B); pass the published esmfold_cache directory or \$HF_HOME/hub"
        }
    }

    if (errors) {
        log.error "Parameter validation failed:\n${errors.join('\n')}"
        System.exit(1)
    }

    if (params.skip_transaid) {
        log.warn "skip_transaid: TransAID is not run; PREDICT scores with its model trained without TransAID features (0.06-0.30 points fewer correct CDSs on five benchmark species)"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
      email
      email_on_fail
      plaintext_email
      outdir
      use_mailx
      files
      counts
      ch_versions

    main:

      if (params.sent_email) {
          EMAIL (
              email,
              email_on_fail,
              plaintext_email,
              outdir,
              use_mailx,
              files,
              counts,
              ch_versions
          )
      }

      workflow.onError {
          log.error "ERROR: Pipeline failed!"
          log.error "ERROR: Please check the following error message:\n${workflow.errorMessage}"
          log.error "ERROR: Refer to github issues: https://github.com/alejandrogzi/HLorf/issues"
      }

      workflow.onComplete {
          log.info "\nPipeline completed successfully!"
      }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    XORF()
}

workflow XORF {
    
    log.info """
    > xorf v${workflow.manifest.version}
    > End-to-end robust and comprehensive ORF prediction pipeline
    > The Hiller Lab at the Senckenberg Research Institute

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

      Regions:   ${params.regions}
      Sequence:  ${params.sequence}
      Database:  ${params.database} (custom: ${params.custom_database}, engine: ${params.engine})
      Outdir:    ${params.outdir}
      Run only:  ${params.run_only_on} (mode: ${params.run_only_mode}, target: ${params.run_only_target})
      Profile:   ${workflow.profile}
    """.stripIndent()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    validateRun()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SAMBA WEIGHTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_samba_weights = Channel.empty()
    if (params.samba_local_weights) {
        ch_samba_weights = Channel.value(
          file(params.samba_local_weights, checkIfExists: true)
        ).map { path -> [ [id : path.baseName ], path ] }
    } else {
        WGET_SAMBA_WEIGHTS(
          Channel.value(
            params.samba_weights
          ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
        )
        ch_samba_weights = WGET_SAMBA_WEIGHTS.out.outfile
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ESMFOLD WEIGHTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_esm_weights = Channel.empty()
    ch_esm_versions = Channel.empty()
    if (params.esm) {
        if (params.esmfold_local_weights) {
            ch_esm_weights = Channel.value(
                file(params.esmfold_local_weights, checkIfExists: true)
            )
        } else {
            DOWNLOAD_ESMFOLD_WEIGHTS(
                file("${projectDir}/../modules/orf/src/esmfold2_fast.py", checkIfExists: true)
            )
            ch_esm_weights = DOWNLOAD_ESMFOLD_WEIGHTS.out.cache
            ch_esm_versions = DOWNLOAD_ESMFOLD_WEIGHTS.out.versions
        }
    } else {
        ch_esm_weights = Channel.value(file("${projectDir}/params.json"))
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PROTEIN DATABASE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_database = Channel.empty()
    ch_database_versions = Channel.empty()
    if (params.engine == 'mmseqs2') {
      // INFO: mmseqs2 never uses the zenodo diamond .dmnd artifact.
      // Default target is params.raw_database (SwissProt FASTA).
      if (params.custom_database && !isProteinFasta(params.custom_database)) {
          error """
          ERROR: custom_database extension not recognized for --engine mmseqs2.
          Please provide a FASTA file:
            .fa/.fasta/.fa.gz/.fasta.gz -> appended to raw_database and indexed with mmseqs
          """.stripIndent()
      }

      WGET_MMSEQS_FASTA(
          Channel.value(params.raw_database)
          .map { it -> [ [ id: 'uniprot_sprot.fasta.gz' ], it ] }
      )

      ch_mmseqs_fasta = WGET_MMSEQS_FASTA.out.outfile
          .map { meta, fasta -> [ [ id: 'mmseqs_db' ], fasta ] }
      ch_database_versions = Channel.empty()

      if (params.custom_database) {
          FASTA_MERGE_MMSEQS(
              Channel.fromPath(params.custom_database, checkIfExists: true)
                .map { it -> [ [ id: it.baseName ], it ] }
                .combine(WGET_MMSEQS_FASTA.out.outfile.map { meta, raw -> raw })
          )
          ch_mmseqs_fasta = FASTA_MERGE_MMSEQS.out.fasta
              .map { meta, fasta -> [ [ id: 'merged_database' ], fasta ] }
          ch_database_versions = FASTA_MERGE_MMSEQS.out.versions
      }

      MMSEQS_CREATEDB(ch_mmseqs_fasta)
      ch_database = MMSEQS_CREATEDB.out.db.map { meta, it -> it }
      ch_database_versions = ch_database_versions.mix(MMSEQS_CREATEDB.out.versions)
    } else if (params.custom_database) {
      // INFO: .dmnd / .dmnd.gz -> replaces the default database entirely
      if (params.custom_database.endsWith('.dmnd')) {
          ch_database = Channel.fromPath(params.custom_database, checkIfExists: true)
      } else if (params.custom_database.endsWith('.dmnd.gz')) {
          GUNZIP_DATABASE(
              Channel.value(
                  [ [id: params.custom_database.tokenize('/')[-1]], params.custom_database ]
              )
          )
          GUNZIP_DATABASE.out.gunzip
            .map { meta, it -> it }
            .set { ch_database }
      // INFO: .fa / .fasta / .fa.gz / .fasta.gz -> downloaded raw database + custom
      // sequences are merged, then reindexed with DIAMOND_MAKEDB
      } else if (isProteinFasta(params.custom_database)) {
          WGET_PROTEIN_DATABASE(
              Channel.value(params.raw_database)
              .map { it -> [ [ id: 'uniprot_sprot.fasta.gz' ], it ] }
          )
          FASTA_MERGE(
              Channel.fromPath(params.custom_database, checkIfExists: true)
                .map { it -> [ [ id: it.baseName ], it ] }
                .combine(WGET_PROTEIN_DATABASE.out.outfile.map { meta, raw -> raw })
          )
          DIAMOND_MAKEDB(
              FASTA_MERGE.out.fasta.map { meta, fasta -> [ [ id: 'merged_database' ], fasta ] },
              [],
              [],
              []
          )
          DIAMOND_MAKEDB.out.db
            .map { meta, it -> it }
            .set { ch_database }
          ch_database_versions = FASTA_MERGE.out.versions.mix(DIAMOND_MAKEDB.out.versions)
      } else {
          error """
          ERROR: custom_database extension not recognized.
          Please provide a custom database in one of these formats:
            .dmnd/.dmnd.gz   -> replaces the default database entirely
            .fa/.fasta/.fa.gz/.fasta.gz -> appended to the default SwissProt database and reindexed
          """.stripIndent()
      }
    } else {
      WGET_PROTEIN_DATABASE(
        Channel.value(params.database)
        .map { it -> [ [ id: 'uniprot_sprot.tar.gz' ], it ] }
      )
      UNTAR(WGET_PROTEIN_DATABASE.out.outfile)
      ch_database = UNTAR.out.contents.map { meta, it -> it }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MAIN WORKFLOW
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // INFO: chr names the chunks of every downstream task; a random value here gave every
    // INFO: task a new hash on each launch, so -resume could not reuse any work
    MAIN (
       Channel.fromPath(params.regions).map { it -> [ [id: it.baseName, chr: stableHash(it)], it ] },
       Channel.fromPath(params.sequence),
       ch_database,
       params.outdir,
       params.chunk_size,
       ch_samba_weights,
       ch_esm_weights,
       params.predict_keep_raw,
       params.selenocysteine_sites,
       params.skip_netstart,
       params.skip_transaid,
       params.rename_deactivate,
       params.do_polishing,
       params.strip_truncations,
       params.skip_joined_concat,
       params.run_only_on,
       params.run_only_mode,
       params.run_only_target,
       ch_database_versions.mix(ch_esm_versions)
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PIPELINE COMPLETION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.use_mailx,
        MAIN.out.files,
        MAIN.out.counts,
        MAIN.out.versions
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// INFO: 6 hex characters derived from the absolute input path -> same value on every launch
def stableHash(path) {
    def digest = java.security.MessageDigest.getInstance('MD5').digest(path.toAbsolutePath().toString().bytes)
    return digest.encodeHex().toString().take(6)
}

def isProteinFasta(path) {
    def p = path.toString().toLowerCase()
    return p.endsWith('.fa') || p.endsWith('.fasta') || p.endsWith('.fa.gz') || p.endsWith('.fasta.gz')
}

def isDiamondDb(path) {
    def p = path.toString().toLowerCase()
    return p.endsWith('.dmnd') || p.endsWith('.dmnd.gz')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (workflow.success) {
        def intermediate_dir = new File(params.outdir as String, '02_merged')
        def intermediate_beds = intermediate_dir.exists() ? (intermediate_dir.listFiles()?.findAll { it.name.endsWith('.bed') } ?: []) : []
        
        def results_dir = new File(params.outdir as String, '04_results')
        def final_beds = results_dir.exists() ? (results_dir.listFiles()?.findAll { it.name.endsWith('.bed') } ?: []) : []
        log.info "Pipeline completed successfully!"
        if (final_beds) {
            log.info "Final predictions: ${final_beds.collect { it.toString() }.join(', ')}"
        } else {
            if (intermediate_beds) {
                log.info "Final predictions: ${intermediate_beds.collect { it.toString() }.join(', ')}"
            }

            log.warn "Pipeline reported success but final bed file was not produced - check that all steps ran"
        }
        log.info "Run time   : ${workflow.duration}"
    } else {
        log.error "Pipeline FAILED — ${workflow.errorMessage}"
    }
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
