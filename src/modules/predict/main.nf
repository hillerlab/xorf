// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREDICT — Integrates BLAST and RNASamba predictions to score and filter ORFs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PREDICT {
    tag "$meta.id:$meta.name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // INFO: pinned, not latest: the model, its cutoff (-1 = bundle cutoff) and --no-transaid change
    // INFO: together with this module, and a cached older latest image would read -1 as "keep every
    // INFO: candidate". Bump the tag with every release that changes modules/predict.
    container 'ghcr.io/hillerlab/orf-predict:v0.2.0'

    input:
    tuple val(meta), path(bed), path(blast), path(samba)

    output:
    tuple val(meta), path("${meta.id}_${meta.name}/${meta.id}_${meta.name}.predictions.bed"), path("${meta.id}_${meta.name}/${meta.id}_${meta.name}.predictions.tsv"), emit: orfs
    tuple val(meta), path("${meta.id}_${meta.name}/*all.predictions.bed"), path("${meta.id}_${meta.name}/*all.predictions.tsv"), emit: raw, optional: true
    tuple val(meta), env(BLAST_PREDICTION_COUNT), env(SAMBA_PREDICTION_COUNT), env(ALL_PREDICTED_REGIONS), env(UNIQUE_PREDICTED_REGIONS), env(KEPT_REGIONS), emit: counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // INFO: negative cutoffs mean the cutoff stored in the model bundle
    def threshold = task.ext.threshold != null ? task.ext.threshold : -1
    def min_score_max_predictions = task.ext.min_score_max_predictions != null ? task.ext.min_score_max_predictions : -1
    def max_predictions = task.ext.max_predictions ?: 3
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.name}"
    """
    # INFO: predict.py < 0.1.0 (v2 model) reads the -1 cutoffs as "keep every candidate"
    predict.py --version 2>&1 | grep -Eq 'Version: (0\\.[1-9]|[1-9])' || {
        echo "ERROR: PREDICT needs predict.py >= 0.1.0 (orf-predict >= v0.2.0); found \$(predict.py --version 2>&1)" >&2
        exit 1
    }

    predict.py \\
    --blast $blast \\
    --samba $samba \\
    --alignments $bed \\
    --outdir ${prefix} \\
    --prefix ${meta.id}_${meta.name} \\
    --threshold $threshold \\
    --min-score-max-predictions $min_score_max_predictions \\
    --max-predictions $max_predictions \\
    $args

    BLAST_PREDICTION_COUNT=\$(wc -l < ${blast})
    SAMBA_PREDICTION_COUNT=\$(wc -l < ${samba})

    ALL_PREDICTED_REGIONS=\$(wc -l ${prefix}/${prefix}.predictions.tsv | awk '{print \$1}')
    UNIQUE_PREDICTED_REGIONS=\$(awk '{print \$4}' ${prefix}/${prefix}.predictions.tsv | sed -E 's/(_ORF|\\.p[0-9]).*//' | sort | uniq | wc -l)
    KEPT_REGIONS=\$(wc -l ${prefix}/${prefix}.predictions.bed | awk '{print \$1}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        predict.py: \$(predict.py --version 2>&1 | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.name}"
    """
    touch ${prefix}
    touch ${prefix}/*bed
    touch ${prefix}/*tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        predict.py: \$(predict.py --version 2>&1 | sed 's/.*Version: //')
    END_VERSIONS
    """
}
