// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLAST — Validates ORFs against a reference database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BLAST {
    tag "$meta.id:$meta.name"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/hillerlab/orf-blast:latest'

    input:
    tuple val(meta), path(bed), path(sequence), path(predictions), path(net)
    each path(database)
    each path(esm_cache)

    output:
    tuple val(meta), path(bed), path("*/*result"), optional: true, emit: blast
    tuple val(meta), env(INITIAL_REGION_COUNT), env(TRANSLATION_COUNT), optional: true, emit: counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def orf_min_len = task.ext.orf_min_len ?: 50
    def orf_min_percent = task.ext.orf_min_percent ?: 0.25
    def upstream = meta.upstream == null ? task.ext.upstream ?: 1000 : meta.upstream
    def downstream = meta.downstream == null ? task.ext.downstream ?: 1000 : meta.downstream
    def keep_temp = task.ext.keep_temp ? "--keep-temp" : ""

    def prefix = task.ext.prefix ?: "${meta.id}_${meta.name}"
    prefix = prefix.replaceAll('\\.', '_')
    def engine = params.engine ?: 'diamond'
    def db = database.isDirectory() ? "${database}/protein" : "${database}"
    def esm_env = params.esm ? """
    export HUGGINGFACE_HUB_CACHE="${esm_cache}"
    export HF_HUB_CACHE="${esm_cache}"
    export HF_HUB_OFFLINE=1
    export TRANSFORMERS_OFFLINE=1
    export HF_HUB_DISABLE_XET=1
    """ : ""
    """
    ${esm_env}
    orf --threads ${task.cpus} blast \\
    --fasta $sequence \\
    --bed $bed \\
    --tai $predictions \\
    --net $net \\
    --outdir . \\
    --orf-min-len $orf_min_len \\
    --orf-min-percent $orf_min_percent \\
    --database $db \\
    --upstream-flank $upstream \\
    --downstream-flank $downstream \\
    --prefix $prefix \\
    --engine $engine \\
    $keep_temp \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-blast: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        diamond: \$(diamond version 2>&1 | sed 's/^.*diamond version //' )
        mmseqs: \$(mmseqs version 2>&1 | sed 's/^.*MMseqs Version: //; s/ .*\$//' | head -n 1)
        psauron: \$(psauron --version 2>&1 | sed 's/^.*psauron //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.name}"
    prefix = prefix.replaceAll('\\.', '_')
    """
    touch ${prefix}.orf
    touch ${prefix}.orf/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-blast: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        diamond: \$(diamond version 2>&1 | sed 's/^.*diamond version //' )
        mmseqs: \$(mmseqs version 2>&1 | sed 's/^.*MMseqs Version: //; s/ .*\$//' | head -n 1)
        psauron: \$(psauron --version 2>&1 | sed 's/^.*psauron //; s/ .*\$//')
    END_VERSIONS
    """
}
