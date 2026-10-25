// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRANSLATION — Translates predicted ORFs and validates translation products
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TRANSLATION {
    tag "$meta.id:$meta.name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // INFO: pinned, not latest: a cached older latest image keeps an orf binary that panics on a CDS
    // INFO: at the first base of a scaffold (fixed in orf 0.0.33). Bump the tag with every release
    // INFO: that changes modules/orf.
    container 'ghcr.io/hillerlab/orf-tai:v0.2.1'

    input:
    tuple val(meta), path(bed), path(sequence)

    output:
    tuple val(meta), path(bed), path(sequence), path("${meta.id}_${meta.name}/*result"), optional: true, emit: predictions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def upstream = meta.upstream == null ? task.ext.upstream ?: 1000 : meta.upstream
    def downstream = meta.downstream == null ? task.ext.downstream ?: 1000 : meta.downstream
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.name}"
    """
    orf tai \\
    --fasta $sequence \\
    --bed $bed \\
    --outdir ${prefix} \\
    -u $upstream \\
    -d $downstream
    
    mv ${prefix}/tai/*result ${prefix}/ && rm -rf ${prefix}/tai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-tai: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        translationai: 0.0.1
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}
    touch ${prefix}/tai
    touch ${prefix}/tai/*result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-tai: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        translationai: 0.0.1
    END_VERSIONS
    """
}
