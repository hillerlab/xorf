// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RNASAMBA — Classifies ORFs as coding or non-coding using machine learning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RNASAMBA {
    tag "$meta.id:$meta.name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/hillerlab/orf-samba:latest'

    input:
    tuple val(meta), path(bed), path(sequence)
    tuple val(_meta), path(weights)

    output:
    tuple val(meta), path("${prefix}/samba/*tsv")      , optional: true, emit: samba
    tuple val(meta), path("${meta.name}.tmp.strip.fa") , optional: true, emit: fasta
    tuple val(meta), path(bed)                        , optional: true, emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def upstream = meta.upstream == null ? task.ext.upstream ?: 1000 : meta.upstream
    def downstream = meta.downstream == null ? task.ext.downstream ?: 1000 : meta.downstream
    prefix = task.ext.prefix ?: "${meta.id}_${meta.name}"
    """
    orf samba \\
    --fasta $sequence \\
    --outdir ${prefix} \\
    --upstream-flank $upstream \\
    --downstream-flank $downstream \\
    --weights $weights \\
    $args

    // INFO: orf writes the stripped FASTA next to its input and the TSV under
    // INFO: ${prefix}/samba. Declare both where they are written instead of
    // INFO: renaming them: rename(2) can return EBUSY on shared filesystems.
    [ -f ${meta.name}.tmp.strip.fa ] || cp $sequence ${meta.name}.tmp.strip.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-samba: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        rnasamba: \$(rnasamba --version 2>&1 | tail -n 1 | sed 's/^rnasamba //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.name}"
    """
    touch ${meta.name}.tmp.strip.fa
    mkdir -p ${prefix}/samba
    touch ${prefix}/samba/${meta.name}.tmp.samba.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-samba: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        rnasamba: \$(rnasamba --version 2>&1 | tail -n 1 | sed 's/^rnasamba //')
    END_VERSIONS
    """
}
