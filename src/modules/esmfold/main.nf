// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOWNLOAD_ESMFOLD_WEIGHTS — Prefetch ESMFold2-Fast + ESMC-6B into a HF hub cache
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DOWNLOAD_ESMFOLD_WEIGHTS {
    label 'process_single'

    conda "${moduleDir}/../blast/environment.yml"
    // INFO: the BLAST image, pinned to the same tag as BLAST (src/modules/blast/main.nf)
    container 'ghcr.io/hillerlab/orf-blast:v0.2.1'

    input:
    path helper

    output:
    path "esmfold_cache", emit: cache
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export HF_HUB_DISABLE_XET=1
    python $helper --prefetch esmfold_cache

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        huggingface_hub: \$(python -c "import huggingface_hub; print(huggingface_hub.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir esmfold_cache

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        huggingface_hub: \$(python -c "import huggingface_hub; print(huggingface_hub.__version__)")
    END_VERSIONS
    """
}
