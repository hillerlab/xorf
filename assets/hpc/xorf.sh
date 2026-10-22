#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# SLURM submission script for xorf
#
# Each array task runs one soft-masking job. Nextflow itself is the
# "main job" — it submits all compute work as child SLURM jobs and only needs
# a small memory footprint.
#
# MANIFEST FILE FORMAT (species_list)
# ─────────────────────────────────────
# A tab-separated file with one run per line, no header:
#
#   <regions>  <sequence>  <database> <selenocysteine_sites>  <prefix>
#
# Example:
#   regions.bed  hg38.2bit  swissprot_vertebrates.dmnd  selenos.bed  hg38
#
# Paths must be absolute. 
#
# USAGE
# ─────
# Edit the four path variables below, then submit with:
#   sbatch --array=1-<N> xorf.sh
# where <N> is the number of lines in your manifest file.
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=XORF
#SBATCH --array=1-10        # set upper bound to number of lines in species_list
#SBATCH -t 2-0
#SBATCH --output=/path/to/logs/%A.%a.out  # MODIFY THIS!
#SBATCH --error=/path/to/logs/%A.%a.err   # MODIFY THIS!
#SBATCH --mem=20G          # memory for the Nextflow process itself (not compute jobs)
#SBATCH -p long            # partition name
#SBATCH -q long            # queue name

# ── Load required modules (adjust to your cluster's module system) ────────────
module load nextflow
module load openjdk

# ── Environment ───────────────────────────────────────────────────────────────
export SLURM_SKIP_EPILOG=1

# Directory where Apptainer caches pulled container images
export NXF_APPTAINER_CACHEDIR=/scratch/$USER/xorf/apptainer

# Optional: pre-build a named SIF to avoid the auto-derived cache filename.
# Build once with:
#   apptainer build $NXF_APPTAINER_CACHEDIR/xorf.sif \
#       ghcr.io/hillerlab/xorf:latest
# Then uncomment:
# export NXF_CONTAINER_IMAGE=$NXF_APPTAINER_CACHEDIR/xorf.sif

# Give Nextflow's JVM enough heap for large runs (thousands of jobs)
export NXF_OPTS="-Xms4g -Xmx16g"

# ── Paths — edit these ────────────────────────────────────────────────────────
species_list="/path/to/manifest.tsv"   # tab-separated manifest (see format above)
working_dir="/path/to/output"          # one subdirectory per genome will be created here
pipeline_dir="/path/to/xorf"       # cloned pipeline repo

# ── Parse manifest line for this array task ───────────────────────────────────
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$species_list")
regions=$(echo "$row" | cut -f1)
sequence=$(echo "$row"   | cut -f2)
database=$(echo "$row" | cut -f3)
selenocysteine_sites=$(echo "$row" | cut -f4)
prefix=$(echo "$row" | cut -f5)


if [[ -z "$regions" || -z "$sequence" ]]; then
    echo "ERROR: could not parse line ${SLURM_ARRAY_TASK_ID} of ${species_list}" >&2
    exit 1
fi

# ── Per-pair working directory ─────────────────────────────────────────────────
run_dir="${working_dir}/${prefix}_XORF"
mkdir -p "${run_dir}/logs"

# ── Write params.json for this pair ───────────────────────────────────────────
# Scientific parameters go here; infrastructure stays in nextflow.config.
cat > "${run_dir}/params.json" <<EOF
{
    "//1": "── Input / Output ──────────────────────────────────────────────────────",
    "regions":  "${regions}",
    "sequence": "${sequence}",
    "database": "${database}",
    "outdir": "${run_dir}/results",

    "//2": "── Selenocysteine masking [optional] ──────────────────────────────────",
    "selenocysteine_sites": "${selenocysteine_sites}",

    "//3": "── ORF prediction [optional] ───────────────────────────────────────────",
    "chunk_size": 10,
    "predict_keep_raw": true,
    "predict_min_score_max_predictions": -1,
    "predict_max_predictions": 3,
    "predict_threshold": -1

}
EOF

# cd into run_dir so each run's .nextflow.log is saved there
cd "$run_dir"

nextflow run "${pipeline_dir}/main.nf" \
    -params-file "${run_dir}/params.json" \
    -profile     apptainer,slurm \
    -w           "${run_dir}/work"
