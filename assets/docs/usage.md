<p align="center">
  <p align="center">
    <img width=100 align="center" src="../figures/xorf.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="../figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="../figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="../figures/hillerlab-light.png"
    >
  </picture>
</p>


  <span>
    <h1 align="center">
        xorf
    </h1>
  </span>

  <span>
    <h2 align="center">
        USER GUIDE
    </h2>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/xorf" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/xorf?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://www.genome.gov/genetics-glossary/Open-Reading-Frame">orf</a> .
        <a href="https://github.com/alejandrogzi/xorf/blob/master/assets/pipeline/xorf.mermaid">pipeline</a> .
        <a href="https://github.com/alejandrogzi/xorf/blob/master/assets/docs/usage.md">usage</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>


This guide explains how to run xorf, what every parameter does, and what you get
back in the results folder. If something is confusing,
let us know in a [GitHub issue](https://github.com/alejandrogzi/xorf/issues).

> [!IMPORTANT]
> - **What is xorf?** An end-to-end pipeline that predicts open reading frames (ORFs)
>  in a genome. You give it a genome and a list of genomic regions; it returns the
>  most confident ORF predictions as BED + TSV files.
> - **How does it work?** Each region is translated and scored with several
>  independent methods — translationAI, RNASamba, PSAURON, TransAID, and
>  BLAST/DIAMOND against a protein database. A machine-learning classifier
>  (XGBoost) combines these signals, and only the ORFs that pass its confidence
>  threshold make it to the final output. NetStart2 can also be run, but the
>  classifier does not use its scores.

---

## Before you start

You need:

| Requirement | Notes |
|---|---|
| Nextflow ≥ 25.04.6 | `nextflow -version` to check |
| Java (JDK 21+) | required by Nextflow |
| Docker **or** Apptainer **or** Singularity | used to run the tools in containers |
| The repository | `git clone https://github.com/hillerlab/xorf.git && cd xorf` |

In the provided environment, make Nextflow available first with
`mamba activate nextflow`.

### Input files

| File | What it is | Accepted formats | Required? |
|---|---|---|---|
| Regions | The genomic regions you want to scan for ORFs | BED, GTF, GFF | **Yes** |
| Genome sequence | The genome those regions come from | FASTA (`.fa`, `.fa.gz`), 2bit | **Yes** |
| Selenocysteine sites | Known TGA codons that code for selenocysteine (see [Selenocysteine masking](#selenocysteine-masking)) | BED | No |
| Protein database | A protein database used by BLAST (DIAMOND `.dmnd`, or FASTA for mmseqs2 / custom sequences) | `.dmnd` / `.dmnd.gz`, or FASTA (`.fa`, `.fasta`, optionally gzipped) | No |

> [!NOTE]
> You do **not** need to provide a protein database. If `custom_database` is left
> empty, xorf automatically downloads UniProt/SwissProt when the pipeline starts
> (the zenodo DIAMOND archive, or the FASTA from `raw_database` when `engine` is `mmseqs2`).
---

## Quick start

There are three ways to run the pipeline. Pick whichever suits you.

### 1. Edit `params.json` (recommended)

Copy the template, fill in your two required inputs, and run:

```bash
nextflow run main.nf -params-file params.json -profile docker
# or, if you use Apptainer / Singularity:
nextflow run main.nf -params-file params.json -profile apptainer
```

The bundled `params.json` template already contains sensible defaults for every
optional setting — most of the time you only need to set `regions` and `sequence`:

```json
{
    "regions":  "/path/to/your/regions.bed",
    "sequence": "/path/to/your/genome.2bit",
    "outdir":   "results"
}
```

### 2. Command-line flags

Set parameters directly on the command line (useful for one-off runs):

```bash
nextflow run main.nf \
    --regions  /path/to/regions.bed \
    --sequence /path/to/genome.2bit \
    --outdir   results \
    -profile   docker
```

> [!NOTE]
> Command-line flags override `params.json`. You can mix both — the JSON file
> sets the base values, flags override individual ones.

### 3. Smoke test (bundled test data)

To check that everything works before running your real data:

```bash
nextflow run main.nf -profile test,apptainer
# or with Docker:
nextflow run main.nf -profile test,docker
```

This runs the pipeline on a small test genome. Results land in `test_results/`.

> [!TIP]
> `nextflow run main.nf --help` prints a summary of the most common parameters.

---

## Parameter reference

Parameters are set in `params.json` or with `--name value` flags. Defaults live
in `nextflow.config` — you should rarely need to touch that file; scientific
settings belong in `params.json`.

### Required inputs

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `regions` | path | — | Genomic regions to scan for ORFs (BED, GTF, or GFF). The pipeline **fails** if this is missing. |
| `sequence` | path | — | Genome sequence the regions refer to (FASTA, gzipped FASTA, or 2bit). The pipeline **fails** if this is missing. |

### Input / output

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `outdir` | path | `./results` | Where all results are written. Set this explicitly per run so different runs don't overwrite each other. |

### Protein database

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `engine` | string | `diamond` | BLAST search engine: `diamond` or `mmseqs2`. |
| `esm` | boolean | `false` | Add mean ESMFold2-Fast pLDDT (0–100) to BLAST output. V1 runs locally on CPU in FP32. When true, `DOWNLOAD_ESMFOLD_WEIGHTS` prefetches the pinned Hugging Face hub cache once and every BLAST task reuses it. |
| `esmfold_local_weights` | path | `null` | Path to a local Hugging Face hub cache containing `models--biohub--ESMFold2-Fast` and `models--biohub--ESMC-6B`. Set this with `esm` to skip `DOWNLOAD_ESMFOLD_WEIGHTS`. |
| `database` | URL | Zenodo UniProt/SwissProt (DIAMOND) | URL the pipeline downloads the protein database from. Only used when `engine` is `diamond` and `custom_database` is not set. |
| `raw_database` | URL | UniProt/SwissProt sequences (`uniprot_sprot.fasta.gz`) | Raw protein FASTA. Used as the default when `engine` is `mmseqs2`, and as the sequences a FASTA `custom_database` is appended to. |
| `custom_database` | path | `null` | **Your own** database. Give a DIAMOND database (`.dmnd`, optionally `.dmnd.gz`) to *replace* the default (diamond engine only), or a FASTA file (`.fa`, `.fasta`, optionally gzipped) to *append* to SwissProt — xorf merges the sequences and rebuilds the engine's database on the fly. |

**When to change this:** leave everything alone for a quick start — the default
database is downloaded automatically. Use `engine` `mmseqs2` to search with
MMseqs2 instead of DIAMOND (faster; uses `raw_database`, not the zenodo `.dmnd`).
Use `custom_database` if you want to:

- **Replace** the default with a curated/species-specific DIAMOND database
  (`.dmnd`), e.g. because you are offline or the download is slow on your
  cluster. Not valid with `engine` `mmseqs2`.
- **Extend** the default with your own proteins (FASTA): any sequence you
  provide is appended to the SwissProt set and the combined set is reindexed
  before BLAST runs. Use `raw_database` to control what your sequences are
  appended to (e.g. full UniProt instead of SwissProt).

Any other `custom_database` extension stops the pipeline with a clear error.

ESMFold2-Fast uses the ESMC-6B backbone (~26 GB hub cache, ~26 GB RAM, CPU
FP32). BLAST is a `process_high_memory` task (4 CPUs, 30 GB, 10 h) so the
ordinary `process_low` defaults are left alone. The pipeline downloads the
pinned checkpoints once (`DOWNLOAD_ESMFOLD_WEIGHTS`) and publishes them to
`${outdir}/00_weights/esmfold_cache`. Set `esmfold_local_weights` to that
directory — or to `$HF_HOME/hub` — to skip the download on later runs. Do not
pass a raw `HF_HOME` (the parent of `hub/`).

### RNASamba

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `samba_weights` | URL | RNASamba `full_length_weights.hdf5` | URL the pipeline downloads the RNASamba model weights from. |
| `samba_local_weights` | path | `null` | Path to a local RNASamba weights file. Set this to use your own weights and skip the download. |

**When to change this:** almost never. `samba_local_weights` is handy on offline
clusters or if you trained custom weights.

### NetStart

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `skip_netstart` | boolean | `true` | Skip the NetStart2 step entirely. NetStart2 predicts translation initiation sites with a neural network; when skipped, the net file's NetStart columns are filled with zeros. |

**When to change this:** NetStart2 is **skipped** by default
(`skip_netstart: true`) because it downloads a model at run time and is the
slowest predictor. The PREDICT model (since xorf 0.2.0) does not use NetStart2
scores, so `skip_netstart: false` only adds the NetStart columns to the raw
tables; calls and scores are the same.

### TransAID

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `skip_transaid` | boolean | `false` | Do not run TransAID. Every TransAID column of the candidate table is then missing (the same value a candidate gets when TransAID did not score it), and PREDICT scores with its second bundled model, trained without TransAID features. |
| `transaid_tis_cutoff` | float | `0.1` | Minimum start-site (TIS) probability for TransAID to pair an ATG with a stop codon. |
| `transaid_tts_cutoff` | float | `0.1` | Minimum stop-site (TTS) probability for TransAID to pair a stop codon. |

**When to change this:** TransAID is the slowest step that runs by default. The default model uses
only whether TransAID paired an ATG with a candidate's stop codon (per candidate
family and per transcript), and calls 0.06–0.30 points more CDSs correctly than the
model without TransAID features on five benchmark species. Set `skip_transaid:
true` when runtime matters more; the pipeline warns when the step is skipped.
TransAID only reports (start, stop) pairs whose probabilities pass both cutoffs,
so lower cutoffs export a strict superset of the default rows (filtering them
back to `> 0.1` reproduces the default output exactly). Both bundled models were
trained on the `0.1` export: keep the cutoffs at their defaults unless you
retrain.

### Selenocysteine masking

Selenocysteine (Sec) is an amino acid coded by the stop codon TGA in specific
contexts. If you don't mask these sites, a real selenocysteine TGA can be
misread as a stop codon — truncating otherwise valid ORFs. The selenocysteine
site file must list these regions.

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `selenocysteine_sites` | path | `null` | BED file of selenocysteine sites. When set, the pipeline **masks** those sites in the genome (TGA → A) and also runs a second, unmasked branch on the regions that overlap them. |
| `run_only_on` | boolean | `false` | Run **only** one branch — either the masked or the unmasked one — instead of both. Requires `selenocysteine_sites`. |
| `run_only_mode` | string | `null` | Which branch to keep: `mask` or `unmask`. Only read when `run_only_on` is `true`. |
| `run_only_target` | string | `null` | Which regions to keep: `intersect` (only regions overlapping selenocysteine sites) or `exclude` (only regions *not* overlapping). Only read when `run_only_on` is `true`. |

The four combinations of `run_only_on: true`:

| `run_only_mode` | `run_only_target` | Result |
|---|---|---|
| `mask` | `intersect` | Predict ORFs on the masked genome, only in regions that contain selenocysteine sites |
| `mask` | `exclude` | Predict ORFs on the masked genome, only in regions *without* selenocysteine sites |
| `unmask` | `intersect` | Predict ORFs on the original genome, only in regions containing selenocysteine sites |
| `unmask` | `exclude` | Predict ORFs on the original genome, only in regions *without* selenocysteine sites |

> [!NOTE]
> Setting `run_only_on: true` without `selenocysteine_sites` (or with an invalid
> `run_only_mode` / `run_only_target`) stops the pipeline with a clear error.

**When to change this:** only if you specifically study selenocysteine genes or
want to compare masked vs. unmasked predictions in a controlled way. The default
(`run_only_on: false`) already gives you both views.

### Chunking & prediction

These parameters control how regions are split up for parallel processing and
how the XGBoost classifier filters ORFs.

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `chunk_size` | integer | `20` | Number of regions per chunk. Regions are split into chunks so they can be processed in parallel. Smaller = more, smaller jobs; larger = fewer, bigger jobs. |
| `predict_keep_raw` | boolean | `false` | Keep the **raw** predictions (every ORF the classifier saw, before filtering) in addition to the final set. See [Outputs](#outputs). |
| `predict_min_score_max_predictions` | float | `-1` | Minimum score (0–1) a prediction needs to be kept. Negative = the cutoff stored with the model (`0.29`; `0.33` with `skip_transaid`). |
| `predict_max_predictions` | integer | `3` | Maximum number of predictions kept per query. |
| `predict_threshold` | float | `-1` | Score cutoff for the `predicted_class` column of the TSV (1 = coding). It does not remove predictions. Negative = the cutoff stored with the model. |

**When to change this:**

- `chunk_size`: on a cluster, tune this to match your queue — too small creates
  thousands of tiny jobs, too large creates memory-hungry ones. The bundled
  `params.json` template uses `10`; the default is `20`.
- `predict_min_score_max_predictions`: the model's cutoff is set so that 0.8 %
  of the training lncRNA transcripts (human, mouse and rat, out-of-fold) get a
  call. On five species not used for training (cow, dog, elephant, platypus,
  harpy eagle), it called:
  - 0.4–1.2 % of the lncRNAs;
  - the annotated CDS for 98.4–99.3 % of coding transcripts.

  Raise it for fewer, higher-confidence ORFs; lower it to call more transcripts
  at the cost of more false calls. Scores of earlier xorf versions (model v2) are
  not comparable: parameter files that still set `0.50` call the annotated CDS
  for 0.55–0.85 % fewer coding transcripts with the current model.
- `predict_max_predictions`: leave at the default unless you need more than 3
  candidate ORFs per query.

### Merging & polishing

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `skip_joined_concat` | boolean | `false` | Skip merging all renamed per-region files into one combined file. When `true`, each region keeps its own renamed BED/TSV as the final output (nothing is merged). |
| `do_polishing` | boolean | `true` | Clean up predictions before the final output: remove **duplicate** predictions (marked `#DU`) and strip **5′-truncated** predictions (`iso-utr`: the transcript starts inside a non-first exon of another prediction). Set to `false` to keep everything untouched. |
| `strip_truncations` | boolean | `true` | With `do_polishing`, drop the predictions `iso-utr` flags as truncated. `false` keeps them in `04_results` while still separating `#DU` duplicates. |

**When to change this:** the default (`do_polishing: true`) is what most people
want — it removes likely artefacts. Turn polishing off if you want to do your own
filtering downstream. `skip_joined_concat` is mostly a technical switch; keep the
default unless you need per-region files as the final result.
Set `strip_truncations: false` when the input regions come from a curated
annotation (e.g. RefSeq, GENCODE) rather than from reads: there, alternative
transcription start sites look like truncations, and stripping them removes
real, annotated CDSs.

### Renaming / output naming

After prediction, xorf renames every ORF with a human-readable, information-rich
identifier so you can tell at a glance what each prediction is. A renamed ID
looks like this:

```
ROOT#SC{score}#DU#OR{N}#TI#IN{M}@protein
```

| Part | Meaning |
|---|---|
| `ROOT` | The original region/query name |
| `UNMASK#` prefix | Added automatically for predictions from the unmasked branch |
| `SC{score}` | ORF score (0–100), from `prob_coding` |
| `DU` | This ORF was flagged as a duplicate |
| `OR{N}` | ORF number |
| `TI` | Translation index present |
| `IN{M}` | Inner number |
| `@protein` | Best BLAST reference hit (protein name) |

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `rename_deactivate` | boolean | `false` | Turn renaming off completely: predictions keep their original, opaque identifiers and are copied through unchanged. |
| `rename_rebase` | boolean | `false` | Replace the original root with a short, collision-free content hash (BLAKE2s). Identical predictions across runs then get identical names — handy for cross-run comparison. |
| `rename_append_orf_score` | boolean | `true` | Append the ORF score (`SC{score}`) to each name. |
| `rename_append_protein_name` | boolean | `true` | Append the best BLAST hit (`@protein`) to each name. |
| `rename_custom_prefix` | string | `null` | Extra prefix prepended to every generated name (e.g. your species or run name). |

> [!WARNING]
> The characters `#` and `@` are structural separators in the generated names —
> they are **not allowed** inside `rename_custom_prefix`. The pipeline errors if
> you use them.

**When to change this:** `rename_rebase` is nice when you compare predictions
across pipeline runs (same predictions → same names, no matter the input file
name). `rename_custom_prefix` is a good way to tag a run (e.g. `hg38`).

### Advanced / rarely touched

These are technical switches. You can safely ignore them.

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `publish_dir_mode` | string | `copy` | How results are copied into `outdir` (`copy`, `symlink`, `link`, `move`). |
| `publish_all` | boolean | `false` | Publish every intermediate file, not just the curated outputs. |
| `help` | boolean | `false` | Print the parameter summary and exit. |
| `help_full` | boolean | `false` | Print the full parameter help and exit. |
| `show_hidden` | boolean | `false` | Include hidden parameters in the help output. |
| `version` | boolean | `false` | Print the pipeline version and exit. |
| `trace_report_suffix` | string | timestamp | Suffix for the Nextflow trace/report file names in `pipeline_info/`. |
| `config_profile_name` | string | `null` | Display name for the active profile (informational). |
| `config_profile_description` | string | `null` | Description of the active profile (informational). |

---

## Outputs

Results are written under `outdir`. A typical successful run looks like this:

```
results/
├── 00_concat/                  Concatenated predictions, one file per region (before renaming)
│   └── raw/                    Raw (unfiltered) predictions — only with predict_keep_raw
├── 01_renamed/                 Renamed predictions, one file per region
├── 02_merged/                  Merged renamed predictions (the full-length set)
│   └── raw/                    Merged raw renamed predictions — only with predict_keep_raw
├── 03_duplicates/              Predictions flagged and removed as duplicates
├── 04_results/                 Final polished predictions (*.hq.bed)
├── XORF_PIPELINE_INFO/         Run summary: per-stage counts + software versions
└── pipeline_info/              Nextflow reports: timeline, trace, report, DAG
```

### The output directories

| Directory | Contents | Notes |
|---|---|---|
| `00_concat/` | `*.bed` + `*.tsv` per region, before renaming | Intermediate step — most users never look here |
| `00_concat/raw/` | Raw predictions per region | Only exists with `predict_keep_raw: true` |
| `01_renamed/` | `*.renamed.bed` + `*.renamed.tsv` per region | Skips when `rename_deactivate: true` |
| `02_merged/` | One merged `*.bed` + `*.tsv` combining all renamed predictions | Skips when `skip_joined_concat: true` |
| `02_merged/raw/` | Merged raw predictions | Only with `predict_keep_raw: true` |
| `03_duplicates/` | `*.duplicates.bed` — predictions removed as duplicates | Only with `do_polishing: true` |
| `04_results/` | `*.hq.bed` — the **final predictions**: duplicates removed, 5′-truncated predictions stripped (unless `strip_truncations: false`) | The files you want, and what `XORF.out.files` emits when `do_polishing: true`. Truncation discards are `XORF.out.truncations`. |
| `XORF_PIPELINE_INFO/` | `XORF_COUNTS/counts.tsv` + `xorf.versions.yml` | See below |
| `pipeline_info/` | `execution_timeline_*.html`, `execution_trace_*.txt`, `execution_report_*.html`, `pipeline_dag_*.html` | Nextflow's own run reports |

> [!TIP]
> If `do_polishing` is off, there is no `03_duplicates` / `04_results` — the
> final predictions are the files from `02_merged` (or from `01_renamed` with
> `skip_joined_concat`).

### The BED + TSV pair

Every prediction appears twice: once in a BED12 file (genomic coordinates) and
once in a TSV file (per-prediction metadata). The `id` column links the two —
after renaming it carries the score, ORF number, duplicate flag, and best BLAST
hit (see [Renaming](#renaming--output-naming)). The TSV includes the classifier
score (`prob_coding`) and the BLAST reference (`blast_reference_id`), among
other features.

### counts.tsv — how many ORFs survived each step

`XORF_PIPELINE_INFO/XORF_COUNTS/counts.tsv` holds one row per region, showing
how many predictions survived each pipeline stage:

| Column | What it counts |
|---|---|
| `id` | Region identifier + run hash (e.g. `regions@K3x9bQ`) |
| `initial` | Regions/ORFs at the start of the BLAST step |
| `translation` | ORFs that translated successfully |
| `transaid` | TransAID predictions |
| `net` | Predictions in the merged NetStart + TransAID net file |
| `blast` | BLAST/DIAMOND candidate hits |
| `samba` | RNASamba candidate predictions |
| `all` | All predicted ORFs (before filtering) |
| `unique` | Unique predicted ORFs |
| `kept` | ORFs kept in the final output |

A quick look at this file tells you where candidates are being lost — for
example, `blast` much smaller than `initial` usually means the protein database
is a poor match for your species.

### xorf.versions.yml

Lists the exact version of every tool used in the run (Nextflow, `orf`,
`predict.py`, DIAMOND, etc.) — useful for reproducing a run later.

---

## Profiles

Profiles select *how* the pipeline runs (containers, executor). Pass them with
`-profile name1,name2` (comma-separated, e.g. `-profile apptainer,slurm`).

| Profile | When to use |
|---|---|
| `local` | Run on your own machine (default) |
| `docker` | Run with Docker containers |
| `apptainer` | Run with Apptainer containers (recommended on clusters) |
| `singularity` | Run with Singularity containers |
| `conda` / `mamba` | Use Conda environments instead of containers |
| `slurm` | Submit jobs to a SLURM cluster (combine with `apptainer` or `docker`) |
| `gpu` | Enable GPU support for container engines (used by some predictors) |
| `test` | Run with the bundled test data |
| `debug` | Extra debugging output and validation |
| `arm`, `podman`, `shifter`, `charliecloud`, `wave`, `gitpod` | Specialised environments |

---

## Running on a SLURM cluster

A helper script is provided at `assets/hpc/xorf.sh`. It runs one Nextflow job per
species as a SLURM job array; each Nextflow job then submits all compute work as
child SLURM jobs.

1. Edit the path variables at the top of the script: the manifest path, the
   output directory, and the pipeline directory.
2. Create a tab-separated manifest, one run per line, no header:

   ```
   <regions>  <sequence>  <database>  <selenocysteine_sites>  <prefix>
   ```

   Example:

   ```
   /path/regions.bed  /path/hg38.2bit  /path/swissprot.dmnd  /path/selenos.bed  hg38
   ```

3. Submit with `sbatch --array=1-<N> xorf.sh`, where `<N>` is the number of
   lines in the manifest.

Each array task generates its own `params.json` and results folder. Partition
routing, array sizes, and resource tiers are configured in `nextflow.config` —
edit there to match your cluster.

---

## FAQ & tips

**Q: I only want the masked (or unmasked) predictions — how?**
Set `run_only_on: true` with `run_only_mode: mask` (or `unmask`) and
`run_only_target: intersect` or `exclude`. See [Selenocysteine masking](#selenocysteine-masking).

**Q: How do I see the raw predictions the classifier saw?**
Set `predict_keep_raw: true`. Raw predictions appear under `00_concat/raw/` and
`02_merged/raw/`.

**Q: The pipeline finished but I have no final BED in `04_results/`?**
`04_results` only exists when `do_polishing: true`. If polishing is off, the
final files are in `02_merged/`. If it's on and the folder is still empty, check
`XORF_PIPELINE_INFO/XORF_COUNTS/counts.tsv` — the `kept` column shows whether
any ORF survived filtering.

**Q: Can I tune how many ORFs come out?**
Yes — `predict_min_score_max_predictions` (lower = more ORFs, more false
positives) and `predict_max_predictions` (how many per query). Polishing
(`do_polishing`) can also be turned off to keep everything.

**Q: Does the pipeline download the protein database every run?**
When you don't provide `custom_database`, xorf downloads the default
UniProt/SwissProt database at run start. On shared clusters, point
`custom_database` at a DIAMOND database (`.dmnd`) you already have to save time
and bandwidth. A FASTA `custom_database` still triggers a download (of the raw
sequences) because the merge step needs them.

**Q: My renamed IDs look cryptic. Can I make them nicer?**
Use `rename_custom_prefix` to prepend a label (e.g. `hg38`), and turn off
`rename_append_orf_score` / `rename_append_protein_name` if you want shorter
names. Remember: no `#` or `@` in your prefix.

**Q: Can I run the same input twice and get identical names?**
Yes — set `rename_rebase: true`. Roots are replaced by content hashes, so
identical predictions produce identical identifiers across runs.

**Q: How do I get identical outputs between two machines?**
`XORF_PIPELINE_INFO/xorf.versions.yml` records every tool version. Keep the
container images and pipeline version identical, and results will be
reproducible.

**Q: I upgraded xorf. Do I need new container images?**
Most modules use `:latest` images, and Docker and Apptainer never pull a cached
`:latest` again. After an upgrade, delete the cached `*-latest` images (or point
`apptainer.cacheDir` at a new directory) so the release's tools are pulled.
PREDICT pins its image to the release (`orf-predict:v0.2.0`), because the model
and its cutoff must match the pipeline; the task stops if it finds an older
`predict.py`.
