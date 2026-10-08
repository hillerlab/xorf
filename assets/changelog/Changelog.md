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
        CHANGELOG
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

All notable changes to this project are documented below.

---

## [0.1.1] - 2026-10-08

### Fixes
- **`predict.py` handles empty inputs**: chunks with no DIAMOND/MMseqs2 hits or RNAsamba rows (e.g. ERCC spike-ins, small rescued chunks) now log and return an empty frame instead of crashing in pandas parsing. Module version bumped to `0.0.24`.

### Infrastructure
- Manifest updated: pipeline version `0.1.1`.

## [0.1.0] - 2026-09-02

### Features
- **CPU ESMFold2-Fast pLDDT on BLAST**: `orf blast --esm` / `-E` (and the pipeline `esm` parameter, default off) folds each BLAST peptide with ESMFold2-Fast on CPU in FP32 and appends mean pLDDT (0–100) as `esm_plddt`, immediately before `blast_reference_id`. Inner `*` residues are rewritten to `X`, a terminal `*` is stripped, and sequences longer than 2048 aa are skipped with a warning and scored 0. No GPU path. The helper is embedded in `orf` (`esmfold2_fast.py`); pins are ESMFold2-Fast `c6c7958…` and ESMC-6B `650e98b…`. Module version bumped to `0.0.30`.
- **`DOWNLOAD_ESMFOLD_WEIGHTS`**: When `esm` is true, the pipeline prefetches the pinned Hugging Face hub cache once (~26 GB: Fast + ESMC-6B backbone) instead of letting every scattered BLAST task download it. The cache is published to `${outdir}/00_weights/esmfold_cache` and is not deleted. BLAST stages that directory and sets `HUGGINGFACE_HUB_CACHE` / `HF_HUB_OFFLINE=1` so the helper loads offline.
- **`esmfold_local_weights`**: Path to an existing hub cache (`models--biohub--ESMFold2-Fast` and `models--biohub--ESMC-6B`). Set this with `esm` to skip the download and reuse a previous publish dir or `$HF_HOME/hub`. Missing paths, files, or a raw `HF_HOME` (parent of `hub/`) fail in `validateRun()`.

### Fixes
- **`predict.py` BLAST reader**: Accepts both the legacy 25-column table and the optional 26-column `esm_plddt` table. `esm_plddt` is kept out of `FEATURES` / `ORDER` (benchmark metadata until the classifier is retrained). The `prefix` split on `.p` is now a regex (`r"\.p"`), which is the documented pandas behaviour. Module version bumped to `0.0.23`.

### Infrastructure
- BLAST container and conda env move to Python 3.12 with CPU-only PyTorch, pinned `esm` (`evolutionaryscale/esm@827ec12…`), and CUDA cuequivariance / NVIDIA wheels uninstalled so ESMFold2 uses the PyTorch fallback.
- BLAST is now `process_high_memory` (4 CPUs, 30 GB, 10 h). The unused 200 GB memory-only overlay is replaced; other `process_low` modules stay at 10 GB / 4 h. The BLAST `time = 2.h` override is removed so the 10 h label actually applies. CPU FP32 folding peaked at ~26 GB RSS in the SCAF8 bench.
- On SLURM, BLAST `scratch` is disabled when `esm` is on so the ~26 GB cache is not copied into every array task. `DOWNLOAD_ESMFOLD_WEIGHTS` is `process_single` with a 12 h wall time.
- New module `src/modules/esmfold/main.nf`; `esmfold2_fast.py --prefetch CACHE_DIR` is the single source of model/revision pins. Help, `params.json`, `usage.md`, and `README.md` document `esm` and `esmfold_local_weights`.
- Manifest updated: pipeline version `0.1.0`.

## [0.0.47] - 2026-09-01

### Features
- **MMseqs2 BLAST engine**: `orf blast --engine mmseqs2` searches with MMseqs2 instead of DIAMOND (still the default). MMseqs2 is asked for the same 11 columns DIAMOND already writes (`query,pident,qlen,tlen,alnlen,qstart,qend,tstart,tend,evalue,target`), so the existing `BlastRecord` parser, best-hit selection, and downstream classifier features are unchanged. Protein search (`--search-type 1`), e-value `1e-10`, and e-value-sorted hits match the DIAMOND filters. Passing a `.dmnd` path with `--engine mmseqs2` fails fast. Module version bumped to `0.0.29`.
- **`--engine` pipeline parameter**: `diamond` (default) or `mmseqs2`. With `mmseqs2` the pipeline never downloads the zenodo DIAMOND archive: it uses `raw_database` (SwissProt FASTA), optionally merges a FASTA `custom_database`, and indexes once. A `.dmnd` `custom_database` with `engine` `mmseqs2` is rejected.

### Fixes
- **`--database` / `--downstream-flank` short-option clash**: both used `-d`, which made `orf blast --help` panic in debug builds. `--database` is now `-D`; `-d` stays `--downstream-flank`, matching the other subcommands.
- **Search thread count**: DIAMOND and MMseqs2 were hardcoded to `--threads 8` while BLAST is a 4-cpu `process_low` task. Both tools now take `orf --threads N` (the BLAST process passes `task.cpus`).

### Infrastructure
- New module `MMSEQS_CREATEDB` (`src/modules/mmseqs/createdb/`): `createdb` + `createindex` once per run, emits a directory so BLAST chunks reuse the same target instead of rebuilding SwissProt on every chunk.
- The BLAST container and conda env install `mmseqs2`; container tests check `mmseqs version`.
- BLAST `versions.yml` records mmseqs; help, `params.json`, and `usage.md` document `engine`.
- Manifest updated: pipeline version `0.0.47`.

## [0.0.46] - 2026-08-26

### Fixes
- **`XORF.out.files` is now the polished HQ**: With `do_polishing` (the default), the subworkflow still ran `DETACH_DUPLICATES` → `ISOTOOLS_TRUNCATION_DETECTOR` → `STRIP_TRUNCATIONS`, but `files` kept pointing at `CONCAT_RENAMED`. Downstream pipelines (xasm polish, the completion email) therefore started as soon as concat finished, in parallel with truncation stripping, and consumed predictions that still contained `#DU` duplicates and 3′UTR-truncated ORFs. `files` is now `STRIP_TRUNCATIONS.out.hq` joined back to the concat TSV, so consumers wait on `04_results`. When polishing is off, `files` is unchanged (concat / renamed output).

### Added
- **`XORF.out.truncations`** and **`XORF.out.duplicates`**: the STRIP discard BED and the DETACH duplicates BED. Empty channels when `do_polishing` is false, so callers can concat first-pass and second-pass discards without a flag check.

### Infrastructure
- Manifest updated: pipeline version `0.0.46`.

## [0.0.45] - 2026-08-18

### Fixes
- **Empty `predictions.tsv` is now truly empty**: When no ORF passes the coding-probability threshold, `predict.py` writes a zero-byte `predictions.tsv` (previously a header-only file) alongside the empty `predictions.bed`, so downstream consumers can no longer mistake a header line for a prediction record. A log line reports the created empty files. Module version bumped to `0.0.22`.

### Infrastructure
- Manifest updated: pipeline version `0.0.45`.

## [0.0.44] - 2026-08-18

### Fixes
- **`predict.py` handles runs without predictions**: When no ORF passes the coding-probability threshold, `predict.py` now writes an empty (header-only) `predictions.tsv` and an empty `predictions.bed` instead of aborting `map_to_blocks`. The `PREDICT` process keeps its non-optional output contract, so runs with zero predictions complete cleanly instead of failing the task on missing output files. Module version bumped to `0.0.21`.
- **Prefix-aware rescue filenames**: The `orf` chunker names the per-chunk rescue BED/FASTA files with the chunk prefix, so rescue files from different chunks can no longer collide in the shared `rescue/` directory. Module version bumped to `0.0.28`.

### Infrastructure
- Manifest updated: pipeline version `0.0.44`.

## [0.0.43] - 2026-08-18

### Fixes
- The `orf` module no longer panics when the PSAURON predictions file is missing or empty: `read_psauron` returns an empty prediction map and logs a `WARN` instead of aborting the run. Module version bumped to `0.0.27`.

### Infrastructure
- Manifest updated: pipeline version `0.0.43`, and `name`/`homePage` now point to `hillerlab/xorf`.

## [0.0.42] - 2026-08-17

### Features
- **Flank rescue for transcripts at sequence edges**: Transcripts whose flanking context cannot be extracted — because they sit within 1000 bp of a chromosome start or end, as happens routinely on unplaced scaffolds — were previously discarded before reaching translationAI. They are now rescued: the chunker retries extraction with 0 bp flanks and routes them through a dedicated rescue branch (tagged `upstream: 0`, `downstream: 0`) into the same candidate, prediction and polishing pipeline as regular transcripts. Flank-aware modules (`TRANSLATION`, `RNASAMBA`, `BLAST`) now read their flanks from the record metadata instead of a single process-wide value, so rescued transcripts are translated, classified and BLASTed without flank-dependent coordinate remapping. The `orf` chunker reports a `WARN` for every rescued record.

### Infrastructure
- The `CHUNKER` process emits two new optional outputs, `rescue_regions` and `rescue_sequences`, for records extracted with 0 bp flanks (written under `tmp/rescue/`), and the `XORF` subworkflow mixes them into `ch_pairs` in every run mode (plain, seleno, `run_only_on` mask/unmask).
- The `RNASAMBA` module no longer assumes a flank-stripped FASTA always exists: when flanks are 0 the input sequence is copied into place instead of failing the `mv`.
- The `orf` chunker gained unit tests covering flank underflow and out-of-bounds extraction with and without rescue.
- Manifest version updated to `0.0.42`.

## [0.0.41] - 2026-08-04

### Breaking Changes
- **`custom_database` is now interpreted by file extension**: Previously any gzipped `custom_database` file was treated as a compressed DIAMOND database, and every other file was passed through to BLAST unchanged. The pipeline now inspects the extension and routes accordingly: `.dmnd` and `.dmnd.gz` replace the default database, FASTA inputs (`.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`) are treated as protein sequences to merge into the reference set, and any other extension aborts the run at startup. Runs relying on an unlisted extension will now fail fast with an explanatory message instead of silently reaching BLAST with a file it cannot use.

### Features
- **Appending custom proteins to the default database**: `custom_database` now also accepts FASTA files. When one is given, xorf downloads the raw UniProt/SwissProt sequences (see `raw_database` below), appends the custom sequences to them, and rebuilds a DIAMOND database on the fly with the new `DIAMOND_MAKEDB` step. The merged database is then used by the BLAST stage, so custom proteins — species-specific proteomes, novel isoforms, experimental sequences — are searched alongside the curated reference instead of replacing it. DIAMOND databases (`.dmnd` / `.dmnd.gz`) continue to replace the default database entirely.
- **`raw_database` parameter**: New parameter, defaulting to the UniProt/SwissProt FASTA download, that controls which raw sequences custom FASTA files are appended to. Point it at your own source (for example the full UniProt release instead of SwissProt) to change the reference set.
- **`--run_only_on`, `--run_only_mode`, `--run_only_target`**: The selenocysteine masking step has always produced two prediction branches — one on the masked genome and one on the original genome. These parameters restrict the run to a single branch. `run_only_on` switches the restriction on, `run_only_mode` selects the branch (`mask` or `unmask`), and `run_only_target` selects which regions to keep (`intersect`, only regions overlapping selenocysteine sites, or `exclude`, only regions that do not). Invalid combinations — `run_only_on` without a `selenocysteine_sites` file, or an unrecognised mode or target — stop the pipeline with a clear error instead of silently misbehaving. The default (`run_only_on: false`) is unchanged: both branches are predicted.
- **User documentation**: A new `assets/docs/usage.md` guide covers quick-start recipes (params file, command line, smoke test), every parameter grouped by purpose, the output layout and `counts.tsv` columns, profiles, running on SLURM clusters via `assets/hpc/xorf.sh`, and a FAQ.

### Infrastructure
- New modules: `DIAMOND_MAKEDB` (`src/modules/diamond/makedb/`) builds DIAMOND databases from FASTA, handling gzipped input and optional taxonomy files; `FASTA_MERGE` (`src/modules/diamond/merge/`) concatenates the raw and custom FASTA files, decompressing gzipped inputs; `BEDTOOLS_EXCLUDE` (`src/modules/bedtools/exclude/`) provides the region exclusion used by `run_only_target`, complementing the existing `BEDTOOLS_INTERSECT`.
- The `XORF` subworkflow now routes chunks through explicit masked and unmasked channels, instantiating `BEDTOOLS_INTERSECT` and `BEDTOOLS_EXCLUDE` separately for each branch (`BEDTOOLS_INTERSECT_MASKED`, `BEDTOOLS_INTERSECT_UNMASKED`, `BEDTOOLS_EXCLUDE_MASKED`, `BEDTOOLS_EXCLUDE_UNMASKED`).
- Database-preparation step versions (`FASTA_MERGE`, `DIAMOND_MAKEDB`) are now recorded in `xorf.versions.yml`.
- All module and subworkflow license headers were switched from Apache-2.0 to GPL-3.0, and the modules that previously carried no header at all (`bedtools/intersect`, `detach`, `gunzip`, `rename`, `strip`, `untar`) received new Hiller Lab copyright blocks.
- The `CHUNKER` stub was fixed to create its scratch directory with `mkdir -p`, so stub runs (`-stub-run`) no longer fail on the chunking step.
- The help output and run banner now document the `run_only_*` parameters and the `custom_database` extension behaviour, and `params.json` was extended with the new parameters.
- `README.md` now links the usage guide; manifest version updated to `0.0.41`.

## [0.0.40] - 2026-07-31

### Breaking Changes
- **Soft-masking of input sequences removed**: Sequence parsing now uppercases every base at read time, so lowercase (soft-masked) input no longer slips through codon matching. Start and stop codon checks in the translation step are case-normalised before comparison, and the translation table is applied to uppercased codons. Pipelines fed with lowercase reference sequences will now produce codon-aware results instead of silently emitting `X` for otherwise valid codons.
- **`orf` bumped to v0.0.26**: The translation component was updated to the case-normalised codon handling described above, replacing the previous soft-mask-dependent behaviour.
- **Full-length predictions now emitted from the renaming stage**: The workflow's `files` output no longer points at the raw concatenation. When renaming is active, the emitted predictions come from the renamed, full-length set; the merged renamed files are returned only when the joined concatenation step is kept enabled. Tooling consuming the `files` output may observe different identifiers and content than in 0.0.39.

### Features
- **`--do_polishing` switch**: New parameter, enabled by default, controlling the polishing stage. When set to `false`, the duplicate-detection and truncation-stripping steps (`DETACH_DUPLICATES`, `ISOTOOLS_TRUNCATION_DETECTOR`, `STRIP_TRUNCATIONS`) are skipped entirely and the full-length predictions pass through untouched, retaining their original duplicates and truncations.
- **`--skip_joined_concat` switch**: New parameter, disabled by default. When enabled, the joined concatenation step (`CONCAT_RENAMED` / `CONCAT_RENAMED_RAW`) is skipped, and the renamed per-region predictions are forwarded directly as the full-length set instead of being merged into a single file first.
- **`skip_netstart` passed as an explicit input**: The `JOIN` process now receives `skip_netstart` as a dedicated argument instead of reading the global parameter directly, making the net-merging behaviour explicit per invocation and removing the last hard dependency on pipeline-level state inside the module.
- **Graceful completion reporting**: When the pipeline completes successfully but no final BED files are found under `04_results`, the completion handler now falls back to reporting the intermediate BED files under `02_merged` before warning, so successful runs without a polished result set are no longer misreported.

### Infrastructure
- Manifest version updated to `0.0.40`.
- `nextflow.config` extended with documented "Joining options" and "Polishing options" sections; the `test` profile enables `skip_joined_concat` and `do_polishing` alongside the existing settings.
- Pipeline banner now credits The Hiller Lab at the Senckenberg Research Institute.
- `.gitattributes` now tags `*.config` files as Nextflow for correct language rendering.
- `assets/pipeline/xorf.mermaid` updated with the `BEDTOOLS_INTERSECT`, `UNMASKED_CANDIDATES`, and `PSAURON` nodes, reflecting the current pipeline topology.

## [0.0.39] - 2026-07-22

### Breaking Changes
- **`rename_deactivate` is now effective**: The `rename_deactivate` flag introduced in 0.0.38 is properly wired end-to-end. When enabled, both the main and raw renaming stages are skipped entirely; the concatenated BED and TSV pass through unchanged, with no identifier rewriting or tag injection. This was previously declared but not enforced in the subworkflow logic.
- **Reworked prediction renaming process**: The renamer was substantially rewritten to parse structured input identifiers (`ROOT_ORF.N[@M][#DU]` or `ROOT.pN[@M][#DU]`) and decompose them into explicit metadata tags. The output identifier now carries the ORF number (`#OR`), translation index (`#TI`), inner number (`#IN`), and duplicate marker (`#DU`) as separate structural tokens rather than embedding them opaquely. The naming layout is now `[PREFIX#...]ROOT[#TAG]...[#SC{SCORE}][#DU][#OR{N}][#TI{N}][#IN{M}][@PROTEIN]`, replacing the previous format.
- **Tag conversion from raw identifiers**: ORF metadata that was previously encoded in the prediction identifier string is now explicitly converted into tags during renaming. The `#DU` duplicate marker, ORF number, translation index, and inner number are parsed from the input id and emitted as independent tags, making downstream filtering and interpretation more straightforward.
- **Duplicate ID detection in the renamer**: The renaming stage now raises an error if two distinct predictions produce the same output identifier, preventing silent data loss from hash collisions or malformed inputs.
- **Region channel uses random hashes**: The region channel metadata now uses a 6-character random alphanumeric hash (`chr` field) instead of a hardcoded `'xorf'` string, avoiding metadata collisions when multiple pipeline runs share the same work directory.

### Features
- **`--rebase` preserves ORF metadata**: When content-derived hashing is enabled, only the parsed root is replaced by the short BLAKE2s hash; ORF metadata (number, translation index, inner number, duplicate marker) is preserved as visible tags in the renamed identifier. Predictions sharing a root therefore share a hash and remain distinguishable by their metadata.
- **Hash computation refactored for performance**: The `short_hashes` function now operates on integer digests and derives the required prefix width from adjacent sorted values in a single pass, avoiding repeated set construction. Peak memory usage is reduced by converting digests in place.
- **`predict.py` deduplicates BED/TSV mapping logic**: The inline merge-and-assign block that mapped prediction coordinates onto BED alignment blocks was extracted into a standalone `_map_predictions_to_bed` function. A pre-write consistency check now verifies that the BED and TSV ID sets match before writing output files, catching mismatches early.
- **Stable sort for prediction ranking**: The `groupby().sort_values()` calls in `predict.py` now use `kind="stable"` and `sort=False` on the groupby, preserving the original input order within groups and improving reproducibility across platforms.
- **Renamer configuration in test profile**: The `test` Nextflow profile now exposes all rename-related parameters (`rename_deactivate`, `rename_rebase`, `rename_append_orf_score`, `rename_append_protein_name`, `rename_custom_prefix`) with explicit defaults, ensuring reproducible test runs.
- **`params.json` extended with output options**: The `params.json` schema now includes the full set of renaming parameters under a documented output options section.

### Infrastructure
- The `CONCAT_RENAMED` and `CONCAT_RENAMED_RAW` modules now propagate metadata from the rename stage through the collect/map chain, enabling the merged channel to carry the correct `meta.id` for downstream directory naming instead of a hardcoded `'renamed'` identifier.
- The subworkflow now passes `rename_deactivate` as a dedicated input to the `XORF` workflow and threads it through both the main and raw renaming branches.
- Added a `randomHash()` utility function in both `main.nf` and the subworkflow to generate consistent random hashes for metadata fields.

## [0.0.38] - 2026-07-20

### Breaking Changes
- **Prediction renaming stage**: Introduced a dedicated `RENAME_PREDICTIONS` step that rewrites the shared `id` column of every BED12 file and its companion TSV with human-readable, information-rich names. Names follow the `[PREFIX#...]ROOT[#TAG]...[#SC{SCORE}][@PROTEIN]` layout, so downstream identifiers now carry masking state, optional tags, the ORF score, and the BLAST reference in a single token. This changes the identifiers emitted by the pipeline and is not backward compatible with names produced by earlier versions.
- **Reworked output directory layout**: The results tree has been reorganised to reflect the new renaming and merging stages. `01_duplicates` and `02_results` have been renumbered to `03_duplicates` and `04_results`, and two new directories were added: `01_renamed` (per-region renamed predictions) and `02_merged` (the merged, renamed BED/TSV, with raw predictions under `02_merged/raw`). Any tooling that consumed the previous fixed paths must be updated.
- **Duplicate detection now consumes renamed output**: `DETACH_DUPLICATES` is fed from the merged, renamed channel (`CONCAT_RENAMED`) instead of the raw concatenation, ensuring the final kept records and truncation stripping operate on the renamed identifiers throughout.

### Features
- **Content-derived identifiers (`--rebase`)**: When `rename_rebase` is enabled, the full original id is replaced by a short, collision-free BLAKE2s hash so that identical predictions across runs receive identical names. The hash width is grown automatically until uniqueness is guaranteed.
- **ORF score and protein annotation**: New `rename_append_orf_score` (default: on) and `rename_append_protein_name` (default: on) options append the normalised ORF score (`prob_coding`, scaled 0-100) and the BLAST reference id to each generated name.
- **Custom prefixes and masking awareness**: Added `rename_custom_prefix` to prepend a user-defined prefix, and the renamer automatically prepends an `UNMASK` prefix for predictions coming from the unmasked branch (propagated via the `masked` metadata field). Prefixes and tags containing the reserved separators `#` or `@` are rejected with a clear error.
- **Opt-out switch**: Added `rename_deactivate` to copy inputs through unchanged, preserving the original identifiers when renaming is not desired.

### Infrastructure
- The `CONCAT` module now additionally emits standalone `bed` and `tsv` channels, allowing the renamer to receive the merged files as separate inputs.
- Added `CONCAT_RENAMED` and `CONCAT_RENAMED_RAW` aliases to merge the renamed predictions (and their raw counterparts) into the new `02_merged` output directory.
- The completion handler now reports final BED files from `04_results`, matching the updated directory layout.

### Removals
- **Dropped the email reporting script**: Removed `src/bin/email.py` and the associated HTML/plaintext summary and SMTP/mailx delivery logic. Run summaries are no longer emailed by the pipeline.

## [0.0.37] - 2026-07-18

### Breaking Changes
- **BLAST/PSAURON/NETSTART2 re-impl**: Refactored the BLAST subworkflow to integrate PSAURON alongside DIAMOND for ORF scoring. The `orf blast` command now runs both tools sequentially, appending PSAURON confidence scores and reference sequence identifiers to each BLAST record. The output format has been extended accordingly.
- **`--skip_netstart` flag**: Added the ability to bypass NETSTART2 execution entirely. When `--skip_netstart` is set, the pipeline skips the NETSTART2 step and fills the merged net file with zeroed values for netstart-derived columns. This addresses cases where NETSTART2 is unavailable or undesirable.
- **NETSTART2 made optional in `orf net`**: The `--netstart` argument is no longer required. When omitted, the net merger receives an empty netstart map, yielding the same zeroed-out behaviour.

### orf v0.0.25
- **PSAURON integration**: The `orf blast` subcommand now runs PSAURON (`psauron`) after DIAMOND. Each ORF prediction is annotated with an in-frame score from PSAURON, and the output includes a `psauron_score` column.
- **Diamond `sseqid` in output**: The DIAMOND output format now includes the `sseqid` (subject sequence ID) column, providing the reference protein identifier for each match.
- **`--keep-temp` flag**: Added the `--prefix` and `--keep-temp` flags to `orf blast`. The `--prefix` flag controls output file naming, while `--keep-temp` retains intermediate files (Diamond raw output, PSAURON output, PEP files) instead of cleaning them up.
- **Graceful empty BLAST handling**: The pipeline no longer panics when DIAMOND produces zero predictions. An empty BLAST table is handled gracefully, preventing hard crashes on sequences with no database hits.
- **Container update**: The BLAST container now pins `diamond=2.2.4` and includes `psauron=1.1.3`.

### Infrastructure
- The `test` profile now uses `custom_database` instead of `database` for consistency.
- The `database` parameter in the `test` profile no longer points to the test dataset directly; instead, `custom_database` is used.
- Fixed a bug where `params.custom_database` was not being resolved correctly in the `GUNZIP_DATABASE` and direct `fromPath` branches.

## [0.0.36] - 2026-07-18

### Breaking Changes
- **CDS mapping fix for edge start/end cases**: Fixed the `map_cds` utility in `orf` (v0.0.24) to correctly handle edge cases where the CDS start or end falls at the boundaries of the alignment. This resolves incorrect ORF classifications for transcripts with start or stop codons exactly at exon junctions.
- **Default protein database retrieval**: Implemented automatic retrieval of the default UniProt/SwissProt protein database when no custom database is provided. The pipeline downloads the database from Zenodo if absent, and caches it for subsequent runs.
- **`custom_database` port**: Exposed the `custom_database` parameter as a first-class input, allowing users to supply a local DIAMOND database instead of relying on the default download.
- **UNMASKED default branch with selenocysteines**: When selenocysteine sites are provided, the pipeline now routes through an UNMASKED default branch, ensuring that selenocysteine masking produces an output with a clean `unmasked` identifier.

### orf v0.0.24
- Inclusive `map_cds` logic updated to correctly process edge positions.

## [0.0.35] - 2026-06-07

### Breaking Changes
- **Inclusive `map_cds`**: Updated `orf` v0.0.23 to make the CDS mapping function inclusive, fixing off-by-one errors that caused certain ORF boundaries to be misclassified.

### Chores
- BLAST processes are now submitted as array jobs for better scalability.

## [0.0.34] - 2026-06-05

### Breaking Changes
- **Net bug fix for non-netstart Transaid predictions**: Fixed a bug in `orf net` v0.0.22 where the net merger would produce incorrect output for transcripts lacking NETSTART2 predictions. Transaid-only predictions now merge correctly into the final net file.

### Features
- Added additional configuration options.
- Fixed documentation typos.
- Selenocysteine masking now outputs FASTA files directly.

## [0.0.33] - 2026-05-26

### Breaking Changes
- **xorf re-implementation [pre-release]**: Major internal rework of the xorf pipeline to improve modularity and data flow.

### Chores
- Updated README documentation.
- Left clear watermarks per job for traceability.

## [0.0.32] - 2026-05-22

### Breaking Changes
- **Detach/Split/Truncation + Arrays**: Implemented support for region detaching, splitting, and truncation. Added array job support for BLAST and other compute-intensive steps, enabling parallel execution across genomic regions.

### Fixes
- Updated container registry paths.
- Fixed version printing.

## [0.0.31] - 2026-05-09

### Breaking Changes
- **Genepred lint + docs**: Added gene prediction linting utilities and expanded documentation.
- **GenomeMask + Selenocysteine masking**: Implemented `genomemask` for selenocysteine masking.
- **Empty output handling**: Empty outputs from pipeline processes are now properly ignored instead of causing downstream failures.
- **Chunker input as tuple channel**: The chunker subworkflow now accepts sequence input as a tuple channel for better metadata propagation.

### Features
- Added `stageIn`/`stageOut` mode configuration options.

### Fixes
- Fixed ghost output emission from the predict step.
- Fixed selenocysteine channel consumption to avoid premature channel closure.
- Fixed WGET call for emitted channel; added support for local database paths.

## [0.0.30] - 2026-04-28

### Breaking Changes
- **Iso-chunk v0.0.21**: Fixed FASTA input handling. Added documentation and a `.gz` handler to prevent future bugs with compressed inputs.
- Switched to `flate2` default backend (`miniz_oxide`) for improved compatibility.

## [0.0.29] - 2026-04-23

### Breaking Changes
- **Chunker `--prefix` implementation**: Added `--prefix` parameter to `orf chunk` (v0.0.20), allowing configurable output file naming.

## [0.0.28] - 2026-04-21

### Breaking Changes
- **Chunker v0.0.19**: Fixed byte-seeking logic to correctly delete empty chunks. Added process tags with names for better traceability.

### Features
- Added `predict_keep_raw` extra argument to preserve raw prediction outputs.

## [0.0.27] - 2026-04-10

### Breaking Changes
- **Predict v0.0.17**: Implemented `concat_raw` for raw output concatenation. Added prediction ID support. Raw TSV and BED output available under the `--keep-raw` (`-K`) flag.

## [0.0.26] - 2026-04-02

### Breaking Changes
- **Chunker v0.0.18**: Implemented graceful handling of empty files when `--ignore-errors` is active. Fixed underflowing error when both flanks are applied to small sequences.

## [0.0.25] - 2026-03-25

### Breaking Changes
- **Predict v0.0.15**: New fixed XGBoost model for TranslationAI predictions, improving accuracy.
- **NetStart offline tokenizer**: Replaced the online tokenizer with an offline version, avoiding maximum request limit errors when queue sizes exceed 500.

## [0.0.24] - 2026-03-19

- Updated prediction model.
- Bug fixes.

## [0.0.23] - 2026-03-12

- Fixed predict headers to include net columns.
- Bug fixes and stability improvements.

## [0.0.22] - 2026-03-05

- Fixed typo in netstart header within predictions.
- Predict v0.0.10.

## [0.0.21] - 2026-02-26

- Fixed ORF coordinates across the board.
- Net shared ORF line tab spacing fix.

## [0.0.20] - 2026-02-19

- Fixed net Dockerfile and NetStart model directory.
- Stability improvements.

## [0.0.19] - 2026-02-12

- Fixed DIAMOND ghost records.
- Predict v0.0.8.

## [0.0.18] - 2026-02-05

- Fixed DIAMOND ghost records (re-fix).
- Predict v0.0.7.

## [0.0.17] - 2026-01-29

- Predict v0.0.6: fixed TranslationAI prefix splitting.
- Bug fixes.

## [0.0.16] - 2026-01-22

- **ORF blast v0.0.13**: Fixed unique TranslationAI NMD type assignment.

## [0.0.15] - 2026-01-15

- NMD detector fixed for negative strand coordinates.
- ORF blast v0.0.12.

## [0.0.14] - 2026-01-08

- Start/stop codon comparison made case-insensitive (ATCGatcgNn handling).

## [0.0.13] - 2026-01-01

- ATCGatcgNn handling across all sequence processing steps.

## [0.0.12] - 2025-12-25

- Detached candidate generation from prediction step.

## [0.0.11] - 2025-12-18

- Invalid codons accepted as X for downstream tolerance.

## [0.0.10] - 2025-12-11

- Fixed chunking slicing logic for edge cases.

## [0.0.9] - 2025-12-04

- Verbose error reporting for TranslationAI failures.

## [0.0.8] - 2025-11-27

- RNAsamba catching empty sequences gracefully.

## [0.0.7] - 2025-11-20

- Fixed exon slicing logic in the chunker.

## [0.0.6] - 2025-11-13

- Fixed chunking bug on extending exonic end.
- ORF v0.0.2.

## [0.0.5] - 2025-11-06

- Sequence slicing logic changed; added `--ignore-errors` flag.

## [0.0.4] - 2025-10-30

- Container builder CI setup.
- Initial release infrastructure.
