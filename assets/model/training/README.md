# Training the ORF classifier

Code behind the two models PREDICT ships in `modules/predict/model/` (xorf 0.2.0 and later):

| Model | Directory | Features | Trees | Cutoff |
|---|---|---:|---:|---:|
| v3.2 (default) | `modules/predict/model/` | 63 | 625 | 0.2943 |
| v3.2 without TransAID (`skip_transaid`) | `modules/predict/model/no_transaid/` | 59 | 650 | 0.3312 |

## Files

- `train_v3.py`: the training CLI.
- `orfclf/`: `features.py` builds every feature (`modules/predict/predict.py` inlines it), `data.py` loads the
  candidate tables and labels them, `model.py` holds the grouped cross-validation, tree count, cutoff, refit and
  bundle I/O, `evaluation.py` the metrics.
- `final_params.json`: the XGBoost hyper-parameters (a 25-trial search on hg38 + mm39), used unchanged for both
  models.
- `train.sh`: the exact commands for both models.
- `requirements.txt`: the library versions the models were trained with.
- `logs/`: the logs of `train.sh` for both models.

## Data

`training_data.tar.gz` (191 MB, md5 `6c0a55a455307f4b90326e736ce7686b`) is not in the repository; it will be
deposited on Zenodo. Unpacked in this directory it creates `data/<species>/`:

- `coding/raw.xorf.tsv` and `noncoding/raw.xorf.tsv`: xorf's raw candidate tables (every candidate ORF with its tool
  scores, from runs with `predict_keep_raw: true`) for coding and lncRNA transcripts;
- `coding/<species>.coding.dataset.latest.bed`: the coding transcripts as BED12, with the annotated CDS as
  thickStart-thickEnd.

| Species | Coding candidates | Coding transcripts | Noncoding candidates |
|---|---:|---:|---:|
| hg38 | 534,758 | 19,248 | 859,350 |
| mm39 | 605,111 | 19,548 | 604,832 |
| HLratNor8 | 890,012 | 29,277 | 419,111 |

A candidate is positive when its span equals its transcript's annotated CDS; every other candidate of a coding
transcript and every candidate of a noncoding transcript is negative. Candidates of the UNMASK branch that repeat a
masked-branch candidate are dropped after the features are built, leaving 3,906,343 rows with 67,990 positives in
116,650 gene groups.

## Run

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
tar -xzf training_data.tar.gz
PYTHON=.venv/bin/python ./train.sh    # or ./train.sh v3.2 / ./train.sh v3.2_no_transaid
```

`train.sh` writes `model_v3.2/` and `model_v3.2_no_transaid/`, each with `orf_classifier_v3.bundle.json` and
`orf_classifier_v3.xgb.json` (copy them to `modules/predict/model/` and `modules/predict/model/no_transaid/`), the
cross-validation report, feature importances, out-of-fold predictions and `training.log`. With the versions in
`requirements.txt` and 8 threads per model (the default, `N_JOBS`), both models come out byte-identical to the
shipped ones (checked 2026-10-26; about 30 min and 9.3 GB of memory each).

## How the models are trained

- Features: per-ORF tool scores and ORF structure, plus context from the candidate's family (same stop codon) and
  transcript; a missing tool score stays missing (NaN). Both models leave out BLAST identity and e-value (self-hits
  leak the label) and the BLAST boundary features. v3.2 also leaves out the per-candidate TransAID scores but keeps
  whether TransAID paired an ATG with the candidate's stop (per family and per transcript); the `skip_transaid`
  model leaves out every TransAID feature.
- Five-fold cross-validation grouped by upper-cased gene symbol (StratifiedGroupKFold), so isoforms and most
  same-named orthologs share a fold, with early stopping on an inner grouped split. The final model is refit on all
  rows with the tree count of the best out-of-fold average precision.
- Cutoff: the lowest score at which, out of fold, the share of lncRNA transcripts with a call averaged over the
  three species is at most 0.8%. The F1-optimal cutoff is kept in the bundle as metadata.
