#!/usr/bin/env bash
# Train the two ORF classifiers shipped in modules/predict/model (xORF >= 0.2.0), from this directory after
# `tar -xzf training_data.tar.gz` (which creates data/):
#   v3.2              the default model (63 features), -> modules/predict/model/
#   v3.2_no_transaid  the skip_transaid model (59 features), -> modules/predict/model/no_transaid/
# Both use hg38, mm39 and HLratNor8 (coding + noncoding), the same hyper-parameters and the specificity anchor
# (at most 0.8% of each species' lncRNA transcripts with a call, out of fold).
#   usage: ./train.sh [v3.2] [v3.2_no_transaid]    (default: both; N_JOBS threads per model, default 8)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
DATA=(
  --dataset hg38:coding:data/hg38/coding/raw.xorf.tsv:data/hg38/coding/hg38.coding.dataset.latest.bed
  --dataset hg38:noncoding:data/hg38/noncoding/raw.xorf.tsv
  --dataset mm39:coding:data/mm39/coding/raw.xorf.tsv:data/mm39/coding/mm39.coding.dataset.latest.bed
  --dataset mm39:noncoding:data/mm39/noncoding/raw.xorf.tsv
  --dataset HLratNor8:coding:data/HLratNor8/coding/raw.xorf.tsv:data/HLratNor8/coding/HLratNor8.coding.dataset.latest.bed
  --dataset HLratNor8:noncoding:data/HLratNor8/noncoding/raw.xorf.tsv
)
PARAMS=(--params-json final_params.json --threshold-objective lnc_rate --target 0.008 --seed 42 --n-jobs "${N_JOBS:-8}")
# dropped from both models: BLAST identity and e-value (self-hit leak) and the BLAST boundary features
NB=blast_pid,neg_log10_blast_evalue,blast_offset,=blast_percentage_aligned,=blast_len_frac_orf,=fam_delta_blast_percentage_aligned,=fam_is_best_blast_percentage_aligned,=tx_delta_blast_percentage_aligned,=tx_rank_blast_percentage_aligned
# the per-candidate TransAID scores; v3.2 keeps the TransAID flags (fam_has_*, famtx_has_transaid, tx_has_transaid)
TA_CAND="=transaid_start_score,=transaid_stop_score,=transaid_integrated_score,=transaid_mean_score,=transaid_start_minus_stop,=fam_delta_transaid_start_score,=fam_is_best_transaid_start_score,=fam_delta_transaid_integrated_score,=fam_is_best_transaid_integrated_score,=tx_delta_transaid_integrated_score,=tx_rank_transaid_integrated_score"
models=("$@")
[ ${#models[@]} -eq 0 ] && models=(v3.2 v3.2_no_transaid)
for m in "${models[@]}"; do
  case $m in
    v3.2) drop=$NB,$TA_CAND ;;
    v3.2_no_transaid) drop=$NB,transaid ;;
    *) echo "unknown model $m (v3.2, v3.2_no_transaid)" >&2; exit 1 ;;
  esac
  out=model_$m
  [ -e "$out" ] && { echo "$out exists; not overwriting" >&2; exit 1; }
  "$PY" train_v3.py "${DATA[@]}" "${PARAMS[@]}" --drop-features "$drop" --outdir "$out"
done
