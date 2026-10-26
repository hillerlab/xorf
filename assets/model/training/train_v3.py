#!/usr/bin/env python3
"""Train the v3 ORF classifier for xorf.

Example
-------
./train_v3.py \
    --dataset hg38:coding:hg38.coding.raw.xorf.tsv:hg38.coding.dataset.latest.bed \
    --dataset hg38:noncoding:hg38.noncoding.raw.xorf.tsv \
    --dataset mm39:coding:mm39.coding.raw.xorf.tsv:mm39.coding.dataset.latest.bed \
    --dataset mm39:noncoding:mm39.noncoding.raw.xorf.tsv \
    --outdir model_v3

Design (see README.md for the reasoning):
  * positives = candidate equal to the annotated CDS; negatives = every other
    candidate of coding transcripts + every candidate of noncoding transcripts;
  * per-ORF + family/transcript context features (orfclf/features.py), NaN kept;
  * StratifiedGroupKFold grouped by upper-cased gene symbol (isoforms and most
    same-named human/mouse orthologs share a fold), early stopping on an inner
    grouped split, threshold chosen on out-of-fold predictions;
  * threshold: by default the specificity anchor (--threshold-objective lnc_rate),
    the lowest cutoff at which at most --target (0.8%) of each species' lncRNA
    transcripts get a call, averaged over species; the F1 cutoff is kept as metadata;
  * final model refit on all rows with the CV tree count.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orfclf import model as M  # noqa: E402
from orfclf.data import Source, load_features  # noqa: E402
from orfclf.evaluation import evaluate, prf, top1_mask  # noqa: E402

V2_THRESHOLD = 0.6418367028236389
log = logging.getLogger("orfclf")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", action="append", required=True,
                    help="species:coding:<xorf raw tsv>:<annotation bed12> or species:noncoding:<xorf raw tsv>")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--params-json", type=Path, help="JSON file overriding XGBoost hyper-parameters")
    ap.add_argument("--drop-features", default="",
                    help="comma-separated substrings; matching features are excluded (=name: exact name only)")
    ap.add_argument("--threshold-objective", choices=["lnc_rate", "f1", "fbeta", "precision", "recall"],
                    default="lnc_rate")
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--target", type=float,
                    help="target for lnc_rate (share of lncRNA transcripts called; default 0.008) "
                         "and precision/recall objectives")
    ap.add_argument("--no-final", action="store_true", help="only cross-validate")
    return ap.parse_args()


def setup_logging(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    log.setLevel(logging.INFO)
    for h in (logging.StreamHandler(sys.stdout), logging.FileHandler(outdir / "training.log", mode="w")):
        h.setFormatter(fmt)
        log.addHandler(h)


def subset_report(name, mask, y, p, tx, thr, prob_v2=None) -> list[dict]:
    rows = []
    yy, pp, tt = y[mask], p[mask], tx[mask]
    if yy.sum() > 0:
        met = evaluate(yy, pp, tt, threshold=thr)
        rows.append({"subset": name, "model": "v3", **_flat(met)})
        if prob_v2 is not None and not np.isnan(prob_v2[mask]).any():
            met = evaluate(yy, prob_v2[mask], tt, threshold=V2_THRESHOLD)
            rows.append({"subset": name, "model": "v2", **_flat(met)})
    else:  # noncoding: every call is a false positive
        for model, scores, t in [("v3", pp, thr), ("v2", None if prob_v2 is None else prob_v2[mask], V2_THRESHOLD)]:
            if scores is None or np.isnan(scores).any():
                continue
            calls = scores >= t
            per_tx = pd.Series(calls).groupby(tt).sum()
            rows.append({
                "subset": name, "model": model, "n": int(len(yy)),
                "at_threshold.threshold": t, "at_threshold.fp": int(calls.sum()),
                "noncoding_tx": int(per_tx.size),
                "noncoding_tx_with_call_frac": float((per_tx > 0).mean()),
                "fp_per_1k_tx": float(1000 * calls.sum() / per_tx.size),
            })
    return rows


def production_rule(y, p, tx, thr, max_predictions) -> dict:
    """xorf keeps candidates with prob >= min_score, top-N per transcript."""
    s = pd.DataFrame({"y": y, "p": p, "tx": tx})
    s = s[s.p >= thr].sort_values("p", ascending=False, kind="stable")
    kept = s.groupby("tx", sort=False).head(max_predictions)
    yhat = np.zeros(len(y), dtype=int)
    yhat[kept.index.to_numpy()] = 1
    return prf(np.asarray(y), yhat)


def _flat(d, prefix=""):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out.update(_flat(v, f"{prefix}{k}."))
        else:
            out[f"{prefix}{k}"] = v
    return out


def main() -> None:
    args = parse_args()
    setup_logging(args.outdir)
    t_start = time.time()
    sources = [Source.parse(s) for s in args.dataset]
    log.info("Datasets: %s", [f"{s.species}:{s.kind}" for s in sources])

    data, X, ids = load_features(sources)
    drop = [s for s in args.drop_features.split(",") if s]
    features = [c for c in X.columns if not any(c == s[1:] if s.startswith("=") else s in c for s in drop)]
    X = X[features]
    y = data["label"].to_numpy().astype(int)
    groups = ids["gene"].to_numpy()
    tx = (data["species"] + ":" + ids["ctx"]).to_numpy()
    log.info("Matrix: %d rows x %d features | positives=%d (%.2f%%) | genes=%d",
             len(X), X.shape[1], y.sum(), 100 * y.mean(), len(np.unique(groups)))

    params = dict(M.DEFAULT_PARAMS)
    if args.params_json:
        params.update(json.loads(args.params_json.read_text()))
    log.info("XGBoost params: %s", params)

    cv = M.cross_validate(X, y, groups, params, args.n_splits, args.n_jobs, args.seed)
    thr_f1 = M.choose_threshold(y, cv.oof, "f1")
    if args.threshold_objective == "lnc_rate":
        is_nc = data["source"].str.endswith(":noncoding").to_numpy()
        nc_tx = np.where(is_nc, (data["species"] + ":" + ids["tx"]).to_numpy(), None)
        target = 0.008 if args.target is None else args.target
        thr_sel = M.lnc_rate_threshold(cv.oof, nc_tx, data["species"].to_numpy(), target)
        log.info("Threshold (lncRNA call rate %.2f%% on out-of-fold): %.4f | per species %s | F1 cutoff %.4f",
                 100 * target, thr_sel["threshold"],
                 {sp: round(100 * r, 3) for sp, r in thr_sel["lnc_rate_by_species"].items()}, thr_f1["threshold"])
    else:
        thr_sel = M.choose_threshold(y, cv.oof, args.threshold_objective, args.beta, args.target)
        log.info("Threshold (%s on out-of-fold): %.4f | P=%.4f R=%.4f",
                 args.threshold_objective, thr_sel["threshold"], thr_sel["precision"], thr_sel["recall"])
    thr_sel["f1"] = thr_f1
    thr = thr_sel["threshold"]

    prob_v2 = pd.to_numeric(data.get("prob_coding"), errors="coerce").to_numpy() if "prob_coding" in data else None
    rows = subset_report("all", np.ones(len(y), bool), y, cv.oof, tx, thr, prob_v2)
    for src in data["source"].unique():
        rows += subset_report(src, (data["source"] == src).to_numpy(), y, cv.oof, tx, thr, prob_v2)
    for sp in data["species"].unique():
        rows += subset_report(f"{sp}:all", (data["species"] == sp).to_numpy(), y, cv.oof, tx, thr, prob_v2)
    report = pd.DataFrame(rows)
    report.to_csv(args.outdir / "cv_report.tsv", sep="\t", index=False)

    prod = {f"top{n}": production_rule(y, cv.oof, tx, thr, n) for n in (1, 2, 3, 5)}
    for k, v in prod.items():
        log.info("CV production rule (p>=thr, %s per transcript): P=%.4f R=%.4f F1=%.4f",
                 k, v["precision"], v["recall"], v["f1"])
    for _, r in report.iterrows():
        if pd.notna(r.get("at_threshold.f1")):
            log.info("CV %-22s %s: AP=%.4f P=%.4f R=%.4f F1=%.4f top1-acc=%.4f",
                     r["subset"], r["model"], r.get("ap", np.nan), r["at_threshold.precision"],
                     r["at_threshold.recall"], r["at_threshold.f1"], r.get("tx_top1_accuracy", np.nan))
        else:
            log.info("CV %-22s %s: FP calls=%d, lncRNA tx with a call=%.4f%%",
                     r["subset"], r["model"], r["at_threshold.fp"], 100 * r["noncoding_tx_with_call_frac"])

    oof = pd.DataFrame({
        "id": data["id"], "species": data["species"], "source": data["source"],
        "tx": ids["ctx"], "gene": ids["gene"], "label": y, "fold": cv.fold, "prob_v3_oof": cv.oof,
        "top1": top1_mask(cv.oof, tx),
    })
    if prob_v2 is not None:
        oof["prob_v2"] = prob_v2
    oof.to_parquet(args.outdir / "oof_predictions.parquet", index=False)
    err = oof.assign(call=oof.prob_v3_oof >= thr)
    err[(err.call) & (err.label == 0)].sort_values("prob_v3_oof", ascending=False).to_csv(
        args.outdir / "oof_false_positives.tsv", sep="\t", index=False)
    err[(~err.call) & (err.label == 1)].sort_values("prob_v3_oof").to_csv(
        args.outdir / "oof_false_negatives.tsv", sep="\t", index=False)

    meta = {
        "threshold_selection": thr_sel,
        "params": params,
        "cv": {
            "n_splits": args.n_splits, "grouping": "upper-cased gene symbol",
            "best_iterations": cv.best_iterations,
            "tree_count": cv.fold_metrics[0]["tree_count"],
            "oof_ap_by_trees": {str(k): v for k, v in cv.fold_metrics[0]["oof_ap_by_trees"].items()},
            "production_rule": prod,
            "report": report.replace({np.nan: None}).to_dict(orient="records"),
        },
        "training_data": {
            "datasets": [f"{s.species}:{s.kind}:{s.tsv}" + (f":{s.bed}" if s.bed else "") for s in sources],
            "rows": int(len(y)), "positives": int(y.sum()),
            "per_source": data.groupby("source")["label"].agg(["size", "sum"]).astype(int).to_dict(orient="index"),
        },
        "seed": args.seed,
    }

    if not args.no_final:
        n_rounds = M.final_rounds(cv)
        log.info("Refitting on all rows with %d trees", n_rounds)
        final = M.fit_final(X, y, params, n_rounds, args.n_jobs, args.seed)
        meta["n_estimators"] = n_rounds
        bundle = M.save_bundle(final, features, thr, meta, args.outdir)
        imp = pd.Series(final.get_booster().get_score(importance_type="total_gain")).reindex(features).fillna(0)
        imp.sort_values(ascending=False).rename("total_gain").to_csv(args.outdir / "feature_importance.tsv", sep="\t")
        log.info("Top features: %s", imp.sort_values(ascending=False).head(12).round(0).to_dict())
        log.info("Saved %s", bundle)
    else:
        (args.outdir / "cv_meta.json").write_text(json.dumps(meta, indent=2, default=M._json_default))
    log.info("Done in %.1f min", (time.time() - t_start) / 60)


if __name__ == "__main__":
    main()
