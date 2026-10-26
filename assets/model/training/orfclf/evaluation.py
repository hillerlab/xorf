"""Metrics for the ORF classifier: ORF-level and transcript-level."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score


def best_f1_threshold(y: np.ndarray, p: np.ndarray) -> tuple[float, float, float, float]:
    """Threshold maximising F1. Returns (threshold, precision, recall, f1)."""
    prec, rec, thr = precision_recall_curve(y, p)
    prec, rec = prec[:-1], rec[:-1]  # last point has no threshold
    f1 = np.where(prec + rec > 0, 2 * prec * rec / (prec + rec), 0.0)
    i = int(np.argmax(f1))
    return float(thr[i]), float(prec[i]), float(rec[i]), float(f1[i])


def prf(y: np.ndarray, yhat: np.ndarray) -> dict[str, float]:
    tp = int(((yhat == 1) & (y == 1)).sum())
    fp = int(((yhat == 1) & (y == 0)).sum())
    fn = int(((yhat == 0) & (y == 1)).sum())
    p = tp / (tp + fp) if tp + fp else 0.0
    r = tp / (tp + fn) if tp + fn else 0.0
    f = 2 * p * r / (p + r) if p + r else 0.0
    return {"precision": p, "recall": r, "f1": f, "tp": tp, "fp": fp, "fn": fn}


def top1_mask(p: np.ndarray, tx: np.ndarray) -> np.ndarray:
    """True for the highest-scoring candidate of each transcript."""
    s = pd.DataFrame({"p": p, "tx": tx})
    idx = s.groupby("tx", sort=False)["p"].idxmax()
    m = np.zeros(len(s), dtype=bool)
    m[idx.to_numpy()] = True
    return m


def evaluate(y: np.ndarray, p: np.ndarray, tx: np.ndarray, threshold: float | None = None) -> dict:
    y = np.asarray(y).astype(int)
    p = np.asarray(p, dtype=float)
    tx = np.asarray(tx)
    out: dict = {
        "n": int(len(y)),
        "positives": int(y.sum()),
        "ap": float(average_precision_score(y, p)),
        "roc_auc": float(roc_auc_score(y, p)),
    }
    thr_opt, p_opt, r_opt, f_opt = best_f1_threshold(y, p)
    out["oracle_best_f1"] = {"threshold": thr_opt, "precision": p_opt, "recall": r_opt, "f1": f_opt}
    if threshold is not None:
        out["at_threshold"] = {"threshold": threshold, **prf(y, (p >= threshold).astype(int))}
        top = top1_mask(p, tx)
        out["top1_and_threshold"] = {"threshold": threshold, **prf(y, (top & (p >= threshold)).astype(int))}
    # transcript level: among transcripts holding a positive, is the arg-max the positive?
    s = pd.DataFrame({"y": y, "p": p, "tx": tx})
    has_pos = s.groupby("tx", sort=False)["y"].transform("max") == 1
    s2 = s[has_pos]
    top = s2.loc[s2.groupby("tx", sort=False)["p"].idxmax()]
    out["tx_top1_accuracy"] = float(top["y"].mean())
    if threshold is not None:
        calls = s.assign(c=(s.p >= threshold)).groupby("tx", sort=False)["c"].sum()
        out["calls_per_tx_mean"] = float(calls.mean())
    return out


def flat(metrics: dict, prefix: str = "") -> dict:
    row = {}
    for k, v in metrics.items():
        if isinstance(v, dict):
            row.update(flat(v, f"{prefix}{k}."))
        else:
            row[f"{prefix}{k}"] = v
    return row
