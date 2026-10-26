"""Training, threshold selection, persistence and inference for the v3 model."""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import GroupShuffleSplit, StratifiedGroupKFold

from .features import build_features

log = logging.getLogger("orfclf")

BUNDLE_FORMAT = "orfclf-v3"
DEFAULT_PARAMS = dict(
    learning_rate=0.05,
    max_depth=8,
    min_child_weight=5,
    subsample=0.8,
    colsample_bytree=0.7,
    reg_lambda=2.0,
    reg_alpha=0.0,
    gamma=0.0,
    tree_method="hist",
    max_bin=256,
    objective="binary:logistic",
    eval_metric="aucpr",
)
MAX_ROUNDS = 5000
EARLY_STOPPING = 150


@dataclass
class CVResult:
    oof: np.ndarray
    fold: np.ndarray
    best_iterations: list[int]
    fold_metrics: list[dict] = field(default_factory=list)


def _model(params: dict, n_estimators: int, n_jobs: int, seed: int, early_stopping: bool) -> xgb.XGBClassifier:
    kw = dict(params)
    kw.update(n_estimators=n_estimators, n_jobs=n_jobs, random_state=seed)
    if early_stopping:
        kw["early_stopping_rounds"] = EARLY_STOPPING
    return xgb.XGBClassifier(**kw)


def fit_early_stopped(X, y, groups, params, n_jobs, seed) -> xgb.XGBClassifier:
    """Fit with early stopping on a gene-grouped 10% split of the training rows.

    The evaluation fold is never used for early stopping, so out-of-fold scores
    and the threshold derived from them are not optimistically biased.
    """
    inner = GroupShuffleSplit(n_splits=1, test_size=0.1, random_state=seed)
    tr, es = next(inner.split(X, y, groups))
    model = _model(params, MAX_ROUNDS, n_jobs, seed, early_stopping=True)
    model.fit(X.iloc[tr], y[tr], eval_set=[(X.iloc[es], y[es])], verbose=False)
    return model


def cross_validate(X, y, groups, params, n_splits, n_jobs, seed) -> CVResult:
    """Gene-grouped CV. Returns out-of-fold scores at one common tree count.

    Early stopping points vary a lot between folds (aucpr plateaus early), so
    instead of taking a median, every fold model scores its test fold at a grid
    of tree counts (multiples of 25) and the count maximising pooled out-of-fold
    AP is used. Every fold is scored at exactly that count; a fold model that
    stopped with fewer trees uses all of its own. The OOF scores (hence the
    threshold) and the final refit then correspond to the same model size.
    """
    from sklearn.metrics import average_precision_score

    skf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    fold = np.full(len(y), -1, dtype=int)
    best, boosters, splits = [], [], []
    for k, (tr, te) in enumerate(skf.split(X, y, groups)):
        t0 = time.time()
        model = fit_early_stopped(X.iloc[tr], y[tr], groups[tr], params, n_jobs, seed + k)
        fold[te] = k
        best.append(int(model.best_iteration) + 1)
        boosters.append(model.get_booster())
        splits.append(te)
        log.info("fold %d/%d: %d train / %d test rows, early stop at %d trees (%.0fs)",
                 k + 1, n_splits, len(tr), len(te), best[-1], time.time() - t0)

    grid = list(range(25, max(b.num_boosted_rounds() for b in boosters) + 1, 25))
    preds = []
    for booster, te in zip(boosters, splits):
        dm = xgb.DMatrix(X.iloc[te])
        n_trees = booster.num_boosted_rounds()
        preds.append({n: booster.predict(dm, iteration_range=(0, min(n, n_trees))) for n in grid})

    def oof_at(n: int) -> np.ndarray:
        oof = np.empty(len(y))
        for te, pr in zip(splits, preds):
            oof[te] = pr[n]
        return oof

    candidates = grid
    scores = {n: average_precision_score(y, oof_at(n)) for n in candidates}
    n_best = max(scores, key=scores.get)
    log.info("OOF AP by tree count: %s", {n: round(s, 5) for n, s in scores.items() if n % 100 == 0 or n == n_best})
    log.info("Chosen common tree count: %d (OOF AP %.5f)", n_best, scores[n_best])
    return CVResult(oof=oof_at(n_best), fold=fold, best_iterations=best,
                    fold_metrics=[{"tree_count": n_best, "oof_ap_by_trees": scores}])


def choose_threshold(y, p, objective="f1", beta=1.0, target=None) -> dict:
    prec, rec, thr = precision_recall_curve(y, p)
    prec, rec = prec[:-1], rec[:-1]
    if objective in ("f1", "fbeta"):
        b2 = 1.0 if objective == "f1" else beta**2
        score = np.where(prec + rec > 0, (1 + b2) * prec * rec / (b2 * prec + rec), 0.0)
        i = int(np.argmax(score))
    elif objective == "precision":
        ok = np.flatnonzero(prec >= target)
        if not len(ok):
            raise ValueError(f"No threshold reaches precision {target}")
        i = int(ok[np.argmax(rec[ok])])
        score = rec
    elif objective == "recall":
        ok = np.flatnonzero(rec >= target)
        i = int(ok[np.argmax(prec[ok])])
        score = prec
    else:
        raise ValueError(objective)
    return {
        "objective": objective, "beta": beta, "target": target,
        "threshold": float(thr[i]), "precision": float(prec[i]),
        "recall": float(rec[i]), "objective_score": float(score[i]),
    }


def final_rounds(cv: CVResult) -> int:
    # The refit uses the tree count the out-of-fold scores (and the threshold)
    # were measured at; the earlier x1.1 for the larger refit set was never
    # validated (OOF AP is flat to 1e-4 around the chosen count).
    return int(cv.fold_metrics[0]["tree_count"])


def lnc_rate_threshold(p, noncoding_tx, species, target: float) -> dict:
    """Specificity-anchored cutoff: the lowest cutoff at which the share of
    noncoding transcripts with a candidate at or above it, averaged over
    species (equal weight), does not exceed ``target``.

    ``p``, ``noncoding_tx`` and ``species`` are per candidate; ``noncoding_tx``
    is the transcript id for noncoding candidates and None elsewhere. A call is
    a transcript-level decision, so this anchors it on what the call costs on
    lncRNAs instead of candidate-level F1, where every alternative start of a
    coding transcript counts as a false positive. On HLratNor8 the rate of
    unseen lncRNAs called tracked the training out-of-fold target closely.
    """
    s = pd.DataFrame({"p": p, "tx": noncoding_tx, "species": species}).dropna(subset=["tx"])
    best = {sp: np.sort(g.groupby("tx")["p"].max().to_numpy()) for sp, g in s.groupby("species")}
    cuts = np.unique(np.concatenate(list(best.values())))
    # share of each species' lncRNA transcripts with max p >= cut, averaged
    rates = np.mean([1 - np.searchsorted(m, cuts, side="left") / len(m) for m in best.values()], axis=0)
    ok = cuts[rates <= target]
    t = float(ok.min()) if len(ok) else float(np.nextafter(cuts.max(), 1))
    by_species = {sp: float((m >= t).mean()) for sp, m in best.items()}
    return {"objective": "lnc_rate", "target": target, "threshold": t,
            "lnc_rate": float(np.mean(list(by_species.values()))), "lnc_rate_by_species": by_species}


def fit_final(X, y, params, n_estimators, n_jobs, seed) -> xgb.XGBClassifier:
    model = _model(params, n_estimators, n_jobs, seed, early_stopping=False)
    model.fit(X, y, verbose=False)
    return model


# ---------------------------------------------------------------------------
# persistence / inference
# ---------------------------------------------------------------------------


def save_bundle(model: xgb.XGBClassifier, features: list[str], threshold: float, meta: dict, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    model_file = outdir / "orf_classifier_v3.xgb.json"
    model.save_model(model_file)
    bundle = {
        "format": BUNDLE_FORMAT,
        "model_file": model_file.name,
        "features": features,
        "threshold": threshold,
        "xgboost_version": xgb.__version__,
        **meta,
    }
    path = outdir / "orf_classifier_v3.bundle.json"
    path.write_text(json.dumps(bundle, indent=2, default=_json_default, allow_nan=False))
    return path


def load_bundle(path: Path) -> tuple[xgb.XGBClassifier, dict]:
    bundle = json.loads(Path(path).read_text())
    if bundle.get("format") != BUNDLE_FORMAT:
        raise ValueError(f"{path} is not an {BUNDLE_FORMAT} bundle")
    model = xgb.XGBClassifier()
    model.load_model(Path(path).parent / bundle["model_file"])
    return model, bundle


def score_candidates(df: pd.DataFrame, model: xgb.XGBClassifier, bundle: dict) -> np.ndarray:
    """Probability of being the real CDS for every row of a candidate table.

    ``df`` must contain all candidates of each transcript (context features are
    computed across them).
    """
    X, _ = build_features(df, with_context=True)
    missing = [c for c in bundle["features"] if c not in X.columns]
    if missing:
        raise ValueError(f"Features missing from input: {missing}")
    return model.predict_proba(X[bundle["features"]])[:, 1]


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))
