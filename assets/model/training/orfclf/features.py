"""Feature engineering for the v3 ORF classifier.

The same code runs at training time and inside xorf's PREDICT step, so every
feature here must be computable from the merged candidate table alone
(``*.xorf.tsv`` rows), with all candidates of a transcript present together.

Three feature blocks:
  * per-ORF features (tool scores, structure). Missing values stay NaN: the
    absence of a TranslationAI/TransAID/BLAST call is that tool's "vote" and
    XGBoost learns a default direction for it.
  * family context: candidates of one transcript sharing a stop codon
    (``<tx>#OR<n>`` plus its ``#IN<k>`` internal starts). This drives
    start-site selection, which is where most errors of the v2 model came from.
  * transcript context: how a candidate compares with every other candidate
    of the same transcript (length, PSAURON, BLAST, AI scores).
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

ID_PREFIX_RE = re.compile(r"^UNMASK#")

# Raw columns read from *.xorf.tsv. netstart_atg_score is deliberately absent:
# NetStart is skipped in production runs (skip_netstart=true), so the column
# is constant 0/NaN there.
NUMERIC_INPUTS = [
    "tai_start_score",
    "tai_stop_score",
    "tai_mean_score",
    "transaid_start_score",
    "transaid_stop_score",
    "transaid_integrated_score",
    "transaid_mean_score",
    "blast_pid",
    "neg_log10_blast_evalue",
    "blast_offset",
    "blast_length",
    "blast_percentage_aligned",
    "psauron_score",
    "rna_score",
    "log_orf_len",
    "relative_orf_start",
    "relative_orf_end",
    "inner_stops",
    "orf_type",
    "nmd_type",
    "blast_match",
]

# Columns that must never become features (identifiers, coordinates, labels,
# and outputs of the previous model).
NON_FEATURES = {
    "chr", "start", "end", "id", "strand", "key", "blast_reference_id",
    "start_codon", "stop_codon", "blast_evalue", "netstart_atg_score",
    "tai_combined_score", "predicted_class", "prob_noncoding", "prob_coding",
    "has_canonical_start", "has_canonical_stop",
}

# Scores compared against the best candidate of the same family / transcript.
# For start-site scores "higher is better"; blast offset uses |offset|, lower
# is better, handled separately.
FAMILY_SCORES = [
    "tai_start_score",
    "tai_mean_score",
    "transaid_start_score",
    "transaid_integrated_score",
    "psauron_score",
    "blast_percentage_aligned",
    "neg_log10_blast_evalue",
]
TRANSCRIPT_SCORES = [
    "psauron_score",
    "log_orf_len",
    "neg_log10_blast_evalue",
    "blast_percentage_aligned",
    "tai_mean_score",
    "transaid_integrated_score",
]


def parse_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Transcript / family / gene identifiers for every candidate.

    * ``tx``: transcript the candidate belongs to. Inside xorf's PREDICT step this
      is the ``prefix`` column; in the renamed tables (``SRGAP2C-202#OR16#IN2``,
      ``UNMASK#Aqp4-213#OR9#IN3``) it is parsed from ``id``.
    * ``ctx``: the candidate set context features are computed over. Rows from
      the selenocysteine-unmasked branch (``UNMASK#``) are scored by a separate
      PREDICT run in production, so they form their own context here too.
    * ``fam``: candidates of one context that share the stop codon (same
      ``relative_orf_end``): the orfipy ORF, its nested ``#IN`` starts and any
      TranslationAI-only ORF ending there.
    * ``gene``: CV grouping key, see gene_from_transcript.
    """
    raw_id = df["id"].astype(str)
    if "prefix" in df.columns:
        tx = df["prefix"].astype(str)
    else:
        core = raw_id.str.replace(ID_PREFIX_RE, "", regex=True).str.split("@").str[0]
        tx = core.str.split("#").str[0]
    out = pd.DataFrame(index=df.index)
    out["tx"] = tx
    out["ctx"] = np.where(raw_id.str.startswith("UNMASK#"), "UNMASK#" + tx, tx)
    stop = pd.to_numeric(df["relative_orf_end"], errors="coerce").astype("Int64").astype(str)
    out["fam"] = out["ctx"] + ":" + stop
    out["gene"] = gene_from_transcript(tx)
    return out


def gene_from_transcript(tx: pd.Series) -> pd.Series:
    """Gene symbol used to group isoforms (and orthologs) in cross-validation.

    ``ABC1-201`` -> ``ABC1``; LNCipedia-style ``lnc-FOO-2_3`` -> ``lnc-FOO-2``.
    Upper-cased so human/mouse orthologs (AQP4/Aqp4) share a group.
    """
    gene = tx.str.replace(r"-\d+$", "", regex=True)
    gene = gene.str.replace(r"_\d+$", "", regex=True)
    return gene.str.upper()


def _start_codon_class(codon: pd.Series) -> pd.DataFrame:
    codon = codon.astype(str).str.upper()
    return pd.DataFrame(
        {
            "start_atg": (codon == "ATG").astype(np.int8),
            "start_ctg": (codon == "CTG").astype(np.int8),
            "start_ttg": (codon == "TTG").astype(np.int8),
            "start_other": (~codon.isin(["ATG", "CTG", "TTG"])).astype(np.int8),
        },
        index=codon.index,
    )


def _stop_codon_class(codon: pd.Series) -> pd.DataFrame:
    codon = codon.astype(str).str.upper()
    return pd.DataFrame(
        {
            "stop_tga": (codon == "TGA").astype(np.int8),
            "stop_taa": (codon == "TAA").astype(np.int8),
            "stop_tag": (codon == "TAG").astype(np.int8),
        },
        index=codon.index,
    )


def per_orf_features(df: pd.DataFrame) -> pd.DataFrame:
    x = pd.DataFrame(index=df.index)
    for col in NUMERIC_INPUTS:
        if col in df.columns:
            x[col] = pd.to_numeric(df[col], errors="coerce").astype(np.float32)
        else:
            x[col] = np.float32(np.nan)

    x["orf_len_nt"] = (x["relative_orf_end"] - x["relative_orf_start"]).astype(np.float32)
    x["abs_blast_offset"] = x["blast_offset"].abs()
    x["log_abs_blast_offset"] = np.log1p(x["abs_blast_offset"])
    # Fraction of the ORF covered by the BLAST alignment (blast_length is in aa).
    x["blast_len_frac_orf"] = (3.0 * x["blast_length"] / x["orf_len_nt"]).where(x["blast_length"] > 0)
    x["tai_start_minus_stop"] = x["tai_start_score"] - x["tai_stop_score"]
    x["transaid_start_minus_stop"] = x["transaid_start_score"] - x["transaid_stop_score"]

    x = pd.concat(
        [
            x,
            _start_codon_class(df.get("start_codon", pd.Series("", index=df.index))),
            _stop_codon_class(df.get("stop_codon", pd.Series("", index=df.index))),
        ],
        axis=1,
    )
    return x


def _group_max(values: pd.Series, groups: pd.Series) -> pd.Series:
    return values.groupby(groups, sort=False).transform("max")


def _group_min(values: pd.Series, groups: pd.Series) -> pd.Series:
    return values.groupby(groups, sort=False).transform("min")


def context_features(x: pd.DataFrame, ids: pd.DataFrame) -> pd.DataFrame:
    """Family- and transcript-relative features. ``x`` is per_orf_features."""
    tx = ids["ctx"]
    fam = ids["fam"]
    c = pd.DataFrame(index=x.index)
    start = x["relative_orf_start"]
    end = x["relative_orf_end"]
    length = x["orf_len_nt"]

    # ---- family: same transcript and stop codon, different starts ----------
    c["fam_size"] = fam.map(fam.value_counts()).astype(np.float32)
    c["fam_start_rank"] = start.groupby(fam, sort=False).rank(method="min") - 1
    c["fam_start_rank_frac"] = c["fam_start_rank"] / (c["fam_size"] - 1).replace(0, np.nan)
    fam_first = _group_min(start, fam)
    c["fam_codons_from_first"] = (start - fam_first) / 3.0
    c["fam_len_ratio"] = length / _group_max(length, fam)

    is_atg = x["start_atg"].astype(bool)
    atg_start = start.where(is_atg)
    first_atg = _group_min(atg_start, fam)
    c["fam_is_first_atg"] = (is_atg & (start == first_atg)).astype(np.int8)
    c["fam_codons_from_first_atg"] = (start - first_atg) / 3.0
    # number of ATG starts upstream of this one within the family
    order = pd.DataFrame({"fam": fam, "start": start, "atg": is_atg.astype(int)})
    order = order.sort_values(["fam", "start"], kind="mergesort")
    order["atg_up"] = order.groupby("fam", sort=False)["atg"].cumsum() - order["atg"]
    c["fam_n_atg_upstream"] = order["atg_up"].reindex(x.index).astype(np.float32)

    for s in FAMILY_SCORES:
        fmax = _group_max(x[s], fam)
        c[f"fam_delta_{s}"] = x[s] - fmax
        c[f"fam_is_best_{s}"] = ((x[s] == fmax) & x[s].notna()).astype(np.int8)
        c[f"fam_has_{s}"] = fmax.notna().astype(np.int8)

    fmin_off = _group_min(x["abs_blast_offset"], fam)
    c["fam_delta_abs_blast_offset"] = x["abs_blast_offset"] - fmin_off
    c["fam_is_best_abs_blast_offset"] = (
        (x["abs_blast_offset"] == fmin_off) & x["abs_blast_offset"].notna()
    ).astype(np.int8)

    # ---- transcript: all candidates of the transcript -----------------------
    c["tx_n_cand"] = tx.map(tx.value_counts()).astype(np.float32)
    c["tx_n_fam"] = fam.groupby(tx, sort=False).transform("nunique").astype(np.float32)
    c["tx_len_rank"] = length.groupby(tx, sort=False).rank(ascending=False, method="min")
    c["tx_len_ratio"] = length / _group_max(length, tx)
    tx_extent = _group_max(end, tx)  # proxy for transcript length
    c["tx_rel_start"] = start / tx_extent
    c["tx_rel_end"] = end / tx_extent
    c["tx_utr3_proxy"] = tx_extent - end
    for s in TRANSCRIPT_SCORES:
        tmax = _group_max(x[s], tx)
        c[f"tx_delta_{s}"] = x[s] - tmax
        c[f"tx_rank_{s}"] = x[s].groupby(tx, sort=False).rank(ascending=False, method="min")

    # family-level aggregates seen from the transcript: is my stop-codon family
    # the best one of the transcript?
    fam_best_ps = _group_max(x["psauron_score"], fam)
    c["famtx_delta_psauron"] = fam_best_ps - _group_max(fam_best_ps, tx)
    fam_best_len = _group_max(length, fam)
    c["famtx_len_ratio"] = fam_best_len / _group_max(fam_best_len, tx)
    fam_blast = _group_max(x["blast_percentage_aligned"], fam)
    c["famtx_delta_blast_aligned"] = fam_blast - _group_max(fam_blast, tx)
    fam_tai = _group_max(x["tai_mean_score"], fam)
    c["famtx_has_tai"] = fam_tai.notna().astype(np.int8)
    fam_ta = _group_max(x["transaid_integrated_score"], fam)
    c["famtx_has_transaid"] = fam_ta.notna().astype(np.int8)
    c["tx_has_transaid"] = _group_max(x["transaid_integrated_score"], tx).notna().astype(np.int8)
    c["tx_has_tai"] = _group_max(x["tai_mean_score"], tx).notna().astype(np.int8)
    c["tx_has_blast"] = _group_max(x["blast_percentage_aligned"], tx).notna().astype(np.int8)
    return c.astype(np.float32)


def build_features(df: pd.DataFrame, with_context: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (feature matrix, id frame). Rows keep df's index."""
    ids = parse_ids(df)
    x = per_orf_features(df)
    if with_context:
        x = pd.concat([x, context_features(x, ids)], axis=1)
    leaked = sorted(set(x.columns) & NON_FEATURES)
    if leaked:
        raise ValueError(f"Non-feature columns leaked into the matrix: {leaked}")
    return x.astype(np.float32), ids
