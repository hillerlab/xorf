"""Loading and labelling xorf candidate tables for training.

Labels:
  * coding transcripts: a candidate is positive when its genomic ORF span
    (``start``/``end`` in the xorf TSV, i.e. thickStart/thickEnd, stop codon
    included) equals the annotated CDS (thickStart/thickEnd) of its transcript in
    the input BED12. Every other candidate of the transcript is negative.
  * noncoding transcripts: every candidate is negative.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .features import build_features, parse_ids

log = logging.getLogger("orfclf")

BED_COLS = ["chr", "tx_start", "tx_end", "name", "score", "strand", "cds_start", "cds_end"]
# Columns that are carried along for reporting but never used as features.
KEEP_META = ["id", "chr", "start", "end", "strand", "start_codon", "prob_coding"]


@dataclass(frozen=True)
class Source:
    species: str
    kind: str  # "coding" | "noncoding"
    tsv: Path
    bed: Path | None = None

    @classmethod
    def parse(cls, spec: str) -> "Source":
        """``species:coding:table.tsv:annotation.bed`` or ``species:noncoding:table.tsv``."""
        parts = spec.split(":")
        if len(parts) == 4 and parts[1] == "coding":
            return cls(parts[0], "coding", Path(parts[2]), Path(parts[3]))
        if len(parts) == 3 and parts[1] == "noncoding":
            return cls(parts[0], "noncoding", Path(parts[2]))
        raise ValueError(f"Bad dataset spec {spec!r}; see --help")


def read_candidates(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    missing = {"id", "start", "end", "relative_orf_start", "relative_orf_end"} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    return df


def read_annotation(path: Path) -> pd.DataFrame:
    bed = pd.read_csv(path, sep="\t", header=None, usecols=range(8), names=BED_COLS)
    if bed["name"].duplicated().any():
        raise ValueError(f"{path}: duplicated transcript names")
    return bed


def label_coding(df: pd.DataFrame, bed: pd.DataFrame) -> pd.Series:
    ids = parse_ids(df)
    ann = bed.set_index("name")
    has_tx = ids["tx"].isin(ann.index)
    if not has_tx.all():
        raise ValueError(f"{(~has_tx).sum()} candidates belong to transcripts absent from the BED")
    cds_start = ids["tx"].map(ann["cds_start"])
    cds_end = ids["tx"].map(ann["cds_end"])
    return ((df["start"] == cds_start) & (df["end"] == cds_end)).astype(np.int8)


def load_features(sources: list[Source]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Candidates, feature matrix and ids, with context features as in production.

    xorf scores the masked-genome branch and the ``UNMASK#`` selenocysteine branch
    in separate PREDICT runs over their complete candidate sets, so features are
    built per species on every candidate first; only then are the UNMASK copies of
    masked candidates dropped, from rows, features and ids alike. Dropping them
    before building features changed the context of the surviving UNMASK rows.
    """
    data = load_sources(sources, drop_duplicates=False)
    feats, ids = [], []
    for _, part in data.groupby("species", sort=False):  # context never crosses species
        X_sp, ids_sp = build_features(part, with_context=True)
        feats.append(X_sp)
        ids.append(ids_sp)
    X = pd.concat(feats).loc[data.index]
    ids = pd.concat(ids).loc[data.index]
    keep = ~duplicate_mask(data)
    if (~keep).any():
        log.info("Dropping %d duplicate candidates (UNMASK copies) after building features", int((~keep).sum()))
    return (data[keep].reset_index(drop=True), X[keep].reset_index(drop=True),
            ids[keep].reset_index(drop=True))


def load_sources(sources: list[Source], drop_duplicates: bool = True) -> pd.DataFrame:
    frames = []
    for src in sources:
        df = read_candidates(src.tsv)
        if src.kind == "coding":
            bed = read_annotation(src.bed)
            df["label"] = label_coding(df, bed)
            n_tx = parse_ids(df)["tx"].nunique()
            n_pos_tx = parse_ids(df[df["label"] == 1])["tx"].nunique()
            log.info(
                "%s coding: %d candidates, %d transcripts, %d positives "
                "(annotated CDS found for %.2f%% of %d annotated transcripts)",
                src.species, len(df), n_tx, int(df["label"].sum()),
                100 * n_pos_tx / len(bed), len(bed),
            )
        else:
            df["label"] = np.int8(0)
            log.info("%s noncoding: %d candidates", src.species, len(df))
        df["species"] = src.species
        df["source"] = f"{src.species}:{src.kind}"
        frames.append(df)
    data = pd.concat(frames, ignore_index=True)
    return deduplicate(data) if drop_duplicates else data


def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """Drop exact duplicate candidates (same species, transcript, span, start codon).

    The selenocysteine-unmasked branch re-emits some ordinary candidates of the
    same transcripts; keeping both would double-count them.
    """
    dup = duplicate_mask(df)
    if dup.any():
        log.info("Dropping %d duplicate candidates (UNMASK copies)", int(dup.sum()))
    return df[~dup].reset_index(drop=True)


def duplicate_mask(df: pd.DataFrame) -> np.ndarray:
    """True for UNMASK candidates identical to a masked candidate (see deduplicate)."""
    ids = parse_ids(df)
    key = pd.DataFrame(
        {
            "species": df["species"],
            "tx": ids["tx"],
            "start": df["start"],
            "end": df["end"],
            "codon": df.get("start_codon", ""),
            "unmask": df["id"].astype(str).str.startswith("UNMASK#"),
        }
    )
    # prefer the ordinary (masked-genome) copy when both exist
    order = key.sort_values("unmask", kind="mergesort").index
    dup = key.loc[order].duplicated(["species", "tx", "start", "end", "codon"], keep="first")
    return dup.reindex(df.index).to_numpy()
