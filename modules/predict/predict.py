#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.23"

import argparse
import logging
from os import PathLike
from pathlib import Path
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import xgboost as xgb

MODEL = "model/orf_classifier_model.json"
MODEL_PATH = Path(__file__).parent / MODEL
UNIQUE_TAI_TAG = "UT"
START_CODONS = ["ATG", "CTG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]
BLAST_COLS: List = [
    "chr",
    "start",
    "end",
    "id",
    "strand",
    "tai_start_score",
    "tai_stop_score",
    "relative_orf_start",
    "relative_orf_end",
    "start_codon",
    "stop_codon",
    "inner_stops",
    "orf_type",
    "nmd_type",
    "netstart_atg_score",
    "transaid_start_score",
    "transaid_stop_score",
    "transaid_integrated_score",
    "blast_pid",
    "blast_evalue",
    "blast_offset",
    "blast_length",
    "blast_percentage_aligned",
    "psauron_score",
    "blast_reference_id"
]
SAMBA_COLS: List = ["prefix", "rna_score"]
TAI_MASKING_NAN_COLS: List = [
    "tai_start_score",
    "tai_stop_score",
    "tai_mean_score",
    "tai_combined_score",
]
BLAST_MASKING_NAN_COLS: List = [
    "blast_pid",
    "neg_log10_blast_evalue",
    "blast_offset",
    "blast_percentage_aligned",
]
FEATURES: List = [
    "tai_start_score",
    "tai_stop_score",
    "tai_mean_score",
    "inner_stops",
    "orf_type",
    "nmd_type",
    "blast_pid",
    "neg_log10_blast_evalue",
    "blast_offset",
    "blast_percentage_aligned",
    "rna_score",
    "blast_match",
    "log_orf_len",
    "has_canonical_start",
    "has_canonical_stop",
    "netstart_atg_score",
    "transaid_start_score",
    "transaid_stop_score",
    "transaid_mean_score",
    "transaid_integrated_score",
]
EXTRA_COLS: List = [
    "start_ai_consensus",
    "stop_ai_consensus",
    "start_stop_symmetry",
    "ai_quality_interaction",
    "orf_quality_score",
]
MISSINGNESS_COLS: List = [
    "rna_score",
    "tai_start_score",
    "tai_stop_score",
    "tai_mean_score",
    "netstart_atg_score",
    "transaid_start_score",
    "transaid_stop_score",
    "transaid_mean_score",
    "transaid_integrated_score",
]
ORDER: List = [
    "tai_start_score",
    "tai_stop_score",
    "orf_type",
    "blast_pid",
    "blast_offset",
    "blast_percentage_aligned",
    "rna_score",
    "tai_mean_score",
    "blast_match",
    "log_orf_len",
    "has_canonical_start",
    "has_canonical_stop",
    "netstart_atg_score",
    "transaid_start_score",
    "transaid_stop_score",
    "transaid_integrated_score",
    # "transaid_mean_score",
    "is_missing_rna_score",
    "is_missing_tai_start_score",
    "is_missing_tai_stop_score",
    "is_missing_tai_mean_score",
    "is_missing_netstart_atg_score",
    "is_missing_transaid_start_score",
    "is_missing_transaid_stop_score",
    "is_missing_transaid_mean_score",
    "is_missing_transaid_integrated_score",
    "is_nmd_target",
    "start_ai_consensus",
    "stop_ai_consensus",
    "start_stop_symmetry",
    "canonical_both",
    "tool_coverage",
    "is_sequence_flawed",
    "orf_quality_score",
    "is_high_quality_orf",
    "length_nmd_interaction",
    "ai_quality_interaction",
]
ORF_TYPE_MAPPING: Dict = {
    "CO": 1,
    "CN": 2,
    "UN": 3,
    "FP": 4,
    "TP": 5,
    "TN": 6,
    "FN": 7,
}
NMD_TYPE_MAPPING: Dict = {
    "NN": 1,
    "WN": 2,
    "SN": 3,
    "UN": 4,
}

log = logging.getLogger(__name__)
logging.basicConfig(encoding="utf-8", level=logging.INFO)


def normalize_missing_sentinels(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert legacy -1 placeholders into NaN so missingness indicators are correct.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to normalize missing sentinels in.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with missing sentinels normalized.

    Example
    -------
    >>> df = pd.DataFrame({
    >>>     'a': [1, 2, 3],
    >>>     'b': [-1, -1, -1]
    >>> })
    >>> df = normalize_missing_sentinels(df)
    >>> print(df.head())
    """
    for col in MISSINGNESS_COLS:
        if col in df.columns:
            df.loc[df[col] == -1, col] = np.nan
    return df


def fill_transaid_mean_from_components(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build transaid_mean_score from start/stop when not present.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to fill transaid_mean_score from components.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with transaid_mean_score filled from components.

    Example
    -------
    >>> df = pd.DataFrame({
    >>>     'transaid_start_score': [1, 2, 3],
    >>>     'transaid_stop_score': [4, 5, 6],
    >>>     'transaid_mean_score': [7, 8, 9]
    >>> })
    >>> df = fill_transaid_mean_from_components(df)
    >>> print(df.head())
    """
    if "transaid_mean_score" not in df.columns:
        df["transaid_mean_score"] = np.nan

    required = {"transaid_start_score", "transaid_stop_score", "transaid_mean_score"}
    if required.issubset(df.columns):
        can_derive = (
            df["transaid_mean_score"].isna()
            & df["transaid_start_score"].notna()
            & df["transaid_stop_score"].notna()
        )
        df.loc[can_derive, "transaid_mean_score"] = (
            df.loc[can_derive, "transaid_start_score"]
            + df.loc[can_derive, "transaid_stop_score"]
        ) / 2.0
    return df


def normalize_discrete_code_column(
    series: pd.Series,
    mapping: Dict[str, int],
    column_name: str,
) -> pd.Series:
    """
    Normalize a discrete code column by mapping values to integers.

    Parameters
    ----------
    series : pd.Series
        The series to normalize.
    mapping : Dict[str, int]
        The mapping from original values to integers.
    column_name : str
        The name of the column being normalized.

    Returns
    -------
    pd.Series
        The normalized series.

    Example
    -------
    >>> series = pd.Series(['A', 'B', 'C', 'D'])
    >>> mapping = {'A': 1, 'B': 2, 'C': 3, 'D': 4}
    >>> normalized_series = normalize_discrete_code_column(series, mapping, 'column_name')
    >>> print(normalized_series)
    """
    numeric = pd.to_numeric(series, errors="coerce")
    normalized = numeric.copy()

    unresolved = normalized.isna()
    if unresolved.any():
        labels = series.astype(str).str.strip().str.upper()
        mapped = labels.map(mapping)
        normalized.loc[unresolved] = mapped.loc[unresolved]

    remaining = normalized.isna()
    if remaining.any():
        unknown = (
            series.loc[remaining]
            .astype(str)
            .str.strip()
            .replace("", "<empty>")
            .unique()
            .tolist()
        )
        log.warning(
            "WARNING: %s has %d unmapped values; examples: %s. They remain NaN.",
            column_name,
            int(remaining.sum()),
            unknown[:8],
        )

    return normalized


def run(args: argparse.Namespace) -> None:
    """
    Executes the main workflow for loading a model, performing predictions,
    and saving the results.

    This function orchestrates the process of loading a RandomForestClassifier
    model, reading and preparing input data from BLAST, TAI, and TOGA files,
    generating predictions, and then saving the final table to a TSV file.
    It also includes a commented-out section for an 'overrule' function,
    indicating potential future functionality.

    Parameters
    ----------
    args : argparse.Namespace
        An object containing command-line arguments, expected to have attributes
        `model` (path to model), `blast` (path to BLAST file),
        `tai` (path to TAI file), `toga` (path to TOGA file), and
        `outdir` (output directory).

    Returns
    -------
    None

    Example
    -------
    >>> import argparse
    >>> # Create a dummy Namespace for demonstration
    >>> # dummy_args = argparse.Namespace(
    >>> #     model="my_model.joblib",
    >>> #     blast="blast.tsv",
    >>> #     tai="tai.tsv",
    >>> #     toga="toga.tsv",
    >>> #     outdir="output_data",
    >>> #     # overrule=False # Include if overrule is an argument
    >>> # )
    >>> # # Assuming necessary files exist and a model can be loaded
    >>> # run(dummy_args)
    >>> # # This would create a 'predictions.tsv' file in 'output_data'
    """
    model = xgb.XGBClassifier()
    model.load_model(args.model)

    log.info("INFO: Model loaded successfully!")
    log.info(f"INFO: Number of features: {model.n_features_in_}")

    table = predict(args.blast, args.samba, model, threshold=args.threshold)

    if len(table) == 0:
        log.info(f"INFO: No predictions found in {args.blast}/{args.samba}, will not map to blocks!")
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        log.info(f"INFO: Creating empty {outdir}/{args.prefix}.predictions.[tsv|bed]")
        Path(f"{outdir}/{args.prefix}.predictions.tsv").touch()
        Path(f"{outdir}/{args.prefix}.predictions.bed").touch()

        return

    log.info(f"INFO: Predictions found: {len(table)}")

    _ = map_to_blocks(
        table,
        args.alignments,
        args.outdir,
        args.prefix,
        args.min_score_max_predictions,
        args.max_predictions,
        keep_raw=args.keep_raw,
    )


def map_to_blocks(
    table: pd.DataFrame,
    alignments: Union[str, PathLike],
    outdir: Union[str, PathLike],
    prefix: str,
    min_score_max_predictions: float,
    max_predictions: int,
    keep_raw: bool = False,
) -> None:
    """
    Maps prediction results to genomic blocks based on alignment data and saves them.

    This function takes a DataFrame of prediction results, merges it with genomic
    alignment data (BED file format), modifies certain columns to represent block
    coordinates, and then saves both the updated prediction table and the new
    BED-formatted table to the specified output directory.

    Parameters
    ----------
    table : pd.DataFrame
        A pandas DataFrame containing prediction results, expected to have
        'blast_id', 'tai_orf_start', and 'tai_orf_end' columns.
    alignments : Union[str, PathLike]
        Path to the genomic alignments file (BED format).
    outdir : Union[str, PathLike]
        The output directory where the modified prediction table and
        the new BED file will be saved.
    prefix : str
        A string prefix to be added to the output filenames.

    Returns
    -------
    None

    Example
    -------
    >>> import pandas as pd
    >>> # Create dummy dataframes and paths for demonstration
    >>> # dummy_table = pd.DataFrame({
    >>> #     'blast_id': ['gene1__orfA', 'gene2__orfB'],
    >>> #     'tai_orf_start': [100, 500],
    >>> #     'tai_orf_end': [200, 600]
    >>> # })
    >>> # with open("dummy_alignments.bed", "w") as f:
    >>> #     f.write("chr1\\t10\\t1000\\tgene1__orfA\\t0\\t+\\n")
    >>> #     f.write("chr1\\t50\\t1200\\tgene2__orfB\\t0\\t+\\n")
    >>> #
    >>> # map_to_blocks(dummy_table, "dummy_alignments.bed", "./output", "test_")
    >>> # # This would create 'test_predictions.tsv' and 'test_predictions.bed'
    >>> # # in the './output' directory.
    """
    log.info(f"INFO: Initial size of table: {len(table)}")
    log.info(f"INFO: Table looks like this:\n{table[['prefix', 'prob_coding']].head()}")
    table = table.copy()
    # table["id"] = [id.split("__")[0] for id in table.blast_id] -> replace by prefix

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info(f"INFO: Reading alignments from {alignments}")
    bed = pd.read_csv(alignments, sep="\t", header=None)
    log.info(f"INFO: Alignments has {len(bed)} rows!")
    log.info(f"INFO: Alignments looks like this:\n{bed[[0, 3]].head()}")

    merged = _map_predictions_to_bed(bed, table)

    log.info(f"INFO: Merged table has {len(merged)} rows!")
    log.info(f"INFO: Merged table looks like this:\n{merged.head()}")

    if keep_raw:
        merged.drop(columns=["prefix", "start", "end", "id"]).to_csv(
            f"{outdir}/{prefix}.all.predictions.bed",
            sep="\t",
            header=False,
            index=False,
        )

        table.drop(columns=["prefix"]).to_csv(
            f"{outdir}/{prefix}.all.predictions.tsv", index=False, header=True, sep="\t"
        )

    table = (
        table[table["prob_coding"] >= min_score_max_predictions]
        .sort_values("prob_coding", ascending=False, kind="stable")
        .groupby("prefix", sort=False)
        .head(max_predictions)
        .reset_index(drop=True)
    )

    # INFO: add #DU tag to non-best predictions
    table["rank"] = table.groupby("prefix", sort=False).cumcount() + 1
    table.loc[table["rank"] > 1, "id"] += "#DU"

    log.info(f"INFO: Final size of table: {len(table)}")
    merged = _map_predictions_to_bed(bed, table)

    table_ids = set(table["id"])
    bed_ids = set(merged[3])
    if table_ids != bed_ids:
        missing_from_bed = sorted(table_ids - bed_ids)[:5]
        missing_from_tsv = sorted(bed_ids - table_ids)[:5]
        raise ValueError(
            "ERROR: Prediction BED/TSV IDs differ before writing; "
            f"missing from BED: {missing_from_bed}; "
            f"missing from TSV: {missing_from_tsv}"
        )

    log.info(f"INFO: Writing predictions to {outdir}/{prefix}.predictions.tsv")
    table.drop(columns=["prefix", "rank"]).to_csv(
        f"{outdir}/{prefix}.predictions.tsv", index=False, header=True, sep="\t"
    )

    log.info(f"INFO: Writing predictions to {outdir}/{prefix}.predictions.bed")
    merged.drop(columns=["prefix", "start", "end", "prob_coding", "id"]).to_csv(
        f"{outdir}/{prefix}.predictions.bed", sep="\t", header=False, index=False
    )

    return


def _map_predictions_to_bed(
    bed: pd.DataFrame,
    table: pd.DataFrame,
) -> pd.DataFrame:
    """Map prediction coordinates and identifiers onto BED alignment blocks."""
    merged = (
        bed.assign(prefix=bed[3])
        .merge(
            table[["prefix", "start", "end", "prob_coding", "id"]],
            on="prefix",
            how="left",
        )
        .dropna(subset=["start"])
        .assign(
            start=lambda df: df["start"].astype(int),
            end=lambda df: df["end"].astype(int),
        )
    )

    merged[3] = merged["id"]
    merged[6] = merged["start"]
    merged[7] = merged["end"]
    return merged


def predict(
    blast: Union[str, PathLike],
    samba: Union[str, PathLike],
    model: xgb.XGBClassifier,
    threshold: float = 0.5,
) -> pd.DataFrame:
    """
    Performs predictions using a loaded RandomForestClassifier model on input data.

    This function reads data from BLAST, TAI, and TOGA files, prepares a query
    table, and then uses the provided model to generate class probabilities.

    Parameters
    ----------
    blast : Union[str, PathLike]
        Path to the BLAST input file.
    tai : Union[str, PathLike]
        Path to the TAI input file.
    toga : Union[str, PathLike]
        Path to the TOGA input file.
    model : RandomForestClassifier
        The pre-trained RandomForestClassifier model used for prediction.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the merged input data along with two
        new columns: 'class0_prob' and 'class1_prob', which
        represent the prediction probabilities for each class.

    Example
    -------
    >>> # Assuming 'blast.tsv', 'tai.tsv', 'toga.tsv' exist and a model is loaded
    >>> # from sklearn.ensemble import RandomForestClassifier
    >>> # dummy_model = RandomForestClassifier() # Replace with your actual loaded model
    >>> # result_df = predict("blast.tsv", "tai.tsv", "toga.tsv", dummy_model)
    >>> # print(result_df.head())
    """
    # INFO: a chunk with no DIAMOND hits (or no RNAsamba rows) has nothing to score: synthetic
    # INFO: spike-ins (ERCC) and small rescued chunks produce empty tables, which pandas cannot parse
    if any(Path(p).stat().st_size == 0 for p in (blast, samba)):
        log.info(f"INFO: {blast} or {samba} is empty, nothing to predict")
        return pd.DataFrame()

    table = read(blast, samba)
    query = table.loc[:, FEATURES]

    for col in query.columns:
        query[col] = pd.to_numeric(query[col], errors="coerce")

    query = add_engineered_features(query)[ORDER]

    probabilities = model.predict_proba(query)
    predictions = (probabilities[:, 1] >= threshold).astype(int)
    log.info("INFO: Finished predictions!")

    table["predicted_class"] = predictions
    table["prob_noncoding"] = probabilities[:, 0]
    table["prob_coding"] = probabilities[:, 1]

    log.info("INFO: Predictions summary:")
    log.info(f"INFO:  Predicted as real ORFs (1): {(predictions == 1).sum()}")
    log.info(f"INFO:  Predicted as spurious (0): {(predictions == 0).sum()}")

    return table


def read(
    blast: Union[str, PathLike, Path],
    samba: Union[str, PathLike, Path],
) -> pd.DataFrame:
    """
    Reads and merges data from BLAST, TAI, and TOGA files into a single DataFrame.

    Parameters
    ----------
    blast : Union[str, PathLike, Path]
        Path to the BLAST input file.
    samba : Union[str, PathLike, Path]
        Path to the RNASamba input file.

    Returns
    -------
    pd.DataFrame
        A merged pandas DataFrame containing data from all three input files.

    Example
    -------
    >>> # Assuming 'blast.tsv', 'samba.tsv' exist
    >>> # merged_data = read("blast.tsv", "samba.tsv")
    """
    blast = read_blast(blast)
    samba = read_samba(samba)

    table = merge_tables(blast, samba)
    return table


def read_blast(path: Union[str, PathLike, Path]) -> pd.DataFrame:
    """
    Reads a BLAST tab-separated file and processes it into a pandas DataFrame.

    This function reads the file using predefined column names (BLAST_COLS)
    and generates a unique 'key' column based on 'blast_id', 'blast_orf_start',
    and 'blast_orf_end'.

    Parameters
    ----------
    path : Union[str, PathLike, Path]
        The file path to the BLAST input file.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the BLAST data with an additional 'key' column.

    Example
    -------
    >>> # Assuming 'my_blast.tsv' exists with appropriate tab-separated data
    >>> # blast_df = read_blast("my_blast.tsv")
    >>> # print(blast_df.head())
    """
    blast = pd.read_csv(path, sep="\t", header=None, names=BLAST_COLS)

    # INFO: needs to be the canonical ID
    # R1_chr1__OR2#NE1 -> R1_chr1 [samba] + R1_chr1:1-10(+) [bed]
    blast["prefix"] = blast["id"].astype(str).str.split("_ORF").str[0]
    blast["prefix"] = blast["prefix"].str.split("\.p").str[0]

    blast["key"] = (
        blast["prefix"]
        + ":"
        + blast["start"].astype(str)
        + "-"
        + blast["end"].astype(str)
        + "("
        + blast["strand"].astype(str)
        + ")"
    )
    return blast


def read_samba(path: Union[str, PathLike, Path]) -> pd.DataFrame:
    """
    Reads a RNASamba tab-separated file and processes it into a pandas DataFrame.

    This function reads the file using predefined column names (SAMBA_COLS)

    Parameters
    ----------
    path : Union[str, PathLike, Path]
        The file path to the RNASamba input file.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the RNASamba data with an additional 'key' column.

    Example
    >>> # Assuming 'my_samba.tsv' exists with appropriate tab-separated data
    >>> # samba_df = read_samba("my_samba.tsv")
    >>> # print(samba_df.head())
    """
    df = pd.read_csv(path, sep="\t", usecols=[0, 1])
    df.columns = SAMBA_COLS
    return df


def merge_tables(
    blast: pd.DataFrame,
    samba: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merges the BLAST and RNASamba DataFrames into a final comprehensive table.

    The merging process involves:
    1. Inner merging BLAST and RNASamba DataFrames on the 'prefix' column.
    2. Extracting chromosome information into 'm_chr' from the merged key.
    3. Creating a 'toga_key' based on chromosome and ORF start/end, considering strand.
    4. Left merging with the TOGA DataFrame on 'toga_key' and TOGA's 'key'.
    5. Adding a binary flag 'toga_overlap_bp' indicating if 'toga_pid' is present.

    Parameters
    ----------
    blast : pd.DataFrame
        The DataFrame containing BLAST data.
    samba : pd.DataFrame
        The DataFrame containing RNASamba data.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame that is the result of merging all three input DataFrames,
        with additional derived columns like 'm_chr', 'toga_key', and 'toga_overlap_bp'.

    Example
    -------
    >>> # Assuming blast_df, samba_df are pre-read DataFrames
    >>> # final_table = merge_tables(blast_df, samba_df)
    >>> # print(final_table.head())
    """
    merged = blast.merge(samba, on="prefix", how="left")

    log.info(f"INFO: Merged df looks like this:\n{merged.head()}")
    log.info(f"INFO: BLAST has {len(blast.prefix.unique())} unique rows!")
    log.info(f"INFO: RNASamba has {len(samba)} rows!")
    log.info(f"INFO: BLAST and RNASamba files have {len(merged)} rows!")

    did_not_get_any_orf = samba[~samba.prefix.isin(blast.prefix)].sort_values(
        "rna_score", ascending=False
    )
    if did_not_get_any_orf.shape[0] > 0:
        log.info(
            f"INFO: These prefixes did not get any ORF prediction:\n{did_not_get_any_orf}"
        )

    not_catched_by_samba = blast[~blast.prefix.isin(samba.prefix)]
    if not_catched_by_samba.shape[0] > 0:
        log.info(
            f"INFO: These prefixes are not caught by samba: {not_catched_by_samba}"
        )

    if len(merged.prefix.unique()) != len(blast.prefix.unique()):
        raise ValueError(
            "ERROR: BLAST and RNASamba files do not match! Some BLAST rows do not have a prediction in RNASamba!"
        )

    merged = normalize_missing_sentinels(merged)
    merged["tai_mean_score"] = (merged.tai_start_score + merged.tai_stop_score) / 2
    merged["blast_match"] = [1 if x > 0 else 0 for x in merged.blast_percentage_aligned]
    merged["log_orf_len"] = np.log1p(
        merged.relative_orf_end - merged.relative_orf_start
    )
    merged["has_canonical_start"] = merged.start_codon.isin(START_CODONS).astype(int)
    merged["has_canonical_stop"] = merged.stop_codon.isin(STOP_CODONS).astype(int)
    merged["neg_log10_blast_evalue"] = -np.log10(merged.blast_evalue)

    mask = (merged.tai_start_score == 0.0) & (merged.tai_stop_score == 0.0)
    merged.loc[
        mask,
        TAI_MASKING_NAN_COLS,
    ] = np.nan

    mask = (
        (merged["blast_evalue"] == 1)
        | (merged["neg_log10_blast_evalue"] < 0)
        | pd.isna(merged["blast_percentage_aligned"])
    )
    merged.loc[
        mask,
        BLAST_MASKING_NAN_COLS,
    ] = np.nan

    merged["orf_type"] = normalize_discrete_code_column(
        merged["orf_type"], ORF_TYPE_MAPPING, "orf_type"
    )
    merged["nmd_type"] = normalize_discrete_code_column(
        merged["nmd_type"], NMD_TYPE_MAPPING, "nmd_type"
    )
    merged = fill_transaid_mean_from_components(merged)

    return merged


def add_missing_indicators(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds indicator columns for missing values in the input DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to add indicator columns to.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with indicator columns added.

    Example
    -------
    >>> df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
    >>> df = add_missing_indicators(df)
    >>> print(df.head())
    """
    for col in MISSINGNESS_COLS:
        if col in df.columns:
            df[f"is_missing_{col}"] = df[col].isna().astype(int)
    return df


def add_engineered_features(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Adds engineered features to the input DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to add engineered features to.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with engineered features added.

    Example
    -------
    >>> df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
    >>> df = add_engineered_features(df)
    >>> print(df.head())
    """
    df = df.copy()
    df = normalize_missing_sentinels(df)
    df = fill_transaid_mean_from_components(df)
    df = add_missing_indicators(df)

    df["blast_pid"] = df["blast_pid"]
    df["blast_offset"] = df["blast_offset"]
    df["blast_percentage_aligned"] = df["blast_percentage_aligned"]

    # Fill NaNs with sentinel -1 for numeric model consumption
    df[MISSINGNESS_COLS] = df[MISSINGNESS_COLS].fillna(-1)

    # INFO: NMD interactions + additional marker
    df["is_nmd_target"] = np.where(df["nmd_type"] == 1, 0, 1).astype(int)

    conditions = [
        # Both present
        (df["is_missing_tai_mean_score"] == 0)
        & (df["is_missing_transaid_mean_score"] == 0),
        # Only tai present
        (df["is_missing_tai_mean_score"] == 0)
        & (df["is_missing_transaid_mean_score"] == 1),
        # Only transaid present
        (df["is_missing_tai_mean_score"] == 1)
        & (df["is_missing_transaid_mean_score"] == 0),
    ]
    choices = [
        df[["tai_start_score", "transaid_start_score", "netstart_atg_score"]].mean(
            axis=1
        ),
        df[["tai_start_score", "netstart_atg_score"]].mean(axis=1),
        df[["transaid_start_score", "netstart_atg_score"]].mean(axis=1),
    ]
    df["start_ai_consensus"] = np.select(conditions, choices, default=0)

    choices = [
        df[["tai_stop_score", "transaid_stop_score"]].mean(axis=1),
        df[["tai_stop_score"]].mean(axis=1),
        df[["transaid_stop_score"]].mean(axis=1),
    ]
    df["stop_ai_consensus"] = np.select(conditions, choices, default=0)

    df["start_stop_symmetry"] = np.where(
        (df["start_ai_consensus"] != 0) & (df["stop_ai_consensus"] != 0),
        1 - np.abs(df["start_ai_consensus"] - df["stop_ai_consensus"]),
        -1,
    )

    df["canonical_both"] = df["has_canonical_start"] * df["has_canonical_stop"]

    df["tool_coverage"] = (
        df["blast_match"]
        + (1 - df["is_missing_tai_mean_score"])
        + (1 - df["is_missing_rna_score"])
        + (1 - df["is_missing_transaid_mean_score"])
        + (1 - df["is_missing_netstart_atg_score"])
    )

    # INFO: is inner_stops == 1, else sequence is flawed
    df["is_sequence_flawed"] = np.where(df["inner_stops"] == 1, 0, 1).astype(int)

    # INFO: trying to create a variable that weights longer-uninterrupted-cannonical-non-NMD
    # ORFs by using 'log_orf_len', 'cannonical_both', 'is_nmd_target', 'is_flawed_sequence'
    min_log = df["log_orf_len"].min()
    max_log = df["log_orf_len"].max()
    print(f"ORF length range: {min_log:.2f}-{max_log:.2f}")
    denom = max(max_log - min_log, 1e-8)
    normalized_log_orf_len = (df["log_orf_len"] - min_log) / denom

    # gives 10% weight to length, and 30%/15%/15% to the other features
    df["orf_quality_score"] = (
        0.2 * normalized_log_orf_len
        + 0.05 * df["canonical_both"]
        + 0.25 * (1 - df["is_missing_tai_mean_score"])
        + 0.20 * (1 - df["is_nmd_target"])
        + 0.40 * (1 - df["is_sequence_flawed"])
    )

    df["is_high_quality_orf"] = np.where(df["orf_quality_score"] > 0.6, 1, 0)

    df["length_nmd_interaction"] = 1 - (df["is_nmd_target"] * normalized_log_orf_len)
    df["ai_quality_interaction"] = df[
        ["rna_score", "start_ai_consensus", "stop_ai_consensus", "orf_quality_score"]
    ].mean(axis=1)

    # Select final feature set
    exclude = {
        "label",
        "negative_type",
        "species",
        "id",
        "transcript_id",
        "chr",
        "start",
        "end",
        "strand",
        "start_codon",
        "stop_codon",
        "key",
        "blast_evalue",
        "blast_length",
        "inner_stops",  # WARN: replaced by is_sequence_flawed
        "nmd_type",  # WARN: replaced by is_nmd_target
        "neg_log10_blast_evalue",
        "transaid_mean_score",
    }
    feature_cols = [c for c in df.columns if c not in exclude]

    return df[feature_cols]


def parse() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    argparse.Namespace

    Example
    -------
    >>> parse()
    """
    parser = argparse.ArgumentParser(
        description="Run Random Forest classifier to predict open-reading-frames"
    )
    parser.add_argument(
        "-b", "--blast", required=True, help="Path to BLAST results file"
    )
    parser.add_argument(
        "-s", "--samba", required=True, help="Path to RNASamba results file"
    )
    parser.add_argument(
        "-a", "--alignments", required=True, help="Path to aligned .bed file"
    )
    parser.add_argument(
        "-m",
        "--model",
        type=str,
        help="Path to ORF model",
        default=MODEL_PATH,
    )
    parser.add_argument(
        "-T",
        "--threshold",
        type=float,
        default=0.5,
        help="Use a non-default threshold for classification",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to outdir",
        default=".",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="Prefix for the output file name",
        default="xgb",
    )
    parser.add_argument(
        "-mm",
        "--min-score-max-predictions",
        type=float,
        default=0.50,
        help="Minimum score to keep a prediction(s), controlled by max_predictions [default: 0.70]",
    )
    parser.add_argument(
        "-mp",
        "--max-predictions",
        type=int,
        default=1,
        help="Maximum number of predictions to keep per query [default: 1]",
    )
    parser.add_argument(
        "-K", "--keep-raw", action="store_true", help="Keep raw predictions"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Version: {__version__}",
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse()
    run(args)
