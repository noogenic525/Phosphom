"""
Step 3 — Data normalisation and cell-signaling pathway analysis.

Normalises kinase rows (one kinase per row) and builds summary tables
for cell-signaling pathway analysis.
"""

from __future__ import annotations

import logging
from typing import Tuple

import pandas as pd

from .motif_db import (
    CELL_SIGNAL_CATEGORY_ORDER,
    CELL_SIGNAL_MAP,
    TOP_KINASES_N,
    _normalize_kinase_key,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------
def normalize_kinase_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise so that each row contains a single kinase.

    Splits comma-separated ``Kinase Name`` (and aligned ``Confidence``)
    values into separate rows, then removes duplicates.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``Kinase Name`` and ``Motif`` columns.
        ``Confidence`` column is handled if present.

    Returns
    -------
    pd.DataFrame
        One kinase per row, duplicates removed.
    """
    required = {"Kinase Name", "Motif"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Required columns are missing: {missing}")

    df = df.copy()

    has_conf = "Confidence" in df.columns

    # Split comma-separated kinase names into lists
    df["Kinase Name"] = df["Kinase Name"].fillna("").astype(str).str.split(",")

    if has_conf:
        df["Confidence"] = df["Confidence"].fillna("").astype(str).str.split(",")

        # Ensure same-length lists (pad confidence with "" if shorter)
        def _align(row):
            k_list = row["Kinase Name"]
            c_list = row["Confidence"]
            while len(c_list) < len(k_list):
                c_list.append("")
            return k_list, c_list[: len(k_list)]

        aligned = df.apply(_align, axis=1, result_type="expand")
        df["Kinase Name"] = aligned[0]
        df["Confidence"] = aligned[1]
        df = df.explode(["Kinase Name", "Confidence"])
        df["Confidence"] = df["Confidence"].str.strip()
    else:
        df = df.explode("Kinase Name")

    df["Kinase Name"] = df["Kinase Name"].str.strip()
    df = df[df["Kinase Name"] != ""]
    df = df.drop_duplicates(subset=["Motif", "Kinase Name"]).reset_index(drop=True)

    logger.info("Normalised to %d rows (one kinase per row)", len(df))
    return df


# ---------------------------------------------------------------------------
# Cell-signaling tables
# ---------------------------------------------------------------------------
def build_cell_signal_tables(
    df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build cell-signaling pathway summary tables.

    Parameters
    ----------
    df : pd.DataFrame
        Normalised DataFrame with ``Kinase Name`` and ``Gene Name``.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        ``(gene_summary, overall_summary, top_kinases_per_signal)``
    """
    required = {"Kinase Name", "Gene Name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Required columns are missing: {missing}")

    temp = df.copy()
    temp["Cell Signal"] = temp["Kinase Name"].apply(
        lambda x: CELL_SIGNAL_MAP.get(_normalize_kinase_key(x), "Unknown")
    )
    temp = temp[temp["Cell Signal"] != "Unknown"].copy()
    temp["Cell Signal"] = pd.Categorical(
        temp["Cell Signal"],
        categories=CELL_SIGNAL_CATEGORY_ORDER,
        ordered=True,
    )
    temp = temp.drop_duplicates(subset=["Gene Name", "Motif", "Kinase Name"])

    # ── Gene-level summary ────────────────────────────────────────────────
    grouped = temp.groupby(["Gene Name", "Cell Signal"], dropna=False)
    gene_counts = grouped["Motif"].nunique().reset_index(name="Count")
    gene_motifs = (
        grouped["Motif"]
        .apply(
            lambda x: ", ".join(
                sorted({str(v).strip() for v in x if str(v).strip()})
            )
        )
        .reset_index(name="Motifs")
    )
    gene_kinases = (
        grouped["Kinase Name"]
        .apply(
            lambda x: ", ".join(
                sorted({str(v).strip() for v in x if str(v).strip()})
            )
        )
        .reset_index(name="Kinases")
    )
    gene_summary = (
        gene_counts.merge(gene_motifs, on=["Gene Name", "Cell Signal"])
        .merge(gene_kinases, on=["Gene Name", "Cell Signal"])
        .sort_values(["Gene Name", "Cell Signal"])
        .reset_index(drop=True)
    )

    # ── Overall summary ───────────────────────────────────────────────────
    overall = (
        temp.drop_duplicates(subset=["Motif", "Cell Signal"])
        .groupby("Cell Signal", dropna=False)["Motif"]
        .nunique()
        .rename_axis("Cell Signal")
        .reset_index(name="Count")
    )
    overall["Cell Signal"] = pd.Categorical(
        overall["Cell Signal"],
        categories=CELL_SIGNAL_CATEGORY_ORDER,
        ordered=True,
    )
    overall = overall.sort_values("Cell Signal").reset_index(drop=True)
    total = overall["Count"].sum()
    overall["Percent"] = (overall["Count"] / total * 100).round(2) if total else 0

    # ── Top kinases per signal ────────────────────────────────────────────
    top_base = temp.drop_duplicates(subset=["Motif", "Cell Signal", "Kinase Name"])
    top_counts = (
        top_base.groupby(["Cell Signal", "Kinase Name"], dropna=False)["Motif"]
        .nunique()
        .reset_index(name="Unique Motifs")
    )
    top_counts["Cell Signal"] = pd.Categorical(
        top_counts["Cell Signal"],
        categories=CELL_SIGNAL_CATEGORY_ORDER,
        ordered=True,
    )
    top_counts = top_counts.sort_values(
        ["Cell Signal", "Unique Motifs", "Kinase Name"],
        ascending=[True, False, True],
    )
    top_kinases = (
        top_counts.groupby("Cell Signal", dropna=False)
        .head(TOP_KINASES_N)
        .copy()
        .reset_index(drop=True)
    )
    top_kinases["Rank"] = top_kinases.groupby("Cell Signal").cumcount() + 1
    top_kinases = top_kinases[["Cell Signal", "Rank", "Kinase Name", "Unique Motifs"]]

    logger.info(
        "Cell-signal tables built: %d gene entries, %d signal categories",
        len(gene_summary),
        len(overall),
    )
    return gene_summary, overall, top_kinases
