"""
Steps 4–5 — Prediction validation (PSP) and F1-score calculation.

Key improvements over v1
------------------------
* Sheet-name matching uses normalised **exact comparison** to avoid false
  positives from substring matching (e.g. CK1 ≠ CK10).
* Kinase-name comparison uses **exact matching** after normalisation,
  eliminating false positives from substring matching (e.g. CK1 ≠ CK12).
* Reference data is read via ``ExcelFile.parse()`` instead of repeated
  ``pd.read_excel()`` calls, opening the file only once.
* F1 calculation uses the same ``_normalize_name_for_matching`` helper
  for consistent name comparison across the pipeline.
* Substring validation requires minimum sequence length to prevent
  short-sequence false positives.
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Tuple

import pandas as pd

from .motif_db import KINASE_KEYS

logger = logging.getLogger(__name__)

# Minimum sequence length for reliable substring matching in validation
_MIN_MATCH_LEN: int = 7


# ---------------------------------------------------------------------------
# Name normalisation helpers
# ---------------------------------------------------------------------------
def _normalize_name_for_matching(name: str) -> str:
    """Normalise a kinase/sheet name for comparison.

    Strips ``& - _ / spaces`` and uppercases so that e.g.
    ``"p38&MAPK"`` and ``"p38 MAPK"`` compare as equal.
    """
    return re.sub(r"[&\s\-_/]", "", str(name)).upper()


# ---------------------------------------------------------------------------
# Reference data loading
# ---------------------------------------------------------------------------
def load_reference_data(ref_file_path: str) -> Dict[str, List[str]]:
    """Load substrate reference sequences from a PhosphoSitePlus Excel file.

    Each worksheet whose name matches a known kinase keyword is read, and
    the ``SITE_+/-7_AA`` column is extracted as a list of cleaned reference
    sequences for that kinase.

    The file is opened once via ``pd.ExcelFile`` and individual sheets are
    read with ``.parse()`` for efficiency.
    """
    excel_reader = pd.ExcelFile(ref_file_path)
    all_sheets = excel_reader.sheet_names
    kinase_ref_data: Dict[str, List[str]] = {}

    for k_name in KINASE_KEYS:
        target_sheet = None
        norm_k = _normalize_name_for_matching(k_name)
        for s in all_sheets:
            norm_s = _normalize_name_for_matching(s)
            if norm_k == norm_s:
                target_sheet = s
                break

        if target_sheet:
            df_ref = excel_reader.parse(target_sheet)
            if "SITE_+/-7_AA" in df_ref.columns:
                seqs = (
                    df_ref["SITE_+/-7_AA"]
                    .astype(str)
                    .str.replace("_", "", regex=False)
                    .str.upper()
                    .str.strip()
                    .unique()
                    .tolist()
                )
                kinase_ref_data[k_name] = seqs
                logger.debug(
                    "Loaded %d ref sequences for %s (sheet: %s)",
                    len(seqs),
                    k_name,
                    target_sheet,
                )
            else:
                logger.warning(
                    "Sheet '%s' is missing the 'SITE_+/-7_AA' column",
                    target_sheet,
                )

    logger.info(
        "Reference data loaded: %d / %d kinases matched to sheets",
        len(kinase_ref_data),
        len(KINASE_KEYS),
    )
    return kinase_ref_data


# ---------------------------------------------------------------------------
# Kinase-name exact matching
# ---------------------------------------------------------------------------
def _kinase_names_match(pred: str, actual: str) -> bool:
    """Return True if *pred* and *actual* refer to the same kinase.

    Uses normalised exact comparison to avoid false positives from
    substring matching (e.g. ``CK1`` must not match ``CK12``).
    """
    return _normalize_name_for_matching(pred) == _normalize_name_for_matching(actual)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def validate_predictions(
    pm_df: pd.DataFrame,
    ref_file_path: str,
) -> Tuple[pd.DataFrame, Dict]:
    """Validate kinase predictions against PhosphoSitePlus reference data.

    For each motif in *pm_df*, checks whether it appears in any reference
    kinase's substrate list (subsequence search with minimum length guard),
    then verifies whether the predicted kinase matches any found reference
    kinase using **exact** name comparison.

    Returns
    -------
    tuple[pd.DataFrame, dict]
        ``(annotated_df, accuracy_summary)``
    """
    kinase_ref_data = load_reference_data(ref_file_path)
    results: List[dict] = []

    for _, row in pm_df.iterrows():
        p_motif = str(row["Motif"]).upper()
        p_preds_raw = (
            str(row["Kinase Name"]) if pd.notna(row["Kinase Name"]) else ""
        )
        p_preds = [p.strip() for p in p_preds_raw.split(",") if p.strip()]

        # Find all reference kinases whose substrates contain this motif
        # Guard: require minimum sequence length to prevent trivial matches
        found_kinases: List[str] = []
        for k_name, ref_seqs in kinase_ref_data.items():
            if any(
                (len(p_motif) >= _MIN_MATCH_LEN and p_motif in s)
                or (len(s) >= _MIN_MATCH_LEN and s in p_motif)
                for s in ref_seqs
            ):
                found_kinases.append(k_name)

        # Check correctness using exact name matching
        is_correct = False
        if found_kinases:
            for pred in p_preds:
                if is_correct:
                    break
                for actual in found_kinases:
                    if _kinase_names_match(pred, actual):
                        is_correct = True
                        break

        results.append(
            {
                "In_PSP": len(found_kinases) > 0,
                "Correct": is_correct,
                "Actual_Kinases_in_PSP": ", ".join(found_kinases),
            }
        )

    res_df = pd.DataFrame(results)
    final_df = pd.concat([pm_df, res_df], axis=1)

    valid_rows = res_df[res_df["In_PSP"]]
    n_correct = int(valid_rows["Correct"].sum()) if not valid_rows.empty else 0
    acc = (n_correct / len(valid_rows) * 100) if not valid_rows.empty else 0.0

    acc_summary = {
        "total_motifs": len(pm_df),
        "matched_psp": len(valid_rows),
        "correct": n_correct,
        "acc_percent": float(acc),
    }
    logger.info(
        "Validation: %d / %d matched PSP, %d correct (ACC %.2f %%)",
        acc_summary["matched_psp"],
        acc_summary["total_motifs"],
        acc_summary["correct"],
        acc_summary["acc_percent"],
    )
    return final_df, acc_summary


# ---------------------------------------------------------------------------
# F1-score calculation
# ---------------------------------------------------------------------------
def _parse_field(s) -> List[str]:
    """Split a comma/semicolon-separated field into trimmed strings."""
    if pd.isna(s):
        return []
    return [p.strip() for p in re.split(r"[;,]", str(s)) if p.strip()]


def calculate_f1_metrics(df: pd.DataFrame) -> Dict[str, float]:
    """Calculate micro-averaged precision, recall, and F1-score.

    For each unique motif, compares predicted kinases (``Kinase Name``)
    against actual kinases (``Actual_Kinases_in_PSP``) using normalised
    exact matching.

    Returns
    -------
    dict
        Keys: ``precision``, ``recall``, ``f1_score``, ``tp``, ``fp``, ``fn``
    """
    total_tp = 0
    total_fp = 0
    total_fn = 0

    for _, group in df.groupby("Motif"):
        # Actual kinases from PSP
        actual_raw = group["Actual_Kinases_in_PSP"].dropna()
        actual_list = _parse_field(actual_raw.iloc[0]) if not actual_raw.empty else []
        actual_set = {_normalize_name_for_matching(a) for a in actual_list}

        # Predicted kinases
        preds_raw = group["Kinase Name"].dropna().astype(str).tolist()
        pred_set: set = set()
        for item in preds_raw:
            for p in _parse_field(item):
                pred_set.add(_normalize_name_for_matching(p))

        total_tp += len(pred_set & actual_set)
        total_fp += len(pred_set - actual_set)
        total_fn += len(actual_set - pred_set)

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    metrics = {
        "precision": float(precision),
        "recall": float(recall),
        "f1_score": float(f1),
        "tp": int(total_tp),
        "fp": int(total_fp),
        "fn": int(total_fn),
    }
    logger.info(
        "F1 metrics: P=%.4f  R=%.4f  F1=%.4f  (TP=%d  FP=%d  FN=%d)",
        precision,
        recall,
        f1,
        total_tp,
        total_fp,
        total_fn,
    )
    return metrics
