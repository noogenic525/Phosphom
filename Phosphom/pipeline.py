"""
Pipeline orchestrator — runs all five steps in sequence.

Steps
-----
1. **Mapping**       – Identify kinase motifs in protein sequences (Word doc)
2. **Extraction**    – Extract structured data from the mapped document
3. **Normalisation** – One kinase per row; build cell-signal tables
4. **Validation**    – Compare predictions against PSP reference
5. **Scoring**       – Accuracy, precision, recall, F1-score

Metadata (version, timestamp, input files, parameters) is recorded in the
output Excel ``metrics`` sheet for reproducibility.
"""

from __future__ import annotations

import datetime
import logging
import os
from typing import Any, Dict, Optional

import pandas as pd

from . import __version__
from .extractor import extract_data_from_docx
from .mapping import start_mapping
from .normalization import build_cell_signal_tables, normalize_kinase_rows
from .validation import calculate_f1_metrics, validate_predictions

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Bundled reference file path
# ---------------------------------------------------------------------------
_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_REF_PATH = os.path.join(_PACKAGE_DIR, "Substrates of protein.xlsx")


def get_default_ref_path() -> str:
    """Return the path to the bundled PhosphoSitePlus reference file.

    Raises
    ------
    FileNotFoundError
        If the bundled reference file is missing.
    """
    if not os.path.isfile(DEFAULT_REF_PATH):
        raise FileNotFoundError(
            f"Bundled reference file not found: {DEFAULT_REF_PATH}\n"
            "Please provide a reference file manually via -r / --ref."
        )
    return DEFAULT_REF_PATH


def run_pipeline(
    word_path: str,
    ref_path: Optional[str] = None,
    output_dir: str = ".",
    base_name: str = "result",
    min_confidence: float = 0.0,
) -> Dict[str, Any]:
    """Execute the full Phosphom analysis pipeline.

    Parameters
    ----------
    word_path : str
        Input Word document with protein sequences and peptide fragments.
    ref_path : str, optional
        Reference Excel file (PhosphoSitePlus substrate data).
        If ``None``, the bundled ``Substrates of protein.xlsx`` is used
        automatically.
    output_dir : str
        Directory for result files (default ``"."``).
    base_name : str
        Base filename prefix for outputs (default ``"result"``).
    min_confidence : float
        Minimum motif-matching confidence score (default 0.0).

    Returns
    -------
    dict
        Keys: ``mapped_docx``, ``csv_path``, ``excel_path``, ``processed``,
        ``acc_summary``, ``f1_metrics``, ``metadata``.
    """
    # Resolve reference file: use bundled default if not provided
    if ref_path is None:
        ref_path = get_default_ref_path()
        logger.info("Using bundled reference file: %s", ref_path)

    os.makedirs(output_dir, exist_ok=True)

    mapped_docx = os.path.join(output_dir, f"{base_name}_mapped.docx")
    csv_path = os.path.join(output_dir, f"{base_name}_normalized.csv")
    excel_path = os.path.join(output_dir, f"{base_name}_validated.xlsx")

    # ── Step 1: Mapping ───────────────────────────────────────────────────
    logger.info("Step 1/5: Mapping phosphorylation sites …")
    processed = start_mapping(word_path, mapped_docx, min_confidence)

    # ── Step 2: Extraction ────────────────────────────────────────────────
    logger.info("Step 2/5: Extracting structured data …")
    df_raw = extract_data_from_docx(mapped_docx)
    if df_raw.empty:
        raise ValueError("No motif data matching the criteria was found.")

    # ── Step 3: Normalisation ─────────────────────────────────────────────
    logger.info("Step 3/5: Normalising kinase data …")
    df_norm = normalize_kinase_rows(df_raw)
    df_norm.to_csv(csv_path, index=False)

    # ── Step 4: Validation ────────────────────────────────────────────────
    logger.info("Step 4/5: Validating against PSP reference …")
    final_df, acc_summary = validate_predictions(df_norm, ref_path)

    # ── Step 5: Scoring ───────────────────────────────────────────────────
    logger.info("Step 5/5: Calculating F1 metrics …")
    f1_metrics = calculate_f1_metrics(final_df)

    # ── Cell-signal analysis ──────────────────────────────────────────────
    gene_signal_df, overall_signal_df, top_kinases_df = build_cell_signal_tables(
        df_norm
    )

    # ── Metadata for reproducibility ──────────────────────────────────────
    metadata: Dict[str, Any] = {
        "phosphom_version": __version__,
        "timestamp": datetime.datetime.now().isoformat(),
        "input_file": os.path.basename(word_path),
        "reference_file": os.path.basename(ref_path),
        "min_confidence": min_confidence,
    }

    metrics_df = pd.DataFrame(
        [
            {
                **acc_summary,
                **f1_metrics,
                "mapped_items": processed,
                **metadata,
            }
        ]
    )

    # ── Write Excel output ────────────────────────────────────────────────
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        final_df.to_excel(writer, index=False, sheet_name="validated_data")
        metrics_df.to_excel(writer, index=False, sheet_name="metrics")

        sig_sheet = "cell_signal"
        gene_signal_df.to_excel(
            writer, index=False, sheet_name=sig_sheet, startrow=0
        )
        ov_start = len(gene_signal_df) + 2
        overall_signal_df.to_excel(
            writer, index=False, sheet_name=sig_sheet, startrow=ov_start
        )
        top_start = ov_start + len(overall_signal_df) + 2
        top_kinases_df.to_excel(
            writer, index=False, sheet_name=sig_sheet, startrow=top_start
        )

    logger.info("Pipeline complete. Results saved to: %s", output_dir)

    return {
        "mapped_docx": mapped_docx,
        "csv_path": csv_path,
        "excel_path": excel_path,
        "processed": processed,
        "acc_summary": acc_summary,
        "f1_metrics": f1_metrics,
        "metadata": metadata,
    }
