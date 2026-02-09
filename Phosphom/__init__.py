"""
Phosphom — Kinase Phosphorylation Site Analysis Toolkit

A bioinformatics pipeline for:
  1. Mapping phosphorylation sites to kinases using motif pattern matching
     with confidence scoring.
  2. Extracting structured data from annotated Word documents.
  3. Normalising kinase rows and building cell-signaling pathway summaries.
  4. Validating predictions against PhosphoSitePlus (PSP) reference data.
     The bundled ``Substrates of protein.xlsx`` is used automatically.
  5. Computing accuracy, precision, recall, and F1-score.

Usage
-----
CLI (reference file is loaded automatically) ::

    python -m Phosphom -w input.docx -o ./results

CLI (with custom reference file) ::

    python -m Phosphom -w input.docx -r custom_ref.xlsx -o ./results

GUI ::

    python -m Phosphom --gui

Programmatic ::

    from Phosphom import run_pipeline
    # Uses bundled reference automatically
    result = run_pipeline("input.docx", output_dir="./results")
    # Or specify a custom reference
    result = run_pipeline("input.docx", ref_path="custom_ref.xlsx")
"""

__version__ = "2.1.0"
__author__ = "Phosphom Team"

from .pipeline import run_pipeline  # noqa: F401
from .pipeline import get_default_ref_path  # noqa: F401

__all__ = ["run_pipeline", "get_default_ref_path", "__version__"]
