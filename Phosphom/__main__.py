"""
CLI & GUI entry point for the Phosphom package.

Usage
-----
    # Launch GUI
    python -m Phosphom --gui

    # Run pipeline via CLI
    python -m Phosphom -w input.docx -r reference.xlsx -o ./results -n analysis

    # With minimum confidence filter and verbose logging
    python -m Phosphom -w input.docx -r ref.xlsx -c 0.5 -v
"""

from __future__ import annotations

import argparse
import logging
import sys

from . import __version__


def main() -> None:
    """Parse arguments and dispatch to CLI pipeline or GUI."""
    parser = argparse.ArgumentParser(
        prog="Phosphom",
        description="Phosphom — Kinase Phosphorylation Site Analysis Tool",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"Phosphom {__version__}",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Launch GUI mode",
    )
    parser.add_argument(
        "-w",
        "--word",
        type=str,
        help="Path to input Word (.docx) file",
    )
    parser.add_argument(
        "-r",
        "--ref",
        type=str,
        default=None,
        help="Path to reference Excel (.xlsx) file (PhosphoSitePlus). "
        "If omitted, the bundled 'Substrates of protein.xlsx' is used.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=".",
        help="Output directory for results (default: current directory)",
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        default="result",
        help="Base name for output files (default: result)",
    )
    parser.add_argument(
        "-c",
        "--confidence",
        type=float,
        default=0.0,
        help="Minimum confidence score (0.0–1.0, default: 0.0)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose logging output (DEBUG level)",
    )
    args = parser.parse_args()

    # Configure logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # ── GUI mode ──────────────────────────────────────────────────────────
    if args.gui:
        from .gui import build_gui

        build_gui()
        return

    # ── CLI mode ──────────────────────────────────────────────────────────
    if not args.word:
        parser.error(
            "--word (-w) is required.  "
            "Use --gui for GUI mode."
        )

    from .pipeline import run_pipeline

    try:
        result = run_pipeline(
            word_path=args.word,
            ref_path=args.ref,
            output_dir=args.output,
            base_name=args.name,
            min_confidence=args.confidence,
        )

        acc = result["acc_summary"]
        f1 = result["f1_metrics"]

        print(f"\n{'=' * 55}")
        print("  Phosphom Analysis Complete")
        print(f"{'=' * 55}")
        print(f"  Mapped items : {result['processed']}")
        print(f"  CSV output   : {result['csv_path']}")
        print(f"  Excel output : {result['excel_path']}")
        print(f"{'─' * 55}")
        print(f"  ACC          : {acc['acc_percent']:.2f} %")
        print(f"  Precision    : {f1['precision']:.4f}")
        print(f"  Recall       : {f1['recall']:.4f}")
        print(f"  F1-score     : {f1['f1_score']:.4f}")
        print(f"  TP / FP / FN : {f1['tp']} / {f1['fp']} / {f1['fn']}")
        print(f"{'=' * 55}\n")

    except Exception as e:
        logging.getLogger(__name__).error("Pipeline failed: %s", e, exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
