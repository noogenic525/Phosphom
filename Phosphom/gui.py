"""
Graphical user interface (tkinter) for the Phosphom pipeline.

Launch via ``python -m Phosphom --gui`` or by calling ``build_gui()``
directly.
"""

from __future__ import annotations

import logging
import os
import threading
import traceback
import tkinter as tk
from tkinter import filedialog, messagebox

from . import __version__
from .pipeline import run_pipeline

logger = logging.getLogger(__name__)


def build_gui() -> None:
    """Build and run the Phosphom tkinter GUI."""
    root = tk.Tk()
    root.title("Phosphom: Kinase Analysis Tool v2.0")
    root.geometry("680x420")
    root.resizable(True, True)

    word_var = tk.StringVar()
    out_var = tk.StringVar()
    base_var = tk.StringVar(value="result")
    conf_var = tk.StringVar(value="0.0")

    # ── File selection callbacks ──────────────────────────────────────────
    def select_word():
        path = filedialog.askopenfilename(
            title="Select Word file",
            filetypes=[("Word files", "*.docx")],
        )
        if path:
            word_var.set(path)
            if not out_var.get():
                out_var.set(os.path.dirname(path))

    def select_out():
        path = filedialog.askdirectory(title="Select output folder")
        if path:
            out_var.set(path)

    # ── Run callback ──────────────────────────────────────────────────────
    def run():
        word_path = word_var.get().strip()
        output_dir = out_var.get().strip()
        base_name = base_var.get().strip() or "result"

        try:
            min_conf = float(conf_var.get().strip())
        except ValueError:
            min_conf = 0.0

        if not word_path or not output_dir:
            messagebox.showwarning(
                "Warning",
                "Please select Word file and output folder.",
            )
            return

        run_btn.config(state="disabled", text="Analyzing...")

        def _task():
            try:
                result = run_pipeline(
                    word_path=word_path,
                    ref_path=None,  # auto-use bundled reference
                    output_dir=output_dir,
                    base_name=base_name,
                    min_confidence=min_conf,
                )
                root.after(0, lambda: _on_success(result))
            except Exception as e:
                logger.error("Pipeline failed:\n%s", traceback.format_exc())
                root.after(0, lambda exc=e: _on_error(exc))

        def _on_success(result):
            run_btn.config(state="normal", text="Start Analysis")
            summary = (
                f"Success!\n\n"
                f"  Mapped items : {result['processed']}\n"
                f"  CSV output   : {os.path.basename(result['csv_path'])}\n"
                f"  Excel output : {os.path.basename(result['excel_path'])}\n\n"
                f"  ACC          : {result['acc_summary']['acc_percent']:.2f} %\n"
                f"  Precision    : {result['f1_metrics']['precision']:.4f}\n"
                f"  Recall       : {result['f1_metrics']['recall']:.4f}\n"
                f"  F1-score     : {result['f1_metrics']['f1_score']:.4f}"
            )
            messagebox.showinfo("Success", summary)

        def _on_error(exc):
            run_btn.config(state="normal", text="Start Analysis")
            messagebox.showerror("Error", f"Pipeline failed:\n{exc}")

        threading.Thread(target=_task, daemon=True).start()

    # ── Layout ────────────────────────────────────────────────────────────
    tk.Label(
        root,
        text="Phosphom: Kinase Analysis Tool",
        font=("Arial", 16, "bold"),
    ).pack(pady=10)

    frame = tk.Frame(root)
    frame.pack(fill="x", padx=20, pady=5)
    frame.columnconfigure(1, weight=1)

    file_rows = [
        ("Word file", word_var, select_word),
        ("Output folder", out_var, select_out),
    ]

    for i, (label, var, cmd) in enumerate(file_rows):
        tk.Label(frame, text=label).grid(row=i, column=0, sticky="w", pady=3)
        tk.Entry(frame, textvariable=var, width=50).grid(
            row=i, column=1, padx=5, sticky="we"
        )
        tk.Button(frame, text="Select", command=cmd, width=6).grid(
            row=i, column=2, padx=2
        )

    row_idx = len(file_rows)

    tk.Label(frame, text="Result name").grid(row=row_idx, column=0, sticky="w", pady=3)
    tk.Entry(frame, textvariable=base_var, width=20).grid(
        row=row_idx, column=1, sticky="w", padx=5
    )

    row_idx += 1
    tk.Label(frame, text="Min confidence").grid(
        row=row_idx, column=0, sticky="w", pady=3
    )
    conf_frame = tk.Frame(frame)
    conf_frame.grid(row=row_idx, column=1, sticky="w", padx=5)
    tk.Entry(conf_frame, textvariable=conf_var, width=8).pack(side="left")
    tk.Label(conf_frame, text="  (0.0 ~ 1.0, default: 0.0)", fg="gray").pack(
        side="left"
    )

    # Start button
    run_btn = tk.Button(
        root,
        text="Start Analysis",
        command=run,
        bg="#2563EB",
        fg="white",
        activebackground="#1D4ED8",
        activeforeground="white",
        width=20,
        height=2,
        font=("Arial", 12, "bold"),
    )
    run_btn.pack(pady=18)

    # Version label
    tk.Label(root, text=f"v{__version__}", fg="gray", font=("Arial", 9)).pack(
        side="bottom", pady=4
    )

    root.mainloop()
