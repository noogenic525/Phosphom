"""
Step 2 — Structured data extraction from mapped Word documents.

Extracts gene names, motif sequences, kinase predictions, and confidence
scores from the annotated Word documents produced by the mapping step.

Header detection
----------------
Supports multiple formats via NCBI accession patterns:
  - NCBI RefSeq : NP_123456, NM_123456, XP_123456, XM_123456
  - UniProt     : sp|P12345|GENE_HUMAN, tr|Q98765|GENE_MOUSE
  - GenBank     : AAB12345
  - FASTA       : lines starting with '>'
  - Descriptive : mRNA, protein, isoform, variant, complete cds …

Gene symbol extraction priority:
  1. UniProt gene name (``sp|ACC|GENE_SPECIES``)
  2. ``GN=`` field in header
  3. NCBI ``[gene=SYMBOL]`` field
  4. Parenthesised gene symbol ``(GENE)``
  5. First word-like token
"""

from __future__ import annotations

import logging
import re
from typing import List, Tuple

import pandas as pd
from docx import Document

from .motif_db import KINASE_KEYS

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Header-detection patterns (NCBI accession & related)
# ---------------------------------------------------------------------------
_ACCESSION_RE = [
    re.compile(r"(?:sp|tr)\|[A-Z0-9]{6,10}\|"),          # UniProt
    re.compile(r"\b(?:NP|NM|XP|XM|YP|NC)_\d{3,}"),       # NCBI RefSeq
    re.compile(r"\b[A-Z]{3}\d{5}(?:\.\d+)?\b"),           # GenBank protein
    re.compile(r"\bgi\|\d+"),                               # NCBI GI (legacy)
]

_HEADER_KEYWORDS: List[str] = [
    "mrna", "protein", "isoform", "variant", "complete cds", "partial cds",
    "transcript", "homo sapiens", "mus musculus", "rattus norvegicus",
    "gene", "chromosome", "predicted", "uncharacterized",
]


def is_header_line(text: str) -> bool:
    """Detect whether *text* is a gene/protein header line.

    Uses NCBI RefSeq, UniProt, GenBank accession patterns, FASTA '>'
    prefix, and biological descriptor keywords.
    """
    if len(text) < 10:
        return False
    # FASTA header
    if text.startswith(">"):
        return True
    # Accession patterns
    for pat in _ACCESSION_RE:
        if pat.search(text):
            return True
    # Keyword-based detection
    text_lower = text.lower()
    return any(kw in text_lower for kw in _HEADER_KEYWORDS)


def extract_gene_symbol(text: str) -> str:
    """Extract a gene symbol from a header line.

    Attempts (in priority order):
      1. UniProt format ``sp|ACC|GENE_SPECIES``
      2. ``GN=GENE`` field
      3. NCBI ``[gene=SYMBOL]`` field
      4. Parenthesised gene symbol ``(GENE)``
      5. First meaningful token
    """
    # 1. UniProt: sp|P12345|GENE_SPECIES
    m = re.search(r"(?:sp|tr)\|[A-Z0-9]+\|(\w+?)_", text)
    if m:
        return m.group(1)

    # 2. GN= field
    m = re.search(r"GN=(\S+)", text)
    if m:
        return m.group(1)

    # 3. NCBI [gene=SYMBOL]
    m = re.search(r"\[gene=(\S+?)\]", text, re.IGNORECASE)
    if m:
        return m.group(1)

    # 4. Gene symbol in parentheses (not purely numeric)
    m = re.search(r"\(([A-Za-z][A-Za-z0-9_-]{0,20})\)", text)
    if m:
        return m.group(1)

    # 5. First word-like token (skip '>' prefix)
    m = re.match(r"[>\s]*([A-Za-z][A-Za-z0-9_-]{1,20})", text)
    if m:
        return m.group(1)

    return "Unknown"


def clean_motif_sequence(text: str) -> str:
    """Remove probability annotations and non-amino-acid characters."""
    text = re.sub(r"\(0\.\d+\)", "", text)
    text = re.sub(r"\(1\)", "", text)
    text = re.sub(r"[^A-Z]", "", text.upper())
    return text.strip()


def extract_kinases_from_info(info_text: str) -> Tuple[str, str]:
    """Extract kinase names and confidence scores from an annotation string.

    Handles both legacy format (``PKA, CaMK2``) and v2 scored format
    (``PKA:0.70, CaMK2:0.65``).

    Returns
    -------
    tuple[str, str]
        ``(comma_separated_kinases, comma_separated_confidences)``
        Confidence is empty-string per kinase when scores are absent.
    """
    found_kinases: List[str] = []
    found_scores: List[str] = []
    upper_info = info_text.upper()

    for key in KINASE_KEYS:
        # Try v2 scored format first: "KINASE:0.XX"
        score_pat = re.compile(re.escape(key) + r":(\d+\.\d+)", re.IGNORECASE)
        score_match = score_pat.search(info_text)
        if score_match:
            if key not in found_kinases:
                found_kinases.append(key)
                found_scores.append(score_match.group(1))
        elif key.upper() in upper_info:
            # Legacy format (no score)
            if key not in found_kinases:
                found_kinases.append(key)
                found_scores.append("")

    return ", ".join(found_kinases), ", ".join(found_scores)


# ---------------------------------------------------------------------------
# Main extraction function
# ---------------------------------------------------------------------------
def extract_data_from_docx(word_file: str) -> pd.DataFrame:
    """Extract structured phosphorylation data from a mapped Word document.

    Parameters
    ----------
    word_file : str
        Path to the mapped/annotated Word (.docx) file.

    Returns
    -------
    pd.DataFrame
        Columns: ``Gene Name``, ``Motif``, ``Kinase Name``, ``Confidence``.
    """
    doc = Document(word_file)
    extracted_data: list = []
    current_gene_symbol: str = "Unknown"

    for para in doc.paragraphs:
        text = para.text.strip()
        if not text:
            continue

        # Gene / protein header line
        if is_header_line(text):
            current_gene_symbol = extract_gene_symbol(text)
            logger.debug("Gene header detected: %s", current_gene_symbol)
            continue

        # Peptide line with annotations
        if "(" in text and ")" in text:
            parts = text.rsplit("(", 1)
            motif_part = parts[0].strip()
            info_part = parts[1].strip() if len(parts) > 1 else ""

            clean_seq = clean_motif_sequence(motif_part)
            kinases, confidences = extract_kinases_from_info(info_part)

            if clean_seq and 3 < len(clean_seq) < 100:
                extracted_data.append(
                    {
                        "Gene Name": current_gene_symbol,
                        "Motif": clean_seq,
                        "Kinase Name": kinases,
                        "Confidence": confidences,
                    }
                )

    df = pd.DataFrame(extracted_data)
    logger.info("Extracted %d motif entries from %s", len(df), word_file)
    return df
