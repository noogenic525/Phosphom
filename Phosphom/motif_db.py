"""
Kinase substrate recognition motif database with confidence scoring.

Each motif entry encapsulates:
  - kinase      : Kinase name
  - pattern     : Pre-compiled regex for the substrate consensus motif
  - specificity : Confidence score (0.0–1.0) reflecting pattern specificity
  - description : Human-readable motif notation (Φ = hydrophobic)

Scoring rationale
-----------------
  0.80–1.00  Highly specific – multiple constrained positions, unique motif
  0.60–0.79  Moderately specific – several constraints
  0.40–0.59  Less specific – few constraints, relatively common
  0.20–0.39  Very loose – needs additional biological evidence (e.g. priming)

Key improvements (v2.1)
-----------------------
  - All regex patterns pre-compiled at module load for performance
  - PKA : Broadened to dual-alternative pattern for better coverage
  - AMPK : Added secondary alternative pattern
  - RSK : Pattern updated per consensus reference
  - CK1 : Merged acidophilic + primed into single alternation pattern (0.45)
  - CaMK2 : Hydrophobic at +1 now optional (specificity lowered to 0.55)
  - GSK-3β, CK1 primed variants remain low confidence (priming uncheckable)
  - Window: ±7 residues (15-mer), matching PSP SITE_+/-7_AA format

References
----------
  Huttlin et al., Cell 2010 – Basophilic kinase motifs
  Dephoure et al., PNAS 2008 – Proline-directed motifs
  PhosphoSitePlus (www.phosphosite.org)
  O'Shea et al., Nat Methods 2013
"""

from __future__ import annotations

import re
from typing import Dict, List, NamedTuple, Tuple


# ---------------------------------------------------------------------------
# Data type
# ---------------------------------------------------------------------------
class MotifEntry(NamedTuple):
    """Single kinase recognition motif entry."""

    kinase: str
    pattern: re.Pattern  # pre-compiled regex
    specificity: float  # 0.0 – 1.0
    description: str


# ---------------------------------------------------------------------------
# Fragment window size
# ---------------------------------------------------------------------------
WINDOW_HALF: int = 7  # ±7 residues → 15-mer window


# ---------------------------------------------------------------------------
# Motif database  (31 kinases, 31 patterns)
# ---------------------------------------------------------------------------
MOTIF_DB: List[MotifEntry] = [
    # ── Basophilic & IIS Signaling ────────────────────────────────────────
    MotifEntry(
        "PKA",
        re.compile(r"([RK][RK].([ST]))|([RK].{1,2}([ST]))"),
        0.60,
        "[R/K][R/K]-X-S/T | [R/K]-X(1-2)-S/T",
    ),
    MotifEntry(
        "AKT",
        re.compile(r"R.[RK]..([ST])(?!P)[LIVMFY]"),
        0.85,
        "R-X-[R/K]-X-X-S/T-!P-Φ",
    ),
    MotifEntry(
        "AMPK",
        re.compile(r"([LIVMF].R.([ST])(?!P).{2}[LIVMF])|([LIVMF].R..([ST])(?!P))"),
        0.70,
        "Φ-X-R-X-S/T-!P-X-X-Φ | Φ-X-R-X-X-S/T-!P",
    ),
    MotifEntry(
        "S6K1",
        re.compile(r"[RK]R..([ST])[LIVMF]"),
        0.75,
        "[R/K]-R-X-X-S/T-Φ",
    ),
    MotifEntry(
        "SGK1",
        re.compile(r"R.R..([ST])[LIV]"),
        0.80,
        "R-X-R-X-X-S/T-[L/I/V]",
    ),
    MotifEntry(
        "RSK",
        re.compile(r"[RK]R.([ST])"),
        0.60,
        "[R/K]-R-X-S/T",
    ),
    MotifEntry(
        "PKC",
        re.compile(r"[RK].([ST])[LIV][RK]"),
        0.75,
        "[R/K]-X-S/T-Φ-[R/K]",
    ),
    MotifEntry(
        "PKD",
        re.compile(r"L.R..([ST])"),
        0.70,
        "L-X-R-X-X-S/T",
    ),
    MotifEntry(
        "LKB1",
        re.compile(r"[LIVM][RK].([ST]).{2}[LIVM]"),
        0.70,
        "Φ-[R/K]-X-S/T-X-X-Φ",
    ),
    # ── Proline-directed (MAPK & Cell Cycle) ──────────────────────────────
    MotifEntry(
        "CDK",
        re.compile(r"([ST])P.[RK]"),
        0.80,
        "S/T-P-X-[K/R]",
    ),
    MotifEntry(
        "Erk",
        re.compile(r"[PLV].([ST])P"),
        0.70,
        "[P/L/V]-X-S/T-P",
    ),
    MotifEntry(
        "JNK",
        re.compile(r"P.([ST])P"),
        0.75,
        "P-X-S/T-P",
    ),
    MotifEntry(
        "p38&MAPK",
        re.compile(r"[LIVMFY].([ST])P"),
        0.60,
        "Φ-X-S/T-P",
    ),
    MotifEntry(
        "mTOR",
        re.compile(r"([ST])P[FLIV]"),
        0.65,
        "S/T-P-Φ",
    ),
    MotifEntry(
        "GSK-3beta",
        re.compile(r"([ST]).{3}[ST]"),
        0.30,
        "S/T-X-X-X-S/T (priming required — high FP without evidence)",
    ),
    MotifEntry(
        "DYRK1A",
        re.compile(r"R..([ST])P"),
        0.75,
        "R-X-X-S/T-P",
    ),
    MotifEntry(
        "HIPK2",
        re.compile(r"([ST])P.K"),
        0.80,
        "S/T-P-X-K",
    ),
    MotifEntry(
        "BUB1",
        re.compile(r"([ST])P.R"),
        0.70,
        "S/T-P-X-R",
    ),
    # ── Acidophilic & DNA Damage ──────────────────────────────────────────
    MotifEntry(
        "CK2",
        re.compile(r"([ST])(?!P).{1,2}[DE]{2,3}"),
        0.75,
        "S/T-!P-X(1-2)-[D/E](2-3)",
    ),
    MotifEntry(
        "CK1",
        re.compile(r"([DE]{1,2}.{1,2}([ST]))|(([ST]).{2,3}([ST]))"),
        0.45,
        "[D/E](1-2)-X(1-2)-S/T | S/T-X(2-3)-S/T (acidophilic + primed)",
    ),
    MotifEntry(
        "PLK1",
        re.compile(r"[DE].([ST])[LIVMF]"),
        0.75,
        "[D/E]-X-S/T-Φ",
    ),
    MotifEntry(
        "ATM&ATR",
        re.compile(r"([ST])Q"),
        0.55,
        "S/T-Q (broad consensus; DNA-PK is more specific)",
    ),
    MotifEntry(
        "DNA-PK",
        re.compile(r"([ST])Q[DE]"),
        0.80,
        "S/T-Q-[D/E]",
    ),
    MotifEntry(
        "IKK",
        re.compile(r"[DE].{1,2}([ST])G.([ST])"),
        0.80,
        "[D/E]-X(1-2)-S-G-X-S",
    ),
    # ── Other Key Kinases ─────────────────────────────────────────────────
    MotifEntry(
        "CaMK2",
        re.compile(r"[RK]..([ST])([LIVMF])?"),
        0.55,
        "[R/K]-X-X-S/T-Φ? (Φ optional)",
    ),
    MotifEntry(
        "Aurora A",
        re.compile(r"[RK].([ST])[LIVMF]"),
        0.65,
        "[R/K]-X-S/T-Φ",
    ),
    MotifEntry(
        "Aurora B",
        re.compile(r"[RK]R.([ST])[LIVMF]"),
        0.75,
        "[R/K]-R-X-S/T-Φ",
    ),
    MotifEntry(
        "Chk1",
        re.compile(r"[LIVMF]R..([ST])"),
        0.65,
        "Φ-R-X-X-S/T",
    ),
    MotifEntry(
        "Chk2",
        re.compile(r"R..([ST])[LIVMF]"),
        0.65,
        "R-X-X-S/T-Φ",
    ),
    MotifEntry(
        "NEK2",
        re.compile(r"[LIVMF]R..([ST])"),
        0.65,
        "Φ-R-X-X-S/T",
    ),
    MotifEntry(
        "MK2",
        re.compile(r"[LIVM].R.([ST])"),
        0.65,
        "[L/I/V/M]-X-R-X-S/T",
    ),
]


# ---------------------------------------------------------------------------
# Auto-derived kinase keys (preserving first-seen order, no duplicates)
# ---------------------------------------------------------------------------
KINASE_KEYS: List[str] = list(dict.fromkeys(e.kinase for e in MOTIF_DB))


# ---------------------------------------------------------------------------
# Cell signaling pathway classification
# ---------------------------------------------------------------------------
CELL_SIGNAL_CATEGORY_ORDER: List[str] = [
    "Growth & Metabolism",
    "Stress Response",
    "DNA Damage Response (DDR)",
    "Cell Cycle & Mitosis",
    "Energy Homeostasis",
    "Second Messenger Signaling",
    "Inflammation & Immunity",
    "Multi-functional",
]

TOP_KINASES_N: int = 5


def _normalize_kinase_key(name: str) -> str:
    """Normalize kinase name for dictionary look-ups (lowercase, no spaces)."""
    return re.sub(r"\s+", "", str(name)).lower()


CELL_SIGNAL_MAP: Dict[str, str] = {
    # Growth & Metabolism
    _normalize_kinase_key("AKT"): "Growth & Metabolism",
    _normalize_kinase_key("mTOR"): "Growth & Metabolism",
    _normalize_kinase_key("S6K1"): "Growth & Metabolism",
    _normalize_kinase_key("SGK1"): "Growth & Metabolism",
    _normalize_kinase_key("GSK-3beta"): "Growth & Metabolism",
    # Stress Response
    _normalize_kinase_key("Erk"): "Stress Response",
    _normalize_kinase_key("RSK"): "Stress Response",
    _normalize_kinase_key("p38&MAPK"): "Stress Response",
    _normalize_kinase_key("JNK"): "Stress Response",
    _normalize_kinase_key("MK2"): "Stress Response",
    # DNA Damage Response
    _normalize_kinase_key("ATM&ATR"): "DNA Damage Response (DDR)",
    _normalize_kinase_key("DNA-PK"): "DNA Damage Response (DDR)",
    _normalize_kinase_key("Chk1"): "DNA Damage Response (DDR)",
    _normalize_kinase_key("Chk2"): "DNA Damage Response (DDR)",
    _normalize_kinase_key("HIPK2"): "DNA Damage Response (DDR)",
    # Cell Cycle & Mitosis
    _normalize_kinase_key("CDK"): "Cell Cycle & Mitosis",
    _normalize_kinase_key("PLK1"): "Cell Cycle & Mitosis",
    _normalize_kinase_key("Aurora A"): "Cell Cycle & Mitosis",
    _normalize_kinase_key("Aurora B"): "Cell Cycle & Mitosis",
    _normalize_kinase_key("NEK2"): "Cell Cycle & Mitosis",
    _normalize_kinase_key("BUB1"): "Cell Cycle & Mitosis",
    # Energy Homeostasis
    _normalize_kinase_key("AMPK"): "Energy Homeostasis",
    _normalize_kinase_key("LKB1"): "Energy Homeostasis",
    # Second Messenger Signaling
    _normalize_kinase_key("PKA"): "Second Messenger Signaling",
    _normalize_kinase_key("PKC"): "Second Messenger Signaling",
    _normalize_kinase_key("PKD"): "Second Messenger Signaling",
    _normalize_kinase_key("CaMK2"): "Second Messenger Signaling",
    # Inflammation & Immunity
    _normalize_kinase_key("IKK"): "Inflammation & Immunity",
    # Multi-functional
    _normalize_kinase_key("CK1"): "Multi-functional",
    _normalize_kinase_key("CK2"): "Multi-functional",
    _normalize_kinase_key("DYRK1A"): "Multi-functional",
}


# ---------------------------------------------------------------------------
# Core matching function
# ---------------------------------------------------------------------------
def identify_kinases(
    full_seq: str,
    index: int,
    min_confidence: float = 0.0,
) -> List[Tuple[str, float]]:
    """Identify potential kinases for a single phosphorylation site.

    Extracts a ±7 residue window (15-mer) centred on *index* in the full
    protein sequence, then tests every motif pattern in ``MOTIF_DB``.

    Parameters
    ----------
    full_seq : str
        Full protein amino-acid sequence (uppercase, no gaps).
    index : int
        0-based position of the phospho-residue (S/T/Y) in *full_seq*.
    min_confidence : float, optional
        Minimum specificity score to include in results (default 0.0).

    Returns
    -------
    list[tuple[str, float]]
        ``(kinase_name, confidence)`` pairs sorted by confidence descending.
        If the same kinase matches via multiple patterns, only the highest
        score is kept.
    """
    start = max(0, index - WINDOW_HALF)
    end = min(len(full_seq), index + WINDOW_HALF + 1)
    fragment = full_seq[start:end]

    hits: Dict[str, float] = {}
    for entry in MOTIF_DB:
        if entry.pattern.search(fragment):
            if entry.specificity >= min_confidence:
                # Keep highest score per kinase
                if entry.kinase not in hits or entry.specificity > hits[entry.kinase]:
                    hits[entry.kinase] = entry.specificity

    return sorted(hits.items(), key=lambda x: -x[1])
