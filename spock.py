#!/usr/bin/env python3
"""
spock.py  (v5 — improved soft mode)

Given a set of per-marker FASTQ.GZ files (produced by extract_ssr_amplicons_3.py),
this script:
  1. Builds a read-length histogram for each of the nine SSR markers.
  2. Applies a minimum-depth filter (default: >= 1000 reads per marker).
  3. Detects peaks, suppresses polymerase stutter (+-2/4/6 bp satellites below
     a stutter-ratio threshold), and calls the one or two true alleles.
  4. Optionally applies per-marker calibration offsets (CE = raw + offset).
  5. Plots size distributions with called alleles annotated on the same scale
     selected for the matrix output (--matrix-scale: raw or ce).
  6. Writes a genotype profile to stdout and to <o>.genotype.tsv.
  7. Writes a count matrix to <o>.count_matrix.tsv:
       rows   = every integer size between --min-length and --max-length
       columns = the nine markers in VIVC order
       cells  = read count at that size (raw or CE-calibrated per --matrix-scale)

Modes for --homozygosity-mode:
  'ambiguous'  — legacy behaviour: close allele pairs (<=4 bp) are reported
                 as the dominant allele only, flagged "ambiguous".
  'soft'       — improved calling: a stutter-aware two-peak model is fit to
                 the on-parity read-length distribution; the second allele is
                 only called when the two-peak model beats a single-peak
                 (homozygous) model by a meaningful margin AND the candidate
                 survives offset-specific anti-stutter thresholds. Designed
                 to recover close alleles (diff 2-6 bp) that the original
                 heuristic cannot.

Plot scaling:
  --matrix-scale raw -> plot x-axis, peak positions, and peak labels are raw
  --matrix-scale ce  -> plot x-axis, peak positions, and peak labels are CE-calibrated

Marker order (VIVC standard):
  VVS2, VVMD5, VVMD7, VVMD25, VVMD27, VVMD28, VVMD32, VrZAG62, VrZAG79

Calibration offsets (CE = raw + offset):
  VVS2: 0, VVMD5: +3, VVMD7: +1, VVMD25: -2, VVMD27: +2,
  VVMD28: -1, VVMD32: +1, VrZAG62: 0, VrZAG79: +1

Column assignment in output:
  _1 = smaller bp size, _2 = larger bp size
  Brackets () mark the weaker peak (lower read count).
  The dominant peak (highest read count) is never bracketed.

Usage:
  python spock.py \\
      --basename 01 \\
      --output   01 \\
      [--min-reads 1000] \\
      [--stutter-ratio 0.10] \\
      [--homozygosity-mode ambiguous|soft] \\
      [--homo-threshold 0.30] \\
      [--min-length 100] \\
      [--max-length 320] \\
      [--matrix-scale raw|ce] \\
      [--soft-model-improvement 0.25] \\
      [--soft-min-h2-ratio 0.08] \\
      [--soft-cross-parity-markers VVMD27,VVMD32]
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import Counter
from typing import Dict, Iterable, List, Optional, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.signal import find_peaks

MARKERS: List[str] = [
    "VVS2", "VVMD5", "VVMD7", "VVMD25",
    "VVMD27", "VVMD28", "VVMD32", "VrZAG62", "VrZAG79",
]

# CE = raw + offset
OFFSETS: Dict[str, int] = {
    "VVS2":    0,
    "VVMD5":   3,
    "VVMD7":   1,
    "VVMD25": -2,
    "VVMD27":  2,
    "VVMD28": -1,
    "VVMD32":  1,
    "VrZAG62": 0,
    "VrZAG79": 1,
}

MARKER_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
]

# ── Empirical stutter profile ───────────────────────────────────────
# Derived from 23 known homozygotes across 9 markers in the training set.
# Values are ratios relative to the main peak height.
# TYPICAL = median observed; MAX = near-upper-bound observed.
# Scaling with dominant size accounts for the fact that very long repeats
# (e.g. VVMD32 at 272 bp) show much heavier stutter (−2 ratio reaching 0.99).

STUTTER_TYPICAL_BASE: Dict[int, float] = {
    -10: 0.03, -8: 0.05, -6: 0.12, -4: 0.30, -2: 0.60, 0: 1.00,
      2: 0.30,   4: 0.11,   6: 0.05,   8: 0.02,  10: 0.01,
}
STUTTER_MAX_BASE: Dict[int, float] = {
    -10: 0.10, -8: 0.15, -6: 0.40, -4: 0.60, -2: 0.85, 0: 1.00,
      2: 0.50,   4: 0.25,   6: 0.15,   8: 0.05,  10: 0.03,
}

# Markers whose reference panel contains alleles at both parities.
# For these, soft mode also considers opposite-parity peaks as possible
# second-allele candidates (with strict additional criteria).
DEFAULT_CROSS_PARITY_MARKERS: Set[str] = {"VVMD27", "VVMD32"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Call SSR alleles from extracted amplicon FASTQ files."
    )
    p.add_argument("--basename", required=True,
                   help="Input basename, e.g. '01' -> reads 01.VVS2.fastq.gz etc.")
    p.add_argument("--output", required=True,
                   help="Output basename for PNG, TSV, and matrix files.")
    p.add_argument("--min-reads", type=int, default=1000,
                   help="Minimum total reads per marker to attempt allele calling (default: 1000).")
    p.add_argument("--stutter-ratio", type=float, default=0.10,
                   help="Peak height fraction below which a +-2/4/6 bp satellite is "
                        "treated as polymerase stutter and suppressed (default: 0.10). "
                        "Used in ambiguous mode only; soft mode uses its own stutter model.")
    p.add_argument("--homozygosity-mode", choices=["ambiguous", "soft"],
                   default="ambiguous",
                   help="How to handle close allele pairs (difference <= 4 bp): "
                        "'ambiguous' reports only the dominant allele; "
                        "'soft' uses a stutter-aware two-peak model to recover "
                        "the second allele when justified by the data "
                        "(default: ambiguous).")
    p.add_argument("--homo-threshold", type=float, default=0.30,
                   help="If the second candidate peak is below this fraction of the "
                        "dominant peak it is treated as homozygous (default: 0.30). "
                        "Used in ambiguous mode only.")
    p.add_argument("--min-length", type=int, default=100,
                   help="Minimum fragment length to include (default: 100).")
    p.add_argument("--max-length", type=int, default=320,
                   help="Maximum fragment length to include (default: 320).")
    p.add_argument("--matrix-scale", choices=["raw", "ce"], default="raw",
                   help="Whether the count matrix row indices (sizes) should be on the "
                        "raw sequencing scale or the CE-calibrated scale (default: raw).")
    # Soft-mode knobs (safe defaults match the training-set optimum)
    p.add_argument("--soft-model-improvement", type=float, default=0.25,
                   help="[soft mode] Minimum relative drop in residual score that the "
                        "two-peak model must achieve over the one-peak model (default: 0.25).")
    p.add_argument("--soft-min-h2-ratio", type=float, default=0.08,
                   help="[soft mode] Minimum fitted height of second allele relative to "
                        "first allele (default: 0.08).")
    p.add_argument("--soft-min-raw-ratio-margin", type=float, default=0.05,
                   help="[soft mode] Required margin by which the candidate's raw ratio "
                        "must exceed the maximum plausible stutter ratio at that offset "
                        "(default: 0.05).")
    p.add_argument("--soft-cross-parity-markers", type=str,
                   default=",".join(sorted(DEFAULT_CROSS_PARITY_MARKERS)),
                   help="[soft mode] Comma-separated list of markers that may carry "
                        "alleles of either parity (default: VVMD27,VVMD32). "
                        "Set to an empty string to disable cross-parity candidate search.")
    p.add_argument("--soft-bracket-threshold", type=float, default=0.30,
                   help="[soft mode] If the fitted secondary/dominant height ratio is "
                        "below this value, the secondary allele is bracketed as "
                        "low-confidence (default: 0.30).")
    return p.parse_args()


# ── I/O helpers ──────────────────────────────────────────────────

def iter_fastq_lengths(path: str) -> Iterable[int]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip()
            fh.readline()
            fh.readline()
            if not seq:
                raise ValueError(f"Incomplete FASTQ record in {path}")
            yield len(seq)


def load_length_counts(path: str) -> Counter:
    c: Counter = Counter()
    for length in iter_fastq_lengths(path):
        c[length] += 1
    return c


def ensure_parent_dir(path: str) -> None:
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


# ── Peak detection ───────────────────────────────────────────────

def suppress_stutter(peaks: List[int], counts: Counter, stutter_ratio: float) -> List[int]:
    stutter_offsets = {2, 4, 6}
    changed = True
    active = set(peaks)
    while changed:
        changed = False
        to_remove = set()
        for pk in list(active):
            for offset in stutter_offsets:
                for parent in (pk - offset, pk + offset):
                    if parent in active and counts[parent] > counts[pk]:
                        if counts[pk] < stutter_ratio * counts[parent]:
                            to_remove.add(pk)
                            break
                if pk in to_remove:
                    break
        if to_remove:
            active -= to_remove
            changed = True
    return sorted(active)


def detect_peaks(
    counts: Counter,
    min_length: int,
    max_length: int,
    stutter_ratio: float,
) -> List[int]:
    xs = list(range(min_length, max_length + 1))
    ys = [counts.get(x, 0) for x in xs]
    global_max = max(ys) if ys else 1
    min_prominence = max(10, 0.05 * global_max)
    raw_indices, _ = find_peaks(ys, prominence=min_prominence, distance=2)
    raw_peaks = [xs[i] for i in raw_indices]
    return suppress_stutter(raw_peaks, counts, stutter_ratio)


# ── Allele calling ───────────────────────────────────────────────
#
# Two separate call paths:
#   * call_alleles_legacy(): original heuristic used when
#                            --homozygosity-mode ambiguous (unchanged).
#   * call_alleles_soft():   improved stutter-aware model used when
#                            --homozygosity-mode soft.
#
# Both return (dominant_raw, secondary_raw, flag) in the raw-size space.

def call_alleles_legacy(
    peaks: List[int],
    counts: Counter,
    homozygosity_mode: str,
    homo_threshold: float,
) -> Tuple[Optional[int], Optional[int], str]:
    """Original ambiguous-mode calling logic (kept intact)."""
    if not peaks:
        return None, None, "no_call"

    peaks_by_count = sorted(peaks, key=lambda p: counts[p], reverse=True)
    dominant = peaks_by_count[0]

    if len(peaks_by_count) == 1:
        return dominant, None, "homozygous"

    secondary = peaks_by_count[1]

    if counts[secondary] < homo_threshold * counts[dominant]:
        return dominant, None, "homozygous"

    diff = abs(dominant - secondary)
    if diff <= 4:
        if homozygosity_mode == "ambiguous":
            return dominant, None, "ambiguous"
        else:
            return dominant, secondary, "soft"

    return dominant, secondary, "heterozygous"


# -- Soft-mode helpers --

def _size_adjusted(base: Dict[int, float], dom_size: int,
                   scale_strength: float) -> Dict[int, float]:
    """Scale the negative-offset stutter ratios upward for large dominant sizes.
    
    Empirically, very long repeats (e.g. VVMD32 at 272 bp) show much heavier
    stutter than mid-range repeats. This applies a linear bump above 240 bp.
    """
    if dom_size <= 240:
        return dict(base)
    scale = 1.0 + scale_strength * (dom_size - 240)
    result = {}
    for off, r in base.items():
        if off < 0:
            result[off] = min(0.99, r * scale)
        else:
            result[off] = r
    return result


def _infer_parity(counts: Counter, min_length: int, max_length: int) -> int:
    """Return 0 or 1 - whichever parity carries the bulk of reads in [min_length, max_length]."""
    even = sum(c for s, c in counts.items()
               if min_length <= s <= max_length and s % 2 == 0)
    odd = sum(c for s, c in counts.items()
              if min_length <= s <= max_length and s % 2 == 1)
    return 0 if even >= odd else 1


def _on_parity(counts: Counter, parity: int,
               min_length: int, max_length: int) -> Counter:
    """Return counts restricted to positions with the given parity and within the length bounds."""
    return Counter({s: c for s, c in counts.items()
                    if min_length <= s <= max_length and s % 2 == parity})


def _fit_two_peak_heights(peak_positions: List[int], observed: Counter,
                          stutter_typical: Dict[int, float],
                          n_iters: int = 5) -> Dict[int, float]:
    """Iterative NNLS-style fit: each peak's 'true' height equals its observed
    count minus the predicted stutter contribution from the other peaks."""
    heights = {p: max(1.0, float(observed.get(p, 0))) for p in peak_positions}
    for _ in range(n_iters):
        new_heights = {}
        for p in peak_positions:
            obs = observed.get(p, 0)
            contrib = 0.0
            for q in peak_positions:
                if q == p:
                    continue
                off = p - q
                contrib += heights[q] * stutter_typical.get(off, 0.0)
            new_heights[p] = max(0.0, obs - contrib)
        heights = new_heights
    return heights


def _residual_score(peak_positions: List[int], observed: Counter,
                    window: List[int],
                    stutter_typical: Dict[int, float]) -> Tuple[float, Dict[int, float]]:
    """Fit heights, build predicted landscape, return (normalized squared residual, heights)."""
    heights = _fit_two_peak_heights(peak_positions, observed, stutter_typical)
    pred = {p: 0.0 for p in window}
    for peak, h in heights.items():
        for off, ratio in stutter_typical.items():
            target = peak + off
            if target in pred:
                pred[target] += h * ratio
    residual = 0.0
    obs_sq_sum = 0.0
    for pos in window:
        obs_val = observed.get(pos, 0)
        pred_val = pred.get(pos, 0.0)
        residual += (obs_val - pred_val) ** 2
        obs_sq_sum += obs_val ** 2
    return residual / max(1.0, obs_sq_sum), heights


def _offset_specific_ok(offset: int, raw_ratio: float,
                        max_stutter_ratio: float, h2_ratio: float,
                        improvement: float, margin: float) -> bool:
    """Offset-specific anti-stutter test. Stutter tails are strongest at -2
    and decay going outward. The raw ratio of a candidate must exceed the
    maximum plausible stutter by a margin, OR the two-peak model must be a
    strong improvement AND the fitted height ratio must be substantial."""
    ratio_cutoff = max_stutter_ratio + margin

    if offset == 0:
        return False  # same as dominant, never

    if offset == -2:
        return raw_ratio >= ratio_cutoff or (h2_ratio >= 0.30 and improvement >= 0.40)
    if offset == +2:
        return raw_ratio >= ratio_cutoff or (h2_ratio >= 0.15 and improvement >= 0.30)
    if offset == -4:
        return raw_ratio >= ratio_cutoff or (h2_ratio >= 0.25 and improvement >= 0.40)
    if offset == +4:
        return raw_ratio >= ratio_cutoff or (h2_ratio >= 0.15 and improvement >= 0.30)
    if offset == -6:
        return raw_ratio >= ratio_cutoff or (h2_ratio >= 0.25 and improvement >= 0.40)
    if offset == +6:
        return raw_ratio >= ratio_cutoff or (h2_ratio >= 0.15 and improvement >= 0.30)
    if abs(offset) in (8, 10):
        if offset < 0:
            return raw_ratio >= ratio_cutoff or (h2_ratio >= 0.20 and improvement >= 0.40)
        return raw_ratio >= ratio_cutoff or (h2_ratio >= 0.15 and improvement >= 0.30)
    # |offset| >= 12: outside stutter ladder, still require meaningful raw ratio
    return raw_ratio >= 0.08


def _is_tail_stutter(candidate: int, dominant: int, counts: Counter,
                     dom_count: int, stutter_typical: Dict[int, float]) -> bool:
    """Reject candidates that sit on a monotonically-increasing stutter tail
    and whose count is near what typical stutter would predict. These are
    very likely part of an extended stutter ladder, not a second allele."""
    offset = candidate - dominant
    if not (-10 <= offset <= -2):
        return False
    # Check that positions between candidate and dominant form a monotone ramp
    prev = counts.get(candidate, 0)
    step = 2
    for p in range(candidate + step, dominant, step):
        cur = counts.get(p, 0)
        if cur <= prev:
            return False
        prev = cur
    # Check whether candidate is close to predicted typical stutter
    expected = stutter_typical.get(offset, 0.0) * dom_count
    return counts.get(candidate, 0) < 1.3 * expected


def call_alleles_soft(
    counts: Counter,
    marker: str,
    min_length: int,
    max_length: int,
    model_improvement: float,
    min_h2_ratio: float,
    min_raw_ratio_margin: float,
    cross_parity_markers: Set[str],
) -> Tuple[Optional[int], Optional[int], str, float]:
    """Improved soft-mode allele caller using a stutter-aware two-peak model.

    Returns (dominant_raw, secondary_raw, flag, h2_ratio).
    flag is one of: 'heterozygous', 'soft', 'homozygous', 'low_evidence', 'no_call'.

    'heterozygous' is returned when the two alleles differ by >6 bp (far-apart pair,
    the easy case). 'soft' is returned when the two alleles are within 6 bp (close
    pair that the legacy algorithm would have marked 'ambiguous'). 'homozygous' is
    returned when the model fits better as a single peak.

    The fourth return value (h2_ratio) is the fitted secondary/dominant height ratio
    and can be used downstream for bracketing low-confidence calls.
    """
    total = sum(counts.values())
    if total == 0:
        return None, None, "no_call", 0.0

    parity = _infer_parity(counts, min_length, max_length)
    op = _on_parity(counts, parity, min_length, max_length)
    if not op or max(op.values()) == 0:
        return None, None, "no_call", 0.0

    dominant = max(op, key=lambda x: op[x])
    dom_count = op[dominant]

    if dom_count < 20:
        return dominant, None, "low_evidence", 0.0

    # Size-adjusted stutter profile
    stutter_typical = _size_adjusted(STUTTER_TYPICAL_BASE, dominant, 0.010)
    stutter_max     = _size_adjusted(STUTTER_MAX_BASE,     dominant, 0.015)

    # Analysis window: span all positions with non-negligible signal, extended
    # by the stutter reach (±10).
    significant = [p for p, c in op.items() if c >= 0.02 * dom_count]
    if not significant:
        return dominant, dominant, "homozygous", 0.0
    w_min = max(min_length, min(significant) - 10)
    w_max = min(max_length, max(significant) + 10)
    window = [p for p in range(w_min, w_max + 1) if p % 2 == parity]

    # Candidate second-allele positions (same parity as dominant)
    min_abs = max(15, 0.05 * dom_count)
    candidates = [p for p in window if p != dominant and op.get(p, 0) >= min_abs]

    # Cross-parity candidates (only for designated markers)
    cross_parity_counts: Counter = Counter()
    if marker in cross_parity_markers:
        cp_parity = 1 - parity
        cp_total = sum(c for s, c in counts.items()
                       if min_length <= s <= max_length and s % 2 == cp_parity)
        # Only consider cross-parity if the opposite parity carries a decent
        # fraction of reads (cross-parity true alleles give >=25% of the main
        # parity's signal; pure off-by-one noise typically gives <25%).
        if cp_total >= 0.25 * sum(op.values()):
            cross_parity_counts = Counter({s: c for s, c in counts.items()
                                           if min_length <= s <= max_length
                                           and s % 2 == cp_parity})
            # Cross-parity candidates must satisfy ALL of:
            #   * absolute count >= 25% of dominant (clearly above off-by-one noise,
            #     which reaches up to ~20% typically, up to ~44% in outliers)
            #   * sit at least 3 bp away from dominant (|offset| >= 3) so we
            #     are not just picking up direct +/-1 noise of the dominant
            #   * be a local maximum on the cross-parity grid
            for p, c in cross_parity_counts.items():
                if c < 0.25 * dom_count:
                    continue
                if abs(p - dominant) < 3:
                    continue
                left = cross_parity_counts.get(p - 2, 0)
                right = cross_parity_counts.get(p + 2, 0)
                if c >= left and c >= right:
                    candidates.append(p)

    # Merge counts (primary + cross-parity) for model fitting
    merged = Counter(op)
    for p, c in cross_parity_counts.items():
        merged[p] = c
    # Expand window to include any cross-parity candidate positions
    if cross_parity_counts:
        ext_window = set(window)
        for p in cross_parity_counts:
            ext_window.add(p)
        window = sorted(ext_window)

    # One-peak (homozygous) model score
    score_hom, _ = _residual_score([dominant], merged, window, stutter_typical)

    # Best two-peak model
    best_score = float("inf")
    best_pos = None
    best_heights: Dict[int, float] = {}
    for cand in candidates:
        score, heights = _residual_score([dominant, cand], merged, window, stutter_typical)
        if score < best_score:
            best_score = score
            best_pos = cand
            best_heights = heights

    if best_pos is None:
        return dominant, dominant, "homozygous", 0.0

    h1 = best_heights.get(dominant, 0.0)
    h2 = best_heights.get(best_pos, 0.0)
    improvement = (score_hom - best_score) / max(1e-9, score_hom)
    h2_over_h1 = h2 / max(1.0, h1)
    raw_ratio = merged[best_pos] / max(1.0, dom_count)
    offset = best_pos - dominant
    cross_parity = (best_pos % 2) != parity

    # Core criteria: the 2-peak model beats the 1-peak model, fitted h2 is substantial
    base_ok = (
        improvement >= model_improvement
        and h2_over_h1 >= min_h2_ratio
        and h2 >= max(30, 0.02 * total)
    )

    # Offset-specific anti-stutter test (only applies to same-parity candidates)
    if cross_parity:
        # Cross-parity candidates bypass the stutter-offset test but require
        # stronger absolute evidence (already enforced in candidate selection).
        off_ok = raw_ratio >= 0.25
    else:
        max_stut = stutter_max.get(offset, 0.0)
        off_ok = _offset_specific_ok(offset, raw_ratio, max_stut,
                                     h2_over_h1, improvement, min_raw_ratio_margin)

    het_ok = base_ok and off_ok

    # Tail guard: candidate lies on a smooth stutter tail and matches typical stutter -> reject
    if het_ok and not cross_parity:
        if _is_tail_stutter(best_pos, dominant, op, dom_count, stutter_typical):
            het_ok = False

    if not het_ok:
        return dominant, dominant, "homozygous", 0.0

    # Distinguish far-apart hets from close hets for the downstream flag;
    # close hets inherit the 'soft' label from the legacy code path.
    if abs(offset) > 6 or cross_parity:
        flag = "heterozygous"
    else:
        flag = "soft"

    return dominant, best_pos, flag, h2_over_h1


# ── Calibration ──────────────────────────────────────────────────

def to_ce_size(raw: Optional[int], marker: str) -> Optional[int]:
    if raw is None:
        return None
    return raw + OFFSETS[marker]


def calibrate_counts(counts: Counter, marker: str) -> Counter:
    """Shift all size keys by the marker's CE offset."""
    offset = OFFSETS[marker]
    if offset == 0:
        return counts
    return Counter({size + offset: count for size, count in counts.items()})


# ── Formatting ───────────────────────────────────────────────────

def format_allele(
    dominant_ce: Optional[int],
    secondary_ce: Optional[int],
    flag: str,
    h2_ratio: float = 0.0,
    bracket_threshold: float = 0.30,
) -> Tuple[str, str]:
    """
    col_1 = smaller CE size, col_2 = larger CE size.
    Brackets mark the weaker (secondary) allele; dominant is never bracketed.

    In soft mode, the secondary allele is bracketed when its fitted height
    ratio is below `bracket_threshold`, signalling a low-confidence call.
    In 'soft' (close-pair) flag, the secondary is always bracketed as before.
    In 'heterozygous' flag, brackets are applied only if h2_ratio is low.
    """
    if flag in ("no_call", "low_depth", "low_evidence") or dominant_ce is None:
        return "", ""

    if secondary_ce is None:
        return str(dominant_ce), ""

    small_ce, large_ce = sorted([dominant_ce, secondary_ce])
    dominant_is_small = dominant_ce <= secondary_ce

    should_bracket = False
    if flag == "soft":
        # Close pairs always bracket the weaker allele (preserved from legacy)
        should_bracket = True
    elif flag == "heterozygous" and h2_ratio > 0 and h2_ratio < bracket_threshold:
        # Low-confidence far-apart het: bracket the weaker allele
        should_bracket = True

    if should_bracket:
        col_1 = str(small_ce) if dominant_is_small else f"({small_ce})"
        col_2 = f"({large_ce})" if dominant_is_small else str(large_ce)
    else:
        col_1 = str(small_ce)
        col_2 = str(large_ce)

    return col_1, col_2


# ── Count matrix output ──────────────────────────────────────────

def write_count_matrix(
    results: List[dict],
    output_path: str,
    min_length: int,
    max_length: int,
    matrix_scale: str,
) -> None:
    """
    Rows = every integer size from min_length to max_length.
    If matrix_scale == 'ce', counts are shifted to CE-equivalent sizes and the
    size column represents CE bp. The row range is adjusted by the min/max
    possible offset so no data is lost.
    Columns = size, then one column per marker in VIVC order.
    """
    ensure_parent_dir(output_path)

    if matrix_scale == "ce":
        counts_by_marker: Dict[str, Counter] = {
            res["marker"]: calibrate_counts(res["counts"], res["marker"])
            for res in results
        }
        all_offsets = list(OFFSETS.values())
        row_min = min_length + min(all_offsets)
        row_max = max_length + max(all_offsets)
        scale_label = "size_ce"
    else:
        counts_by_marker = {res["marker"]: res["counts"] for res in results}
        row_min = min_length
        row_max = max_length
        scale_label = "size_raw"

    with open(output_path, "wt") as fh:
        fh.write("\t".join([scale_label] + MARKERS) + "\n")
        for size in range(row_min, row_max + 1):
            row = [str(size)]
            for marker in MARKERS:
                row.append(str(counts_by_marker.get(marker, Counter()).get(size, 0)))
            fh.write("\t".join(row) + "\n")

    print(f"INFO  Count matrix ({matrix_scale} scale) written to {output_path}", file=sys.stderr)


# ── Plotting ─────────────────────────────────────────────────────

def make_plot(
    results: List[dict],
    output_png: str,
    basename: str,
    min_length: int,
    max_length: int,
    matrix_scale: str,
) -> None:
    fig, ax = plt.subplots(figsize=(16, 7))
    handles = []

    if matrix_scale == "ce":
        plot_min = min_length + min(OFFSETS.values())
        plot_max = max_length + max(OFFSETS.values())
        axis_label = "PCR product length (bp, CE-calibrated scale)"
        title_suffix = "(peak positions and labels on CE-calibrated scale)"
    else:
        plot_min = min_length
        plot_max = max_length
        axis_label = "PCR product length (bp, raw sequencing scale)"
        title_suffix = "(peak positions and labels on raw sequencing scale)"

    for i, res in enumerate(results):
        marker = res["marker"]
        color = MARKER_COLORS[i % len(MARKER_COLORS)]

        if matrix_scale == "ce":
            counts_plot = calibrate_counts(res["counts"], marker)
            dominant_plot = res["dominant_ce"]
            secondary_plot = res["secondary_ce"]
        else:
            counts_plot = res["counts"]
            dominant_plot = res["dominant_raw"]
            secondary_plot = res["secondary_raw"]

        xs_plot = sorted(
            x for x in counts_plot if plot_min <= x <= plot_max and counts_plot[x] >= 5
        )
        if not xs_plot:
            continue
        ys_plot = [counts_plot[x] for x in xs_plot]
        ax.plot(xs_plot, ys_plot, color=color, linewidth=1.5,
                marker="o", markersize=3, alpha=0.85)

        for plot_allele, is_dominant in [
            (dominant_plot, True),
            (secondary_plot, False),
        ]:
            if plot_allele is not None and plot_min <= plot_allele <= plot_max:
                h = counts_plot.get(plot_allele, 0)
                if not is_dominant and res["flag"] == "soft":
                    label_txt = f"({plot_allele})"
                elif not is_dominant and res["flag"] == "ambiguous":
                    label_txt = f"~{plot_allele}"
                else:
                    label_txt = str(plot_allele)
                ax.annotate(
                    label_txt,
                    xy=(plot_allele, h),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center", fontsize=7, color=color, fontweight="bold",
                )
                ax.axvline(plot_allele, color=color, linestyle="--",
                           linewidth=0.8, alpha=0.5)

        patch = mpatches.Patch(color=color, label=f"{marker} ({res['flag']})")
        handles.append(patch)

    major_tick_start = (min_length // 10) * 10
    if major_tick_start == min_length:
        major_tick_start -= 10
    minor_tick_start = (min_length // 2) * 2
    if minor_tick_start > plot_min:
        minor_tick_start -= 2

    ax.set_xlim(plot_min, plot_max)
    ax.set_xlabel(axis_label)
    ax.set_ylabel("Read count")
    ax.set_title(f"SSR size distributions – {os.path.basename(basename)}\n{title_suffix}")
    ax.set_xticks(list(range(major_tick_start, plot_max + 1, 10)))
    ax.set_xticks(list(range(minor_tick_start, plot_max + 1, 2)), minor=True)
    ax.tick_params(axis="x", which="major", length=6)
    ax.tick_params(axis="x", which="minor", length=3)
    for x in range(minor_tick_start, plot_max + 1, 2):
        ax.axvline(x, linestyle=":", linewidth=0.3, alpha=0.25, color="grey", zorder=0)
    ax.legend(handles=handles, loc="upper left",
              bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0, frameon=True)
    plt.tight_layout()
    ensure_parent_dir(output_png)
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"INFO  Plot written to {output_png}", file=sys.stderr)


# ── Main ─────────────────────────────────────────────────────────

def main() -> int:
    args = parse_args()

    # Parse cross-parity marker set
    cross_parity_markers = {m.strip() for m in args.soft_cross_parity_markers.split(",")
                            if m.strip()}
    unknown = cross_parity_markers - set(MARKERS)
    if unknown:
        print(f"WARN  Unknown markers in --soft-cross-parity-markers: {unknown}", file=sys.stderr)
        cross_parity_markers &= set(MARKERS)

    results: List[dict] = []

    for marker in MARKERS:
        path = f"{args.basename}.{marker}.fastq.gz"
        if not os.path.exists(path):
            print(f"WARN  {path} not found – skipping {marker}.", file=sys.stderr)
            results.append({
                "marker": marker, "counts": Counter(), "total_reads": 0,
                "dominant_raw": None, "secondary_raw": None,
                "dominant_ce": None,  "secondary_ce": None,
                "flag": "no_call", "h2_ratio": 0.0,
            })
            continue

        counts = load_length_counts(path)
        total = sum(counts.values())

        if total < args.min_reads:
            print(f"WARN  {marker}: only {total} reads (< {args.min_reads} minimum) – no call.",
                  file=sys.stderr)
            results.append({
                "marker": marker, "counts": counts, "total_reads": total,
                "dominant_raw": None, "secondary_raw": None,
                "dominant_ce": None,  "secondary_ce": None,
                "flag": "low_depth", "h2_ratio": 0.0,
            })
            continue

        peaks = detect_peaks(counts, args.min_length, args.max_length, args.stutter_ratio)
        h2_ratio = 0.0

        if args.homozygosity_mode == "soft":
            dominant_raw, secondary_raw, flag, h2_ratio = call_alleles_soft(
                counts, marker,
                args.min_length, args.max_length,
                args.soft_model_improvement,
                args.soft_min_h2_ratio,
                args.soft_min_raw_ratio_margin,
                cross_parity_markers,
            )
        else:
            dominant_raw, secondary_raw, flag = call_alleles_legacy(
                peaks, counts, args.homozygosity_mode, args.homo_threshold
            )

        dominant_ce  = to_ce_size(dominant_raw,  marker)
        secondary_ce = to_ce_size(secondary_raw, marker)

        results.append({
            "marker": marker, "counts": counts, "total_reads": total,
            "dominant_raw": dominant_raw, "secondary_raw": secondary_raw,
            "dominant_ce":  dominant_ce,  "secondary_ce":  secondary_ce,
            "flag": flag, "h2_ratio": h2_ratio,
        })

        col_1, col_2 = format_allele(dominant_ce, secondary_ce, flag,
                                     h2_ratio, args.soft_bracket_threshold)
        print(f"INFO  {marker:<10} {flag:<14} CE alleles: {col_1} / {col_2}  "
              f"[{total} reads, peaks_raw={peaks}]", file=sys.stderr)

    make_plot(results, f"{args.output}.size_distribution.png",
              args.basename, args.min_length, args.max_length, args.matrix_scale)

    write_count_matrix(results, f"{args.output}.count_matrix.tsv",
                       args.min_length, args.max_length, args.matrix_scale)

    output_tsv = f"{args.output}.genotype.tsv"
    ensure_parent_dir(output_tsv)
    header_cols = ["sample"]
    for m in MARKERS:
        header_cols += [f"{m}_1", f"{m}_2"]
    header_cols += ["flags"]

    row_cols = [os.path.basename(args.basename)]
    flag_notes = []
    for res in results:
        c1, c2 = format_allele(res["dominant_ce"], res["secondary_ce"], res["flag"],
                               res.get("h2_ratio", 0.0), args.soft_bracket_threshold)
        row_cols += [c1, c2]
        if res["flag"] not in ("heterozygous", "homozygous", "soft"):
            flag_notes.append(f"{res['marker']}:{res['flag']}")
    row_cols.append(";".join(flag_notes) if flag_notes else "OK")

    with open(output_tsv, "wt") as fh:
        fh.write("\t".join(header_cols) + "\n")
        fh.write("\t".join(row_cols) + "\n")
    print(f"INFO  Genotype TSV written to {output_tsv}", file=sys.stderr)

    print("\n=== GENOTYPE PROFILE ===")
    print(f"Sample: {os.path.basename(args.basename)}")
    print(f"{'Marker':<12} {'_1 (small)':>12} {'_2 (large)':>12}  {'Status'}")
    print("-" * 55)
    for res in results:
        c1, c2 = format_allele(res["dominant_ce"], res["secondary_ce"], res["flag"],
                               res.get("h2_ratio", 0.0), args.soft_bracket_threshold)
        print(f"{res['marker']:<12} {c1:>12} {c2:>12}  {res['flag']}")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
