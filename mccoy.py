#!/usr/bin/env python3
"""
mccoy.py

Compares a sample genotype profile (produced by ssr_genotype_sample.py)
against the VIVC SSR reference database to identify the best-matching variety.

Matching strategy:
  - Only clearly heterozygous loci (both alleles present, no brackets/tildes)
    are used as anchor loci for initial filtering.
  - A candidate is only reported if at least --min-anchor-matches anchor loci
    have both alleles matching within tolerance (default: 3).
  - Passing candidates are then scored across all callable loci for ranking.

Reference column names (as they appear in vivc_ssr-profiles.tsv):
  VVS2_1/2, VVMD5_1/2, VVMD7_1/2, VVMD25_1/2, VVMD27_1/2,
  VVMD28_1/2, VVMD32_1/2, VRZAG62_1/2, VRZAG79_1/2

Usage:
  python mccoy.py \
      --genotype  01.genotype.tsv \
      --reference vivc_ssr-profiles.tsv \
      --output    01 \
      [--top-n 10] \
      [--match-tolerance 1] \
      [--min-anchor-matches 3]
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, List, Optional, Tuple

# Internal marker names used throughout this script
MARKERS: List[str] = [
    "VVS2", "VVMD5", "VVMD7", "VVMD25",
    "VVMD27", "VVMD28", "VVMD32", "VrZAG62", "VrZAG79",
]

# Map from internal marker name -> column prefix in the reference TSV
# The reference uses uppercase VRZAG62 / VRZAG79 (no lowercase 'r')
REF_COL: Dict[str, str] = {
    "VVS2":    "VVS2",
    "VVMD5":   "VVMD5",
    "VVMD7":   "VVMD7",
    "VVMD25":  "VVMD25",
    "VVMD27":  "VVMD27",
    "VVMD28":  "VVMD28",
    "VVMD32":  "VVMD32",
    "VrZAG62": "VRZAG62",
    "VrZAG79": "VRZAG79",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Match a sample SSR genotype against the VIVC reference database."
    )
    p.add_argument("--genotype", required=True,
                   help="Genotype TSV produced by ssr_genotype_sample.py.")
    p.add_argument("--reference", required=True,
                   help="VIVC SSR profiles TSV (vivc_ssr-profiles.tsv).")
    p.add_argument("--output", required=True,
                   help="Output basename; writes <output>.variety_match.tsv.")
    p.add_argument("--top-n", type=int, default=10,
                   help="Number of top candidates to report (default: 10).")
    p.add_argument("--match-tolerance", type=int, default=1,
                   help="+-bp tolerance when comparing alleles (default: 1).")
    p.add_argument("--min-anchor-matches", type=int, default=3,
                   help="Minimum number of heterozygous anchor loci that must "
                        "match (both alleles) for a candidate to be reported "
                        "(default: 3).")
    return p.parse_args()


def parse_allele(s: str) -> Optional[int]:
    s = s.strip().lstrip("(~").rstrip(")")
    if not s:
        return None
    try:
        return int(s)
    except ValueError:
        return None


def is_bracketed(s: str) -> bool:
    s = s.strip()
    return s.startswith("(") or s.startswith("~")


def alleles_match(a: Optional[int], b: Optional[int], tol: int) -> bool:
    if a is None or b is None:
        return False
    return abs(a - b) <= tol


def score_locus(
    s1: Optional[int], s2: Optional[int],
    r1: Optional[int], r2: Optional[int],
    tol: int,
) -> Tuple[bool, bool]:
    """Returns (is_comparable, both_alleles_match)."""
    if r1 is None or s1 is None:
        return False, False

    ref_hom  = r2 is None or r1 == r2
    samp_hom = s2 is None

    if samp_hom and ref_hom:
        return True, alleles_match(s1, r1, tol)
    if samp_hom and not ref_hom:
        return True, alleles_match(s1, r1, tol) or alleles_match(s1, r2, tol)
    if not samp_hom and ref_hom:
        return True, alleles_match(s1, r1, tol) and alleles_match(s2, r1, tol)
    # both heterozygous
    return True, (
        (alleles_match(s1, r1, tol) and alleles_match(s2, r2, tol)) or
        (alleles_match(s1, r2, tol) and alleles_match(s2, r1, tol))
    )


def load_genotype(
    path: str,
) -> Tuple[str, Dict[str, Tuple[Optional[int], Optional[int]]], List[str]]:
    """
    Returns (sample_name, genotype_dict, anchor_markers).
    anchor_markers = markers with both alleles clearly called (no brackets/tildes).
    """
    with open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty genotype file: {path}")
    row = rows[0]
    sample = row.get("sample", os.path.basename(path))

    genotype: Dict[str, Tuple[Optional[int], Optional[int]]] = {}
    anchor_markers: List[str] = []

    for marker in MARKERS:
        raw_1 = row.get(f"{marker}_1", "")
        raw_2 = row.get(f"{marker}_2", "")
        a1 = parse_allele(raw_1)
        a2 = parse_allele(raw_2)
        genotype[marker] = (a1, a2)
        if (a1 is not None and a2 is not None
                and not is_bracketed(raw_1)
                and not is_bracketed(raw_2)):
            anchor_markers.append(marker)

    return sample, genotype, anchor_markers


def load_reference(path: str) -> List[dict]:
    """
    Load the VIVC reference TSV. Column names for the markers follow the
    pattern <REF_COL[marker]>_1 and <REF_COL[marker]>_2.
    """
    entries = []
    with open(path, "rt", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # Normalise header keys to uppercase for robust lookup
        fieldnames_upper = {k.upper(): k for k in (reader.fieldnames or [])}

        for row in reader:
            # Re-key the row with uppercase keys for lookup
            row_upper = {k.upper(): v for k, v in row.items()}

            entry: dict = {"name": row_upper.get("PRIME NAME", "").strip()}
            for marker in MARKERS:
                col = REF_COL[marker]
                raw_1 = row_upper.get(f"{col}_1", "")
                raw_2 = row_upper.get(f"{col}_2", "")
                a1 = parse_allele(raw_1)
                a2 = parse_allele(raw_2)
                # Ensure a1 <= a2
                if a1 is not None and a2 is not None and a1 > a2:
                    a1, a2 = a2, a1
                entry[marker] = (a1, a2)
            entries.append(entry)
    return entries


def evaluate_entry(
    sample_gt: Dict[str, Tuple[Optional[int], Optional[int]]],
    anchor_markers: List[str],
    ref_entry: dict,
    tol: int,
    min_anchor: int,
) -> Optional[Tuple[float, int, int, int, int]]:
    """
    Returns (score, full_matched, full_comparable, full_mismatched, anchor_matched)
    or None if the anchor threshold is not met.
    """
    anchor_matched = 0
    for marker in anchor_markers:
        s1, s2 = sample_gt[marker]
        r1, r2 = ref_entry[marker]
        comparable, both_match = score_locus(s1, s2, r1, r2, tol)
        if comparable and both_match:
            anchor_matched += 1

    if anchor_matched < min_anchor:
        return None

    full_matched = 0
    full_comparable = 0
    for marker in MARKERS:
        s1, s2 = sample_gt[marker]
        if s1 is None:
            continue
        r1, r2 = ref_entry[marker]
        comparable, both_match = score_locus(s1, s2, r1, r2, tol)
        if comparable:
            full_comparable += 1
            if both_match:
                full_matched += 1

    full_mismatched = full_comparable - full_matched
    score = full_matched / full_comparable if full_comparable > 0 else 0.0
    return score, full_matched, full_comparable, full_mismatched, anchor_matched


def main() -> int:
    args = parse_args()

    sample_name, sample_gt, anchor_markers = load_genotype(args.genotype)
    print(f"INFO  Sample: {sample_name}", file=sys.stderr)
    print(f"INFO  Anchor markers (clearly heterozygous): "
          f"{', '.join(anchor_markers) if anchor_markers else 'none'}", file=sys.stderr)
    print(f"INFO  Minimum anchor matches required: {args.min_anchor_matches}", file=sys.stderr)

    if len(anchor_markers) < args.min_anchor_matches:
        print(f"WARN  Only {len(anchor_markers)} anchor markers available; "
              f"minimum required is {args.min_anchor_matches}. "
              f"Consider lowering --min-anchor-matches.", file=sys.stderr)

    ref_entries = load_reference(args.reference)
    print(f"INFO  Loaded {len(ref_entries)} reference entries.", file=sys.stderr)

    scored = []
    for entry in ref_entries:
        result = evaluate_entry(
            sample_gt, anchor_markers, entry,
            args.match_tolerance, args.min_anchor_matches,
        )
        if result is not None:
            score, matched, comparable, mismatched, anchor_matched = result
            scored.append((score, matched, comparable, mismatched, anchor_matched, entry["name"]))

    scored.sort(key=lambda x: (-x[0], x[3], -x[2], -x[4]))
    top = scored[:args.top_n]

    output_tsv = f"{args.output}.variety_match.tsv"
    with open(output_tsv, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["rank", "variety_name", "score_pct",
                         "matched_loci", "comparable_loci", "mismatched_loci",
                         "anchor_loci_matched"])
        for rank, (score, matched, comparable, mismatched, anchor_matched, name) in enumerate(top, 1):
            writer.writerow([rank, name, f"{score*100:.1f}",
                             matched, comparable, mismatched, anchor_matched])
    print(f"INFO  Variety match TSV written to {output_tsv}", file=sys.stderr)

    print(f"\n=== VARIETY MATCH RESULTS ===")
    print(f"Sample: {sample_name}  |  "
          f"Anchor markers: {len(anchor_markers)}  |  "
          f"Candidates passing anchor filter: {len(scored)}")
    print(f"{'Rank':<6} {'Score':>8} {'Full match':>12} {'Mismatches':>11} "
          f"{'Anchors':>8}  Variety")
    print("-" * 72)
    for rank, (score, matched, comparable, mismatched, anchor_matched, name) in enumerate(top, 1):
        print(f"{rank:<6} {score*100:>7.1f}% {matched:>6}/{comparable:<5} "
              f"{mismatched:>11} {anchor_matched:>8}  {name}")
    print()

    if top:
        best_score = top[0][0]
        best_name  = top[0][5]
        if best_score == 1.0:
            print(f"\u2713  Perfect match: {best_name}")
        elif best_score >= 0.8:
            print(f"\u26a0  Best match at {best_score*100:.1f}%: {best_name}")
        else:
            print(f"\u2717  No high-confidence match (best: {best_score*100:.1f}% – {best_name}).")
    else:
        print(f"\u2717  No candidates passed the anchor filter "
              f"(>= {args.min_anchor_matches} heterozygous loci matching both alleles).")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
