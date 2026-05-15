#!/usr/bin/env python3
"""
kirk.py

Extract PCR products (including primer sequences) from Oxford Nanopore reads
for one or more predefined SSR loci.

If --locus is not specified, all predefined loci are processed.

For each processed marker, the script writes:
  - <output_basename>.<MARKER>.fastq.gz
  - <output_basename>.<MARKER>.summary.tsv   (if --summary true)

The extracted amplicon includes both primers.

Usage examples
--------------
Process all loci from one FASTQ:
    python kirk.py \
        --input sample.fastq.gz \
        --output results/sample1

Process only one locus:
    python kirk.py \
        --locus VrZAG62 \
        --input sample.fastq.gz \
        --output results/sample1

Disable summary files:
    python kirk.py \
        --input sample.fastq.gz \
        --output results/sample1 \
        --summary false
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class Locus:
    name: str
    forward: str
    reverse: str
    min_len: int
    max_len: int


LOCI: Dict[str, Locus] = {
    "VVMD28": Locus(
        name="VVMD28",
        forward="AACAATTCAATGAAAAGAGAGAGAGAGA",
        reverse="TCATCAATTTCGTATCTCTATTTGCTG",
        min_len=210,
        max_len=279,
    ),
    "VVMD5": Locus(
        name="VVMD5",
        forward="CTAGAGCTACGCCAATCCAA",
        reverse="TATACCAAAAATCATATTCCTAAA",
        min_len=215,
        max_len=270,
    ),
    "VVS2": Locus(
        name="VVS2",
        forward="CAGCCCGTAAATGTATCCATC",
        reverse="AAATTCAAAATTCTAATTCAACTGG",
        min_len=123,
        max_len=165,
    ),
    "VrZAG79": Locus(
        name="VrZAG79",
        forward="AGATTGTGGAGGAGGGAACAAACCG",
        reverse="TGCCCCCATTTTCAAACTCCCTTCC",
        min_len=230,
        max_len=270,
    ),
    "VrZAG62": Locus(
        name="VrZAG62",
        forward="GGTGAAATGGGCACCGAACACACGC",
        reverse="CCATGTCTCTCCTCAGCTTCTCAGC",
        min_len=181,
        max_len=220,
    ),
    "VVMD27": Locus(
        name="VVMD27",
        forward="GTACCAGATCTGAATACATCCGTAAGT",
        reverse="ACGGGTATAGAGCAAACGGTGT",
        min_len=165,
        max_len=210,
    ),
    "VVMD7": Locus(
        name="VVMD7",
        forward="AGAGTTGCGGAGAACAGGAT",
        reverse="CGAACCTTCACACGCTTGAT",
        min_len=180,
        max_len=270,
    ),
    "VVMD25": Locus(
        name="VVMD25",
        forward="TTCCGTTAAAGCAAAAGAAAAAGG",
        reverse="TTGGATTTGAAATTTATTGAGGGG",
        min_len=229,
        max_len=275,
    ),
    "VVMD32": Locus(
        name="VVMD32",
        forward="TATGATTTTTTAGGGGGGTGAGG",
        reverse="GGAAAGATGGGATGACTCGC",
        min_len=239,
        max_len=273,
    ),
}


@dataclass
class Match:
    start: int
    end: int
    mismatches: int


@dataclass
class Candidate:
    orientation: str
    f_match: Match
    r_match: Match
    total_mismatches: int
    product_len: int
    midpoint_distance: float
    extracted_seq: str
    extracted_qual: str


def str2bool(value: str) -> bool:
    value_lower = value.lower()
    if value_lower in {"true", "t", "yes", "y", "1"}:
        return True
    if value_lower in {"false", "f", "no", "n", "0"}:
        return False
    raise argparse.ArgumentTypeError(
        f"Invalid boolean value: {value!r}. Use true/false, yes/no, 1/0."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract PCR products for one or all predefined SSR loci from a FASTQ/FASTQ.GZ file. "
            "The extracted amplicon includes both primers."
        )
    )
    parser.add_argument(
        "--locus",
        default=None,
        choices=sorted(LOCI.keys()),
        help="Optional locus name to extract. If omitted, all loci are processed.",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input FASTQ file (can be .fastq, .fq, .fastq.gz, or .fq.gz).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help=(
            "Output basename. Files will be written as "
            "<basename>.<MARKER>.fastq.gz and optionally <basename>.<MARKER>.summary.tsv"
        ),
    )
    parser.add_argument(
        "--summary",
        type=str2bool,
        default=True,
        help="Whether to write per-marker summary TSV files (default: true).",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=2,
        help="Maximum mismatches allowed per primer (default: 2). No indels are allowed.",
    )
    parser.add_argument(
        "--slack",
        type=int,
        default=20,
        help=(
            "Allowed extra slack added to the expected product length range on both sides "
            "(default: 20)."
        ),
    )
    parser.add_argument(
        "--min-read-length",
        type=int,
        default=0,
        help="Skip reads shorter than this length before searching (default: 0).",
    )
    return parser.parse_args()


def open_text(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def iter_fastq(path: str) -> Iterable[Tuple[str, str, str, str]]:
    with open_text(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Incomplete FASTQ record encountered in {path}")
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def count_mismatches(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def find_matches_no_indels(seq: str, primer: str, max_mismatches: int, seed_len: int = 8) -> List[Match]:
    matches: List[Match] = []
    p_len = len(primer)
    s_len = len(seq)
    if s_len < p_len:
        return matches

    primer = primer.upper()
    seq = seq.upper()

    # Use an exact seed from the 5' end of the primer to find candidate positions quickly
    seed = primer[:seed_len]
    start = 0
    seen = set()

    while True:
        pos = seq.find(seed, start)
        if pos == -1:
            break

        if pos not in seen and pos + p_len <= s_len:
            window = seq[pos:pos + p_len]
            mm = count_mismatches(window, primer)
            if mm <= max_mismatches:
                matches.append(Match(start=pos, end=pos + p_len, mismatches=mm))
            seen.add(pos)

        start = pos + 1

    return matches


def build_candidates_for_orientation(
    seq: str,
    qual: str,
    forward_primer: str,
    reverse_primer_rc: str,
    min_len: int,
    max_len: int,
    max_mismatches: int,
    orientation_label: str,
) -> List[Candidate]:
    candidates: List[Candidate] = []
    f_matches = find_matches_no_indels(seq, forward_primer, max_mismatches, seed_len=8)
    r_matches = find_matches_no_indels(seq, reverse_primer_rc, max_mismatches, seed_len=8)

    midpoint = (min_len + max_len) / 2.0

    for f_match in f_matches:
        for r_match in r_matches:
            if r_match.start <= f_match.start:
                continue

            product_len = r_match.end - f_match.start
            if product_len < min_len or product_len > max_len:
                continue

            extracted_seq = seq[f_match.start:r_match.end]
            extracted_qual = qual[f_match.start:r_match.end]

            candidates.append(
                Candidate(
                    orientation=orientation_label,
                    f_match=f_match,
                    r_match=r_match,
                    total_mismatches=f_match.mismatches + r_match.mismatches,
                    product_len=product_len,
                    midpoint_distance=abs(product_len - midpoint),
                    extracted_seq=extracted_seq,
                    extracted_qual=extracted_qual,
                )
            )

    return candidates


def choose_best_candidate(candidates: List[Candidate]) -> Optional[Candidate]:
    if not candidates:
        return None
    return sorted(
        candidates,
        key=lambda c: (
            c.total_mismatches,
            c.midpoint_distance,
            c.product_len,
            c.f_match.start,
        ),
    )[0]


def prepare_output_path(path: str) -> None:
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)


def process_locus(
    input_path: str,
    output_basename: str,
    locus: Locus,
    max_mismatches: int,
    slack: int,
    min_read_length: int,
    write_summary: bool,
) -> Dict[str, object]:
    output_fastq = f"{output_basename}.{locus.name}.fastq.gz"
    summary_path = f"{output_basename}.{locus.name}.summary.tsv"

    prepare_output_path(output_fastq)
    if write_summary:
        prepare_output_path(summary_path)

    forward = locus.forward.upper()
    reverse_rc = revcomp(locus.reverse.upper())
    min_len = locus.min_len - slack
    max_len = locus.max_len + slack

    total_reads = 0
    extracted_reads = 0
    nohit_reads = 0
    too_short_reads = 0

    summary_handle = None
    try:
        if write_summary:
            summary_handle = open(summary_path, "wt")
            summary_handle.write(
                "\t".join(
                    [
                        "read_id",
                        "locus",
                        "status",
                        "orientation",
                        "product_start",
                        "product_end",
                        "product_length",
                        "forward_primer_mismatches",
                        "reverse_primer_mismatches",
                        "total_primer_mismatches",
                        "input_read_length",
                    ]
                ) + "\n"
            )

        with gzip.open(output_fastq, "wt") as out_fq:
            for header, seq, plus, qual in iter_fastq(input_path):
                total_reads += 1
                read_id = header[1:].split()[0] if header.startswith("@") else header.split()[0]

                if len(seq) < min_read_length:
                    too_short_reads += 1
                    nohit_reads += 1
                    if write_summary and summary_handle is not None:
                        summary_handle.write(
                            f"{read_id}\t{locus.name}\ttoo_short\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t{len(seq)}\n"
                        )
                    continue

                candidates: List[Candidate] = []

                candidates.extend(
                    build_candidates_for_orientation(
                        seq=seq,
                        qual=qual,
                        forward_primer=forward,
                        reverse_primer_rc=reverse_rc,
                        min_len=min_len,
                        max_len=max_len,
                        max_mismatches=max_mismatches,
                        orientation_label="forward_read",
                    )
                )

                seq_rc = revcomp(seq)
                qual_rc = qual[::-1]
                candidates.extend(
                    build_candidates_for_orientation(
                        seq=seq_rc,
                        qual=qual_rc,
                        forward_primer=forward,
                        reverse_primer_rc=reverse_rc,
                        min_len=min_len,
                        max_len=max_len,
                        max_mismatches=max_mismatches,
                        orientation_label="reverse_read",
                    )
                )

                best = choose_best_candidate(candidates)

                if best is None:
                    nohit_reads += 1
                    if write_summary and summary_handle is not None:
                        summary_handle.write(
                            f"{read_id}\t{locus.name}\tno_valid_product\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t{len(seq)}\n"
                        )
                    continue

                extracted_reads += 1

                new_header = (
                    f"@{read_id} locus={locus.name} orientation={best.orientation} "
                    f"len={best.product_len} mismatches={best.total_mismatches}"
                )
                out_fq.write(new_header + "\n")
                out_fq.write(best.extracted_seq + "\n")
                out_fq.write("+\n")
                out_fq.write(best.extracted_qual + "\n")

                if write_summary and summary_handle is not None:
                    summary_handle.write(
                        "\t".join(
                            [
                                read_id,
                                locus.name,
                                "extracted",
                                best.orientation,
                                str(best.f_match.start),
                                str(best.r_match.end),
                                str(best.product_len),
                                str(best.f_match.mismatches),
                                str(best.r_match.mismatches),
                                str(best.total_mismatches),
                                str(len(seq)),
                            ]
                        ) + "\n"
                    )
    finally:
        if summary_handle is not None:
            summary_handle.close()

    return {
        "locus": locus.name,
        "output_fastq": output_fastq,
        "summary_path": summary_path if write_summary else None,
        "total_reads": total_reads,
        "extracted_reads": extracted_reads,
        "nohit_reads": nohit_reads,
        "too_short_reads": too_short_reads,
    }


def main() -> int:
    args = parse_args()

    loci_to_process: List[Locus]
    if args.locus is None:
        loci_to_process = [LOCI[name] for name in sorted(LOCI.keys())]
    else:
        loci_to_process = [LOCI[args.locus]]

    print(
        f"[INFO] Processing {len(loci_to_process)} marker(s): "
        + ", ".join(locus.name for locus in loci_to_process),
        file=sys.stderr,
    )

    all_results = []
    for locus in loci_to_process:
        print(f"[INFO] Processing locus {locus.name}", file=sys.stderr)
        result = process_locus(
            input_path=args.input,
            output_basename=args.output,
            locus=locus,
            max_mismatches=args.max_mismatches,
            slack=args.slack,
            min_read_length=args.min_read_length,
            write_summary=args.summary,
        )
        all_results.append(result)

        print(
            f"[INFO]   Input reads: {result['total_reads']}\n"
            f"[INFO]   Extracted PCR products: {result['extracted_reads']}\n"
            f"[INFO]   Reads without valid product: {result['nohit_reads']}\n"
            f"[INFO]   Output FASTQ: {result['output_fastq']}"
            + (
                f"\n[INFO]   Summary TSV: {result['summary_path']}"
                if result["summary_path"] is not None
                else ""
            ),
            file=sys.stderr,
        )

    print("[INFO] Done.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
