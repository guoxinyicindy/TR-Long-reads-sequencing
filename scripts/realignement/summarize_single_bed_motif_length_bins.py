#!/usr/bin/env python3
"""
Summarize motif length bins for a single BED file.

BED column 4 is expected to be motif by default.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


BIN_LABELS = ["1bp", "2-10bp", "11-30bp", "31-50bp", "51-100bp", ">100bp"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count motif length bins in a single BED file."
    )
    parser.add_argument(
        "bed",
        help="Input BED file.",
    )
    parser.add_argument(
        "--motif-column",
        type=int,
        default=4,
        help="1-based motif column in the BED file. Default: 4.",
    )
    parser.add_argument(
        "--out-tsv",
        help="Optional output TSV. Default: print to stdout.",
    )
    return parser.parse_args()


def motif_len_bin(motif: str) -> str:
    motif_len = len(motif)
    if motif_len == 1:
        return "1bp"
    if 2 <= motif_len <= 10:
        return "2-10bp"
    if 11 <= motif_len <= 30:
        return "11-30bp"
    if 31 <= motif_len <= 50:
        return "31-50bp"
    if 51 <= motif_len <= 100:
        return "51-100bp"
    return ">100bp"


def summarize_bed(path: Path, motif_index: int) -> tuple[Counter[str], int]:
    counts: Counter[str] = Counter()
    total = 0
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#") or row[0].startswith("track"):
                continue
            if len(row) <= motif_index:
                raise ValueError(
                    f"{path}:{line_number} has {len(row)} columns; "
                    f"motif column {motif_index + 1} is not available"
                )
            counts[motif_len_bin(row[motif_index])] += 1
            total += 1
    return counts, total


def write_summary(path: Path, motif_column: int, out_fh) -> None:
    motif_index = motif_column - 1
    if motif_index < 0:
        raise ValueError("--motif-column must be >= 1")
    if not path.exists():
        raise FileNotFoundError(f"BED file not found: {path}")

    counts, total = summarize_bed(path, motif_index)
    writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["bed_file", "bin", "count", "pct"])
    for label in BIN_LABELS:
        count = counts[label]
        pct = f"{100.0 * count / total:.6f}" if total else "NA"
        writer.writerow([path.name, label, count, pct])
    writer.writerow([path.name, "total", total, "100.000000" if total else "NA"])


def main() -> None:
    args = parse_args()
    path = Path(args.bed)

    if args.out_tsv:
        out_path = Path(args.out_tsv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8", newline="") as out_fh:
            write_summary(path, args.motif_column, out_fh)
        print(f"wrote motif length summary: {out_path}", file=sys.stderr)
    else:
        write_summary(path, args.motif_column, sys.stdout)


if __name__ == "__main__":
    main()
