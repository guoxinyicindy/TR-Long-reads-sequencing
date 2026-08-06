#!/usr/bin/env python3
"""
Filter read-reference Levenshtein TSV outputs to loci present in a BED file.

The BED and read_ref_lev TSV files are matched by (chrom, start, end), using
standard 0-based half-open BED coordinates. Summary total rows are not kept,
because they describe the unfiltered input.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


LocusKey = tuple[str, str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Keep only read_ref_lev summary/per-read rows whose chrom/start/end "
            "are present in an input BED file."
        )
    )
    parser.add_argument("--bed", required=True, help="BED file of loci to keep.")
    parser.add_argument(
        "--summary-tsv",
        required=True,
        help="Input read_ref_lev.ont_pacbio.summary.tsv before filtering.",
    )
    parser.add_argument(
        "--per-read-tsv",
        help="Optional input read_ref_lev.ont_pacbio.per_read.tsv before filtering.",
    )
    parser.add_argument(
        "--out-summary-tsv",
        required=True,
        help="Filtered summary TSV output.",
    )
    parser.add_argument(
        "--out-per-read-tsv",
        help="Filtered per-read TSV output. Required if --per-read-tsv is provided.",
    )
    return parser.parse_args()


def load_bed_keys(bed_path: str) -> set[LocusKey]:
    keys: set[LocusKey] = set()
    with open(bed_path, "r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise SystemExit(f"{bed_path}:{line_no} has fewer than 3 BED columns")
            chrom, start, end = fields[0], fields[1], fields[2]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError as exc:
                raise SystemExit(f"{bed_path}:{line_no} has non-integer start/end") from exc
            if end_i <= start_i:
                raise SystemExit(f"{bed_path}:{line_no} has end <= start")
            keys.add((chrom, str(start_i), str(end_i)))
    if not keys:
        raise SystemExit(f"No loci found in BED: {bed_path}")
    return keys


def require_columns(path: str, fieldnames: list[str] | None, required: set[str]) -> None:
    missing = required.difference(fieldnames or [])
    if missing:
        raise SystemExit(f"{path} is missing columns: {', '.join(sorted(missing))}")


def row_key(row: dict[str, str]) -> LocusKey:
    return row["chrom"], row["start"], row["end"]


def filter_summary(summary_tsv: str, out_tsv: str, keep_keys: set[LocusKey]) -> tuple[int, int]:
    in_rows = 0
    kept_rows = 0
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        open(summary_tsv, "r", encoding="utf-8", newline="") as in_handle,
        open(out_path, "w", encoding="utf-8", newline="") as out_handle,
    ):
        reader = csv.DictReader(in_handle, delimiter="\t")
        require_columns(summary_tsv, reader.fieldnames, {"scope", "chrom", "start", "end"})
        writer = csv.DictWriter(
            out_handle,
            fieldnames=reader.fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            in_rows += 1
            if row["scope"] != "locus":
                continue
            if row_key(row) not in keep_keys:
                continue
            writer.writerow(row)
            kept_rows += 1
    return in_rows, kept_rows


def filter_per_read(per_read_tsv: str, out_tsv: str, keep_keys: set[LocusKey]) -> tuple[int, int]:
    in_rows = 0
    kept_rows = 0
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        open(per_read_tsv, "r", encoding="utf-8", newline="") as in_handle,
        open(out_path, "w", encoding="utf-8", newline="") as out_handle,
    ):
        reader = csv.DictReader(in_handle, delimiter="\t")
        require_columns(per_read_tsv, reader.fieldnames, {"chrom", "start", "end"})
        writer = csv.DictWriter(
            out_handle,
            fieldnames=reader.fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            in_rows += 1
            if row_key(row) not in keep_keys:
                continue
            writer.writerow(row)
            kept_rows += 1
    return in_rows, kept_rows


def main() -> None:
    args = parse_args()
    if args.per_read_tsv and not args.out_per_read_tsv:
        raise SystemExit("--out-per-read-tsv is required when --per-read-tsv is provided")

    keep_keys = load_bed_keys(args.bed)
    summary_in, summary_kept = filter_summary(
        args.summary_tsv,
        args.out_summary_tsv,
        keep_keys,
    )
    print(
        (
            f"[filter-read-ref-lev] BED loci={len(keep_keys)} "
            f"summary_rows={summary_in} summary_kept={summary_kept} "
            f"out_summary={args.out_summary_tsv}"
        ),
        file=sys.stderr,
    )

    if args.per_read_tsv:
        per_read_in, per_read_kept = filter_per_read(
            args.per_read_tsv,
            args.out_per_read_tsv,
            keep_keys,
        )
        print(
            (
                f"[filter-read-ref-lev] per_read_rows={per_read_in} "
                f"per_read_kept={per_read_kept} "
                f"out_per_read={args.out_per_read_tsv}"
            ),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
