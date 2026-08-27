#!/usr/bin/env python3
"""Summarize motif-length-1 TRF BED records by exact homopolymer motif.

Expected input is the BED-like format written by merge_window_trf_dat_to_bed.py:

    chrom  bed_start  bed_end  trf_start  trf_end  period  ...  consensus_pattern

By default, the script keeps records whose period / motif length in column 6 is
1, then counts the consensus motif in column 17. It writes:

1. counts by motif for all records, chromosome, haplotype, chromosome/haplotype,
   and input file;
2. repeat-length distributions for each motif, using BED end - start;
3. the input BED file list.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


BASE_ORDER = ("A", "C", "G", "T")
TOTAL_KEY = "__total__"
HAP_RE = re.compile(r"(?:^|_)(MATERNAL|PATERNAL)(?:_|\.|$)", re.IGNORECASE)


@dataclass(frozen=True)
class Motif1Record:
    source: str
    chrom: str
    haplotype: str
    motif: str
    repeat_len_bp: int


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def iter_bed_files(input_dir: Path, recursive: bool) -> Iterable[Path]:
    patterns = ("*.bed", "*.bed.gz")
    for pattern in patterns:
        glob_pattern = f"**/{pattern}" if recursive else pattern
        yield from sorted(path for path in input_dir.glob(glob_pattern) if path.is_file())


def parse_haplotype(*values: str) -> str:
    for value in values:
        match = HAP_RE.search(value)
        if match:
            return match.group(1).upper()
    return "UNKNOWN"


def parse_int(value: str) -> int:
    return int(float(value))


def normalize_motif(value: str) -> str:
    motif = value.strip().upper()
    if not motif:
        return "EMPTY"
    return motif


def read_motif1_records(
    bed_path: Path,
    motif_len_col: int,
    motif_col: int,
    skip_header: bool,
) -> Iterable[Motif1Record]:
    motif_len_idx = motif_len_col - 1
    motif_idx = motif_col - 1
    start_idx = 1
    end_idx = 2
    max_idx = max(motif_len_idx, motif_idx, start_idx, end_idx)

    with open_text(bed_path) as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if skip_header and line_no == 1 and (
                len(fields) <= end_idx or not fields[start_idx].lstrip("-").isdigit()
            ):
                continue
            if len(fields) <= max_idx:
                raise ValueError(
                    f"{bed_path}:{line_no} has {len(fields)} columns; "
                    f"need at least column {max_idx + 1}"
                )
            try:
                motif_len = parse_int(fields[motif_len_idx])
                start0 = parse_int(fields[start_idx])
                end0 = parse_int(fields[end_idx])
            except ValueError as exc:
                raise ValueError(f"{bed_path}:{line_no} has non-numeric coordinate or motif length") from exc
            if motif_len != 1:
                continue

            chrom = fields[0]
            yield Motif1Record(
                source=bed_path.name,
                chrom=chrom,
                haplotype=parse_haplotype(chrom, bed_path.name),
                motif=normalize_motif(fields[motif_idx]),
                repeat_len_bp=end0 - start0,
            )


def add_count(
    counts_by_group: dict[tuple[str, str], Counter[str]],
    group_name: str,
    group_value: str,
    motif: str,
) -> None:
    counts = counts_by_group.setdefault((group_name, group_value), Counter())
    counts[TOTAL_KEY] += 1
    counts[motif] += 1


def motif_sort_key(motif: str) -> tuple[int, str]:
    if motif in BASE_ORDER:
        return (BASE_ORDER.index(motif), motif)
    return (len(BASE_ORDER), motif)


def write_counts_summary(
    path: Path,
    counts_by_group: dict[tuple[str, str], Counter[str]],
) -> None:
    motifs = sorted(
        {motif for counts in counts_by_group.values() for motif in counts if motif != TOTAL_KEY},
        key=motif_sort_key,
    )

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["group", "value", "total_motif1", "motif", "count", "fraction"])
        for group_name, group_value in sorted(counts_by_group):
            counts = counts_by_group[(group_name, group_value)]
            total = counts[TOTAL_KEY]
            for motif in motifs:
                count = counts.get(motif, 0)
                fraction = count / total if total else 0.0
                writer.writerow([group_name, group_value, total, motif, count, f"{fraction:.6g}"])


def write_length_distribution(
    path: Path,
    length_counts: dict[tuple[str, int], int],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["motif", "repeat_len_bp", "count"])
        for motif, repeat_len_bp in sorted(length_counts, key=lambda item: (motif_sort_key(item[0]), item[1])):
            writer.writerow([motif, repeat_len_bp, length_counts[(motif, repeat_len_bp)]])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Count motif-length-1 / homopolymer records in TRF BED files."
    )
    parser.add_argument("input_dir", type=Path, help="Directory containing TRF BED files.")
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        help="Output prefix, e.g. results/trf_motif1_homopolymers",
    )
    parser.add_argument(
        "--motif-len-col",
        type=int,
        default=6,
        help="1-based column containing TRF period / motif length. Default: 6.",
    )
    parser.add_argument(
        "--motif-col",
        type=int,
        default=17,
        help="1-based column containing TRF consensus pattern / motif. Default: 17.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search for BED files recursively under input_dir.",
    )
    parser.add_argument(
        "--skip-header",
        action="store_true",
        help="Skip a first-line header if the BED start column is non-numeric.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {args.input_dir}")
    if args.motif_len_col < 1 or args.motif_col < 1:
        raise SystemExit("--motif-len-col and --motif-col must be >= 1")

    bed_files = list(iter_bed_files(args.input_dir, args.recursive))
    if not bed_files:
        raise SystemExit(f"No .bed or .bed.gz files found in {args.input_dir}")

    counts_by_group: dict[tuple[str, str], Counter[str]] = {}
    length_counts: dict[tuple[str, int], int] = Counter()
    total_records = 0

    for bed_file in bed_files:
        for record in read_motif1_records(
            bed_file,
            motif_len_col=args.motif_len_col,
            motif_col=args.motif_col,
            skip_header=args.skip_header,
        ):
            total_records += 1
            add_count(counts_by_group, "all", "all", record.motif)
            add_count(counts_by_group, "chromosome", record.chrom, record.motif)
            add_count(counts_by_group, "haplotype", record.haplotype, record.motif)
            add_count(
                counts_by_group,
                "chromosome_haplotype",
                f"{record.chrom}|{record.haplotype}",
                record.motif,
            )
            add_count(counts_by_group, "file", record.source, record.motif)
            length_counts[(record.motif, record.repeat_len_bp)] += 1

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    counts_path = prefix.with_suffix(".motif1_counts.summary.tsv")
    lengths_path = prefix.with_suffix(".motif1_repeat_length_distribution.tsv")
    inputs_path = prefix.with_suffix(".input_files.txt")

    write_counts_summary(counts_path, counts_by_group)
    write_length_distribution(lengths_path, length_counts)
    with inputs_path.open("w", encoding="utf-8") as handle:
        for bed_file in bed_files:
            handle.write(str(bed_file) + "\n")

    print(f"Read {len(bed_files)} BED files")
    print(f"Found {total_records} motif-length-1 records")
    print(f"Wrote {counts_path}")
    print(f"Wrote {lengths_path}")
    print(f"Wrote {inputs_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
