#!/usr/bin/env python3
"""
Summarize per-locus read coverage by motif length bins.

Coverage is the n_reads value from read_ref_lev.summary.tsv locus rows. The
per-read TSV is used as a consistency check against those summary counts. By
default, locus_id is treated as the motif. If locus_id is not the motif, pass
--bed so BED column 4 can be used as the motif by chrom/start/end.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


BIN_ORDER = ["1bp", "2-10bp", "11-30bp", "31-50bp", "51-100bp", ">100bp"]
TOOL_ORDER = ["ONT", "PacBio", "realigned"]
LocusKey = tuple[str, str, str]
LocusToolKey = tuple[str, str, str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize average per-locus read coverage by motif length bins."
    )
    parser.add_argument(
        "--summary-tsv",
        required=True,
        help="Input read_ref_lev.summary.tsv with scope=locus rows and n_reads.",
    )
    parser.add_argument(
        "--per-read-tsv",
        required=True,
        help="Input read_ref_lev.per_read.tsv, used to check read counts.",
    )
    parser.add_argument(
        "--bed",
        help=(
            "Optional BED file used to assign motifs by chrom/start/end. "
            "BED column 4 is used by default."
        ),
    )
    parser.add_argument(
        "--motif-column",
        type=int,
        default=4,
        help="1-based motif column in --bed. Default: 4.",
    )
    parser.add_argument(
        "--out-tsv",
        required=True,
        help="Output motif-bin coverage summary TSV.",
    )
    parser.add_argument(
        "--out-length-1-10-tsv",
        help="Optional output TSV grouped by exact motif length from 1bp to 10bp.",
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


def fmt(value: float | int | None) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{value:.6f}"


def mean_or_none(values: list[int]) -> float | None:
    return statistics.fmean(values) if values else None


def median_or_none(values: list[int]) -> float | None:
    return statistics.median(values) if values else None


def min_or_none(values: list[int]) -> int | None:
    return min(values) if values else None


def max_or_none(values: list[int]) -> int | None:
    return max(values) if values else None


def require_columns(path: str, fieldnames: list[str] | None, required: set[str]) -> None:
    missing = required.difference(fieldnames or [])
    if missing:
        raise SystemExit(f"{path} is missing columns: {', '.join(sorted(missing))}")


def load_motif_map(bed_path: str | None, motif_column: int) -> dict[LocusKey, str]:
    if bed_path is None:
        return {}
    motif_index = motif_column - 1
    if motif_index < 0:
        raise SystemExit("--motif-column must be >= 1")

    motif_by_locus: dict[LocusKey, str] = {}
    with open(bed_path, "r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#") or row[0].startswith("track"):
                continue
            if len(row) < 3:
                raise SystemExit(f"{bed_path}:{line_no} has fewer than 3 BED columns")
            if len(row) <= motif_index:
                raise SystemExit(f"{bed_path}:{line_no} does not have motif column {motif_column}")
            motif_by_locus[(row[0], str(int(row[1])), str(int(row[2])))] = row[motif_index]
    return motif_by_locus


def row_locus_key(row: dict[str, str]) -> LocusKey:
    return row["chrom"], str(int(row["start"])), str(int(row["end"]))


def row_locus_tool_key(row: dict[str, str]) -> LocusToolKey:
    chrom, start, end = row_locus_key(row)
    return chrom, start, end, row["tool"]


def row_motif(row: dict[str, str], motif_by_locus: dict[LocusKey, str]) -> str | None:
    if motif_by_locus:
        return motif_by_locus.get(row_locus_key(row))
    return row.get("locus_id")


def load_per_read_counts(per_read_tsv: str) -> Counter[LocusToolKey]:
    counts: Counter[LocusToolKey] = Counter()
    with open(per_read_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "start", "end", "tool"}
        require_columns(per_read_tsv, reader.fieldnames, required)
        for row in reader:
            counts[row_locus_tool_key(row)] += 1
    return counts


def load_summary_coverages(
    summary_tsv: str,
    motif_by_locus: dict[LocusKey, str],
) -> tuple[
    dict[tuple[str, str], list[int]],
    dict[tuple[int, str], list[int]],
    dict[LocusToolKey, int],
    int,
]:
    coverages_by_bin_tool: dict[tuple[str, str], list[int]] = defaultdict(list)
    coverages_by_length_tool: dict[tuple[int, str], list[int]] = defaultdict(list)
    summary_counts: dict[LocusToolKey, int] = {}
    skipped_missing_motif = 0

    with open(summary_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"scope", "locus_id", "chrom", "start", "end", "tool", "n_reads"}
        require_columns(summary_tsv, reader.fieldnames, required)
        for row in reader:
            if row["scope"] != "locus":
                continue
            motif = row_motif(row, motif_by_locus)
            if motif_by_locus and motif is None:
                skipped_missing_motif += 1
                continue
            try:
                n_reads = int(row["n_reads"])
            except ValueError as exc:
                raise SystemExit(f"{summary_tsv} has non-integer n_reads: {row['n_reads']}") from exc
            if n_reads < 0:
                raise SystemExit(f"{summary_tsv} has negative n_reads: {n_reads}")
            tool = row["tool"]
            motif_value = motif or row["locus_id"]
            summary_counts[row_locus_tool_key(row)] = n_reads
            coverages_by_bin_tool[(motif_len_bin(motif_value), tool)].append(n_reads)
            motif_length = len(motif_value)
            if 1 <= motif_length <= 10:
                coverages_by_length_tool[(motif_length, tool)].append(n_reads)
    return coverages_by_bin_tool, coverages_by_length_tool, summary_counts, skipped_missing_motif


def tool_sort_key(tool: str) -> tuple[int, str]:
    try:
        return TOOL_ORDER.index(tool), tool
    except ValueError:
        return len(TOOL_ORDER), tool


def count_mismatched_loci(
    summary_counts: dict[LocusToolKey, int],
    per_read_counts: Counter[LocusToolKey],
) -> int:
    mismatches = 0
    for key, summary_count in summary_counts.items():
        if per_read_counts.get(key, 0) != summary_count:
            mismatches += 1
    return mismatches


def write_output(
    coverages_by_bin_tool: dict[tuple[str, str], list[int]],
    out_tsv: str,
) -> None:
    tools = sorted({tool for _, tool in coverages_by_bin_tool}, key=tool_sort_key)
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "motif_bin",
        "tool",
        "n_loci",
        "n_loci_with_reads",
        "n_loci_without_reads",
        "total_reads",
        "mean_reads_per_locus",
        "median_reads_per_locus",
        "min_reads_per_locus",
        "max_reads_per_locus",
        "mean_reads_per_covered_locus",
        "median_reads_per_covered_locus",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for motif_bin in BIN_ORDER:
            for tool in tools:
                values = coverages_by_bin_tool.get((motif_bin, tool), [])
                covered_values = [value for value in values if value > 0]
                writer.writerow(
                    [
                        motif_bin,
                        tool,
                        len(values),
                        len(covered_values),
                        len(values) - len(covered_values),
                        sum(values),
                        fmt(mean_or_none(values)),
                        fmt(median_or_none(values)),
                        fmt(min_or_none(values)),
                        fmt(max_or_none(values)),
                        fmt(mean_or_none(covered_values)),
                        fmt(median_or_none(covered_values)),
                    ]
                )


def write_length_1_10_output(
    coverages_by_length_tool: dict[tuple[int, str], list[int]],
    out_tsv: str,
) -> None:
    tools = sorted({tool for _, tool in coverages_by_length_tool}, key=tool_sort_key)
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "motif_length",
        "motif_length_label",
        "tool",
        "n_loci",
        "n_loci_with_reads",
        "n_loci_without_reads",
        "total_reads",
        "mean_reads_per_locus",
        "median_reads_per_locus",
        "min_reads_per_locus",
        "max_reads_per_locus",
        "mean_reads_per_covered_locus",
        "median_reads_per_covered_locus",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for motif_length in range(1, 11):
            for tool in tools:
                values = coverages_by_length_tool.get((motif_length, tool), [])
                covered_values = [value for value in values if value > 0]
                writer.writerow(
                    [
                        motif_length,
                        f"{motif_length}bp",
                        tool,
                        len(values),
                        len(covered_values),
                        len(values) - len(covered_values),
                        sum(values),
                        fmt(mean_or_none(values)),
                        fmt(median_or_none(values)),
                        fmt(min_or_none(values)),
                        fmt(max_or_none(values)),
                        fmt(mean_or_none(covered_values)),
                        fmt(median_or_none(covered_values)),
                    ]
                )


def main() -> None:
    args = parse_args()
    motif_by_locus = load_motif_map(args.bed, args.motif_column)
    per_read_counts = load_per_read_counts(args.per_read_tsv)
    (
        coverages_by_bin_tool,
        coverages_by_length_tool,
        summary_counts,
        skipped_missing_motif,
    ) = load_summary_coverages(
        args.summary_tsv,
        motif_by_locus,
    )
    write_output(coverages_by_bin_tool, args.out_tsv)
    if args.out_length_1_10_tsv:
        write_length_1_10_output(coverages_by_length_tool, args.out_length_1_10_tsv)

    mismatches = count_mismatched_loci(summary_counts, per_read_counts)
    if skipped_missing_motif:
        print(
            f"warning: skipped {skipped_missing_motif} summary rows without BED motif mapping",
            file=sys.stderr,
        )
    if mismatches:
        print(
            (
                "warning: "
                f"{mismatches} locus/tool rows have n_reads different from per-read row counts"
            ),
            file=sys.stderr,
        )
    print(
        (
            "[read-coverage-by-motif-bin] "
            f"locus_tool_rows={len(summary_counts)} "
            f"out={args.out_tsv}"
            + (
                f" out_length_1_10={args.out_length_1_10_tsv}"
                if args.out_length_1_10_tsv
                else ""
            )
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
