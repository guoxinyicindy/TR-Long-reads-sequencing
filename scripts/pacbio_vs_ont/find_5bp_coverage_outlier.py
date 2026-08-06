#!/usr/bin/env python3
"""
Find a high-coverage exact-motif locus and print 5 bp mean read coverage after
excluding that locus.

The preferred path uses locus-level or per-read detail, because that can identify
the exact coordinates and remove the same locus from both ONT and PacBio. If the
current detail files do not contain the requested outlier, the script falls back
to the existing aggregate coverage TSVs and prints which assumptions are needed.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


DEFAULT_SUMMARY = "realignment/read_ref_lev.summary.tsv"
DEFAULT_PER_READ = "realignment/read_ref_lev.per_read.tsv"
DEFAULT_LENGTH_TSV = "pacbio_vs_ont/read_ref_lev.ont_pacbio.filtered.coverage_by_motif_length_1_10bp.tsv"
DEFAULT_EXACT_TSV = "pacbio_vs_ont/read_ref_lev.ont_pacbio.coverage_by_exact_5bp_motif.tsv"
DEFAULT_TOOLS = ("ONT", "PacBio")


@dataclass(frozen=True)
class LocusRow:
    motif: str
    chrom: str
    start: str
    end: str
    tool: str
    n_reads: int

    @property
    def key(self) -> tuple[str, str, str]:
        return self.chrom, self.start, self.end

    @property
    def coord(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Locate an exact-motif coverage outlier and print adjusted 5 bp mean coverage."
    )
    parser.add_argument("--summary-tsv", default=DEFAULT_SUMMARY)
    parser.add_argument(
        "--per-read-tsv",
        default=DEFAULT_PER_READ,
        help=(
            "Optional per-read TSV used as a secondary check. If the file does not "
            "exist, the script skips this check."
        ),
    )
    parser.add_argument("--coverage-by-motif-length-tsv", default=DEFAULT_LENGTH_TSV)
    parser.add_argument("--coverage-by-exact-motif-tsv", default=DEFAULT_EXACT_TSV)
    parser.add_argument("--motif", default="AATGG")
    parser.add_argument("--motif-length", type=int, default=5)
    parser.add_argument("--tool", default="PacBio")
    parser.add_argument("--coverage", type=int, default=25468)
    parser.add_argument("--tools", nargs="+", default=list(DEFAULT_TOOLS))
    parser.add_argument(
        "--exclude-coverage-above",
        type=int,
        help=(
            "Filter out an entire locus if any selected tool has n_reads greater "
            "than this value, then print the filtered 5 bp mean read coverage."
        ),
    )
    return parser.parse_args()


def require_columns(path: str, fieldnames: list[str] | None, required: set[str]) -> None:
    missing = required.difference(fieldnames or [])
    if missing:
        raise SystemExit(f"{path} is missing columns: {', '.join(sorted(missing))}")


def read_summary_rows(path: str, tools: set[str], motif_length: int) -> list[LocusRow]:
    rows: list[LocusRow] = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames, {"scope", "locus_id", "chrom", "start", "end", "tool", "n_reads"})
        for row in reader:
            motif = row["locus_id"].upper()
            if row["scope"] != "locus" or row["tool"] not in tools or len(motif) != motif_length:
                continue
            rows.append(
                LocusRow(
                    motif=motif,
                    chrom=row["chrom"],
                    start=str(int(row["start"])),
                    end=str(int(row["end"])),
                    tool=row["tool"],
                    n_reads=int(row["n_reads"]),
                )
            )
    return rows


def read_per_read_counts(path: str, tools: set[str], motif_length: int) -> list[LocusRow]:
    counts: Counter[tuple[str, str, str, str, str]] = Counter()
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames, {"locus_id", "chrom", "start", "end", "tool"})
        for row in reader:
            motif = row["locus_id"].upper()
            if row["tool"] not in tools or len(motif) != motif_length:
                continue
            key = motif, row["chrom"], str(int(row["start"])), str(int(row["end"])), row["tool"]
            counts[key] += 1

    return [
        LocusRow(motif=motif, chrom=chrom, start=start, end=end, tool=tool, n_reads=count)
        for (motif, chrom, start, end, tool), count in counts.items()
    ]


def find_target(rows: list[LocusRow], motif: str, tool: str, coverage: int) -> LocusRow | None:
    motif = motif.upper()
    for row in rows:
        if row.motif == motif and row.tool == tool and row.n_reads == coverage:
            return row
    return None


def print_available_motif_rows(label: str, rows: list[LocusRow], motif: str) -> None:
    motif_rows = [row for row in rows if row.motif == motif.upper()]
    if not motif_rows:
        print(f"{label}: no {motif.upper()} rows found")
        return

    print(f"{label}: {motif.upper()} rows found")
    print("motif\ttool\tcoord\tn_reads")
    for row in sorted(motif_rows, key=lambda r: (r.coord, r.tool)):
        print(f"{row.motif}\t{row.tool}\t{row.coord}\t{row.n_reads}")


def print_exact_adjusted_mean(rows: list[LocusRow], outlier: LocusRow, tools: list[str]) -> None:
    print()
    print("Adjusted 5 bp mean read coverage after excluding the located locus:")
    print("tool\tn_loci_before\ttotal_reads_before\tn_loci_after\ttotal_reads_after\tmean_after_excluding_locus")
    for tool in tools:
        values_before = [row.n_reads for row in rows if row.tool == tool]
        values_after = [row.n_reads for row in rows if row.tool == tool and row.key != outlier.key]
        mean_after = statistics.fmean(values_after) if values_after else None
        mean_text = "NA" if mean_after is None else f"{mean_after:.6f}"
        print(
            f"{tool}\t{len(values_before)}\t{sum(values_before)}\t"
            f"{len(values_after)}\t{sum(values_after)}\t{mean_text}"
        )


def excluded_locus_keys(
    rows: list[LocusRow],
    threshold: int,
) -> set[tuple[str, str, str]]:
    return {row.key for row in rows if row.n_reads > threshold}


def print_excluded_loci(rows: list[LocusRow], excluded_keys: set[tuple[str, str, str]], tools: list[str], threshold: int) -> None:
    print(f"Excluded 5 bp loci with any selected-tool n_reads > {threshold}:")
    if not excluded_keys:
        print("None")
        return

    print("motif\tcoord\texcluded_tools\tall_tool_n_reads")
    for key in sorted(excluded_keys, key=lambda item: (item[0], int(item[1]), int(item[2]))):
        locus_rows = [row for row in rows if row.key == key and row.tool in tools]
        if not locus_rows:
            continue
        motif = locus_rows[0].motif
        coord = locus_rows[0].coord
        locus_rows = sorted(locus_rows, key=lambda row: tools.index(row.tool) if row.tool in tools else len(tools))
        excluded_tools = ",".join(f"{row.tool}={row.n_reads}" for row in locus_rows if row.n_reads > threshold)
        all_counts = ",".join(f"{row.tool}={row.n_reads}" for row in locus_rows)
        print(f"{motif}\t{coord}\t{excluded_tools}\t{all_counts}")


def print_filtered_mean(rows: list[LocusRow], tools: list[str], threshold: int) -> None:
    excluded_keys = excluded_locus_keys(rows, threshold)
    print_excluded_loci(rows, excluded_keys, tools, threshold)
    print()
    print("Filtered 5 bp mean read coverage:")
    print("tool\tn_loci_before\ttotal_reads_before\tn_loci_after\ttotal_reads_after\tmean_after_filtering")
    for tool in tools:
        values_before = [row.n_reads for row in rows if row.tool == tool]
        values_after = [row.n_reads for row in rows if row.tool == tool and row.key not in excluded_keys]
        mean_after = statistics.fmean(values_after) if values_after else None
        mean_text = "NA" if mean_after is None else f"{mean_after:.6f}"
        print(
            f"{tool}\t{len(values_before)}\t{sum(values_before)}\t"
            f"{len(values_after)}\t{sum(values_after)}\t{mean_text}"
        )


def read_length_5_aggregates(path: str, tools: set[str]) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            path,
            reader.fieldnames,
            {"motif_length", "tool", "n_loci", "total_reads", "mean_reads_per_locus", "max_reads_per_locus"},
        )
        for row in reader:
            if row["motif_length"] == "5" and row["tool"] in tools:
                rows[row["tool"]] = row
    return rows


def read_exact_motif_aggregates(path: str, motif: str, tools: set[str]) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            path,
            reader.fieldnames,
            {"motif", "tool", "n_loci", "total_reads", "mean_reads_per_locus", "max_reads_per_locus"},
        )
        for row in reader:
            if row["motif"].upper() == motif.upper() and row["tool"] in tools:
                rows[row["tool"]] = row
    return rows


def print_aggregate_fallback(
    length_rows: dict[str, dict[str, str]],
    exact_rows: dict[str, dict[str, str]],
    tools: list[str],
    target_tool: str,
    target_coverage: int,
) -> None:
    print()
    print("Aggregate fallback result:")
    print(
        "Because the current locus/per-read detail does not contain the target outlier, "
        "the exact coordinates cannot be recovered from aggregate TSVs alone."
    )
    print(
        "For PacBio, the 25468 value is directly removed. For ONT, the script removes "
        "the max AATGG ONT coverage as the paired-locus estimate."
    )
    print()
    print("tool\tn_loci_before\ttotal_reads_before\tremoved_reads\tn_loci_after\ttotal_reads_after\tmean_after_excluding_outlier")

    for tool in tools:
        if tool not in length_rows:
            print(f"{tool}\tNA\tNA\tNA\tNA\tNA\tNA")
            continue
        n_loci = int(length_rows[tool]["n_loci"])
        total_reads = int(length_rows[tool]["total_reads"])
        if tool == target_tool:
            removed_reads = target_coverage
        elif tool in exact_rows:
            removed_reads = int(exact_rows[tool]["max_reads_per_locus"])
        else:
            removed_reads = int(length_rows[tool]["max_reads_per_locus"])
        n_after = n_loci - 1
        total_after = total_reads - removed_reads
        mean_after = total_after / n_after if n_after > 0 else float("nan")
        print(f"{tool}\t{n_loci}\t{total_reads}\t{removed_reads}\t{n_after}\t{total_after}\t{mean_after:.6f}")


def main() -> None:
    args = parse_args()
    if args.exclude_coverage_above is not None and args.exclude_coverage_above < 0:
        raise SystemExit("--exclude-coverage-above must be >= 0")

    tools = list(args.tools)
    tool_set = set(tools)

    summary_rows = read_summary_rows(args.summary_tsv, tool_set, args.motif_length)
    per_read_rows: list[LocusRow] = []
    per_read_path = Path(args.per_read_tsv)
    if args.per_read_tsv and per_read_path.exists():
        per_read_rows = read_per_read_counts(args.per_read_tsv, tool_set, args.motif_length)

    if args.exclude_coverage_above is not None:
        print("Filter mode:")
        print(f"motif_length={args.motif_length}\texclude_coverage_above={args.exclude_coverage_above}")
        print()
        print_filtered_mean(summary_rows, tools, args.exclude_coverage_above)
        return

    print("Search target:")
    print(f"motif={args.motif.upper()}\ttool={args.tool}\tcoverage={args.coverage}\tmotif_length={args.motif_length}")
    print()

    outlier = find_target(summary_rows, args.motif, args.tool, args.coverage)
    source = "summary TSV"
    if outlier is None and per_read_rows:
        outlier = find_target(per_read_rows, args.motif, args.tool, args.coverage)
        source = "per-read TSV"

    if outlier is not None:
        print(f"Located outlier from {source}:")
        print("motif\ttool\tcoord\tn_reads")
        print(f"{outlier.motif}\t{outlier.tool}\t{outlier.coord}\t{outlier.n_reads}")
        rows_for_mean = summary_rows if source == "summary TSV" else per_read_rows
        print_exact_adjusted_mean(rows_for_mean, outlier, tools)
        return

    print("Target outlier was not found in the current detail files.")
    print_available_motif_rows("summary TSV check", summary_rows, args.motif)
    if per_read_rows:
        print_available_motif_rows("per-read TSV check", per_read_rows, args.motif)
    else:
        print(f"per-read TSV check: skipped because {args.per_read_tsv} was not found")

    length_rows = read_length_5_aggregates(args.coverage_by_motif_length_tsv, tool_set)
    exact_rows = read_exact_motif_aggregates(args.coverage_by_exact_motif_tsv, args.motif, tool_set)
    print_aggregate_fallback(length_rows, exact_rows, tools, args.tool, args.coverage)


if __name__ == "__main__":
    main()
