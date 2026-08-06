#!/usr/bin/env python3
"""
Summarize per-locus read coverage by exact motif sequence.

This is useful for questions like: among 5 bp motifs, what is the average read
coverage for AAAAC, AAAAG, and each other motif sequence?
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


LocusKey = tuple[str, str, str]
TOOL_ORDER = ["ONT", "PacBio", "realigned"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize average read coverage by exact motif sequence."
    )
    parser.add_argument(
        "--summary-tsv",
        required=True,
        help="Input read_ref_lev.summary.tsv with scope=locus rows and n_reads.",
    )
    parser.add_argument(
        "--bed",
        help=(
            "Optional BED file used to assign motifs by chrom/start/end. "
            "BED column 4 is used by default. If omitted, locus_id is treated as the motif."
        ),
    )
    parser.add_argument(
        "--motif-column",
        type=int,
        default=4,
        help="1-based motif column in --bed. Default: 4.",
    )
    parser.add_argument(
        "--motif-length",
        type=int,
        default=5,
        help="Exact motif length to keep. Default: 5.",
    )
    parser.add_argument(
        "--keep-case",
        action="store_true",
        help="Do not normalize motifs to uppercase before grouping.",
    )
    parser.add_argument(
        "--out-tsv",
        required=True,
        help="Output TSV with one row per motif/tool.",
    )
    parser.add_argument(
        "--tools",
        nargs="+",
        default=["ONT", "PacBio"],
        help="Tools to include. Default: ONT PacBio.",
    )
    parser.add_argument(
        "--print-outlier-adjusted-mean",
        action="store_true",
        help=(
            "Print to stdout the max-coverage locus for the selected motif length "
            "and the mean coverage after excluding that locus per tool."
        ),
    )
    parser.add_argument(
        "--outlier-motif",
        help="Optional exact motif sequence to inspect. Default: inspect all motifs of --motif-length.",
    )
    parser.add_argument(
        "--exclude-coverage-above",
        type=int,
        help=(
            "Exclude an entire locus from the summary if any selected tool has "
            "n_reads greater than this value. Excluded loci are printed to stdout."
        ),
    )
    return parser.parse_args()


def require_columns(path: str, fieldnames: list[str] | None, required: set[str]) -> None:
    missing = required.difference(fieldnames or [])
    if missing:
        raise SystemExit(f"{path} is missing columns: {', '.join(sorted(missing))}")


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


def normalize_motif(motif: str, keep_case: bool) -> str:
    motif = motif.strip()
    return motif if keep_case else motif.upper()


def load_motif_map(
    bed_path: str | None,
    motif_column: int,
    keep_case: bool,
) -> dict[LocusKey, str]:
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
            key = row[0], str(int(row[1])), str(int(row[2]))
            motif_by_locus[key] = normalize_motif(row[motif_index], keep_case)
    return motif_by_locus


def row_locus_key(row: dict[str, str]) -> LocusKey:
    return row["chrom"], str(int(row["start"])), str(int(row["end"]))


def row_motif(row: dict[str, str], motif_by_locus: dict[LocusKey, str], keep_case: bool) -> str | None:
    if motif_by_locus:
        return motif_by_locus.get(row_locus_key(row))
    return normalize_motif(row.get("locus_id", ""), keep_case)


def tool_sort_key(tool: str) -> tuple[int, str]:
    try:
        return TOOL_ORDER.index(tool), tool
    except ValueError:
        return len(TOOL_ORDER), tool


def load_coverages(
    summary_tsv: str,
    motif_by_locus: dict[LocusKey, str],
    motif_length: int,
    keep_case: bool,
    tools: set[str],
) -> tuple[dict[tuple[str, str], list[int]], int, int]:
    coverages: dict[tuple[str, str], list[int]] = defaultdict(list)
    skipped_missing_motif = 0
    kept_locus_tool_rows = 0

    with open(summary_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"scope", "locus_id", "chrom", "start", "end", "tool", "n_reads"}
        require_columns(summary_tsv, reader.fieldnames, required)

        for row in reader:
            if row["scope"] != "locus":
                continue
            if row["tool"] not in tools:
                continue
            motif = row_motif(row, motif_by_locus, keep_case)
            if motif_by_locus and motif is None:
                skipped_missing_motif += 1
                continue
            if motif is None or len(motif) != motif_length:
                continue
            try:
                n_reads = int(row["n_reads"])
            except ValueError as exc:
                raise SystemExit(f"{summary_tsv} has non-integer n_reads: {row['n_reads']}") from exc
            if n_reads < 0:
                raise SystemExit(f"{summary_tsv} has negative n_reads: {n_reads}")
            coverages[(motif, row["tool"])].append(n_reads)
            kept_locus_tool_rows += 1

    return coverages, kept_locus_tool_rows, skipped_missing_motif


def load_locus_rows(
    summary_tsv: str,
    motif_by_locus: dict[LocusKey, str],
    motif_length: int,
    keep_case: bool,
    tools: set[str],
    outlier_motif: str | None,
) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    target_motif = normalize_motif(outlier_motif, keep_case) if outlier_motif else None

    with open(summary_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"scope", "locus_id", "chrom", "start", "end", "tool", "n_reads"}
        require_columns(summary_tsv, reader.fieldnames, required)

        for row in reader:
            if row["scope"] != "locus" or row["tool"] not in tools:
                continue
            motif = row_motif(row, motif_by_locus, keep_case)
            if motif is None or len(motif) != motif_length:
                continue
            if target_motif is not None and motif != target_motif:
                continue
            n_reads = int(row["n_reads"])
            rows.append(
                {
                    "motif": motif,
                    "tool": row["tool"],
                    "chrom": row["chrom"],
                    "start": str(int(row["start"])),
                    "end": str(int(row["end"])),
                    "locus_id": row["locus_id"],
                    "n_reads": n_reads,
                }
            )
    return rows


def build_coverages_from_locus_rows(
    rows: list[dict[str, str | int]],
    exclude_coverage_above: int | None,
) -> tuple[dict[tuple[str, str], list[int]], int, list[tuple[tuple[str, str, str], list[dict[str, str | int]]]]]:
    coverages: dict[tuple[str, str], list[int]] = defaultdict(list)

    rows_by_locus: dict[tuple[str, str, str], list[dict[str, str | int]]] = defaultdict(list)
    for row in rows:
        key = (str(row["chrom"]), str(row["start"]), str(row["end"]))
        rows_by_locus[key].append(row)

    excluded_loci: list[tuple[tuple[str, str, str], list[dict[str, str | int]]]] = []
    excluded_keys: set[tuple[str, str, str]] = set()
    if exclude_coverage_above is not None:
        for key, locus_rows in rows_by_locus.items():
            if any(int(row["n_reads"]) > exclude_coverage_above for row in locus_rows):
                excluded_keys.add(key)
                excluded_loci.append((key, locus_rows))

    kept_locus_tool_rows = 0
    for row in rows:
        key = (str(row["chrom"]), str(row["start"]), str(row["end"]))
        if key in excluded_keys:
            continue
        coverages[(str(row["motif"]), str(row["tool"]))].append(int(row["n_reads"]))
        kept_locus_tool_rows += 1

    excluded_loci.sort(key=lambda item: (item[0][0], int(item[0][1]), int(item[0][2])))
    return coverages, kept_locus_tool_rows, excluded_loci


def print_excluded_loci(
    excluded_loci: list[tuple[tuple[str, str, str], list[dict[str, str | int]]]],
    threshold: int,
) -> None:
    print(f"Excluded loci with any selected-tool n_reads > {threshold}:")
    if not excluded_loci:
        print("None")
        return

    print(
        "\t".join(
            [
                "chrom",
                "start",
                "end",
                "motif",
                "excluded_tools",
                "all_tool_n_reads",
            ]
        )
    )
    for (chrom, start, end), rows in excluded_loci:
        motif = str(rows[0]["motif"])
        triggering = [
            f"{row['tool']}={row['n_reads']}"
            for row in sorted(rows, key=lambda row: tool_sort_key(str(row["tool"])))
            if int(row["n_reads"]) > threshold
        ]
        all_counts = [
            f"{row['tool']}={row['n_reads']}"
            for row in sorted(rows, key=lambda row: tool_sort_key(str(row["tool"])))
        ]
        print(
            "\t".join(
                [
                    chrom,
                    start,
                    end,
                    motif,
                    ",".join(triggering),
                    ",".join(all_counts),
                ]
            )
        )


def print_outlier_adjusted_mean(rows: list[dict[str, str | int]]) -> None:
    if not rows:
        print("No matching locus rows found.")
        return

    max_row = max(rows, key=lambda row: int(row["n_reads"]))
    outlier_key = (max_row["chrom"], max_row["start"], max_row["end"])

    print("Max-coverage locus:")
    print(
        "\t".join(
            [
                "motif",
                "tool",
                "chrom",
                "start",
                "end",
                "locus_id",
                "n_reads",
            ]
        )
    )
    print(
        "\t".join(
            [
                str(max_row["motif"]),
                str(max_row["tool"]),
                str(max_row["chrom"]),
                str(max_row["start"]),
                str(max_row["end"]),
                str(max_row["locus_id"]),
                str(max_row["n_reads"]),
            ]
        )
    )
    print()
    print("Mean read coverage after excluding that locus:")
    print("\t".join(["tool", "n_loci", "total_reads", "mean_reads_per_locus"]))

    for tool in sorted({str(row["tool"]) for row in rows}, key=tool_sort_key):
        values = [
            int(row["n_reads"])
            for row in rows
            if row["tool"] == tool and (row["chrom"], row["start"], row["end"]) != outlier_key
        ]
        mean_value = statistics.fmean(values) if values else None
        print(
            "\t".join(
                [
                    tool,
                    str(len(values)),
                    str(sum(values)),
                    fmt(mean_value),
                ]
            )
        )


def write_output(coverages: dict[tuple[str, str], list[int]], out_tsv: str, motif_length: int) -> None:
    motifs = sorted({motif for motif, _tool in coverages})
    tools = sorted({tool for _motif, tool in coverages}, key=tool_sort_key)
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "motif",
        "motif_length",
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
        for motif in motifs:
            for tool in tools:
                values = coverages.get((motif, tool), [])
                if not values:
                    continue
                covered_values = [value for value in values if value > 0]
                writer.writerow(
                    [
                        motif,
                        motif_length,
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
    if args.motif_length <= 0:
        raise SystemExit("--motif-length must be > 0")

    motif_by_locus = load_motif_map(args.bed, args.motif_column, args.keep_case)
    if args.exclude_coverage_above is not None and args.exclude_coverage_above < 0:
        raise SystemExit("--exclude-coverage-above must be >= 0")

    rows = load_locus_rows(
        args.summary_tsv,
        motif_by_locus,
        args.motif_length,
        args.keep_case,
        set(args.tools),
        None,
    )
    coverages, kept_locus_tool_rows, excluded_loci = build_coverages_from_locus_rows(
        rows,
        args.exclude_coverage_above,
    )
    skipped_missing_motif = 0
    write_output(coverages, args.out_tsv, args.motif_length)

    if args.exclude_coverage_above is not None:
        print_excluded_loci(excluded_loci, args.exclude_coverage_above)
        print()

    if args.print_outlier_adjusted_mean:
        outlier_rows = load_locus_rows(
            args.summary_tsv,
            motif_by_locus,
            args.motif_length,
            args.keep_case,
            set(args.tools),
            args.outlier_motif,
        )
        if args.exclude_coverage_above is not None:
            excluded_keys = {key for key, _locus_rows in excluded_loci}
            outlier_rows = [
                row
                for row in outlier_rows
                if (str(row["chrom"]), str(row["start"]), str(row["end"])) not in excluded_keys
            ]
        print_outlier_adjusted_mean(outlier_rows)

    if skipped_missing_motif:
        print(
            f"warning: skipped {skipped_missing_motif} summary rows without BED motif mapping",
            file=sys.stderr,
        )
    print(
        (
            "[coverage-by-exact-motif] "
            f"motif_length={args.motif_length} "
            f"motifs={len({motif for motif, _tool in coverages})} "
            f"locus_tool_rows={kept_locus_tool_rows} "
            f"excluded_loci={len(excluded_loci)} "
            f"out={args.out_tsv}"
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
