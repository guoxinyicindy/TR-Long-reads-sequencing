#!/usr/bin/env python3
"""
Calculate per-locus standard error of read-reference Levenshtein distances.

Distances are read from read_ref_lev.per_read.tsv. If a summary TSV is provided,
all locus/tool rows from the summary are retained in the per-locus output, even
when there are zero reads or only one read. Motif bins are assigned from BED
column 4 by default, or from locus_id when --bed is not provided.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


BIN_ORDER = ["1bp", "2-10bp", "11-30bp", "31-50bp", "51-100bp", ">100bp"]
TOOL_ORDER = ["ONT", "PacBio", "realigned"]
LocusKey = tuple[str, str, str]
LocusToolKey = tuple[str, str, str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-locus standard error of Levenshtein distance and "
            "summarize it by motif length bins."
        )
    )
    parser.add_argument(
        "--per-read-tsv",
        required=True,
        help="Input read_ref_lev.per_read.tsv, filtered or unfiltered.",
    )
    parser.add_argument(
        "--summary-tsv",
        help=(
            "Optional read_ref_lev.summary.tsv. If provided, locus/tool rows "
            "with zero reads are kept in the per-locus output."
        ),
    )
    parser.add_argument(
        "--bed",
        help=(
            "Optional BED file used to assign motif bins by chrom/start/end. "
            "BED column 4 is used as motif by default."
        ),
    )
    parser.add_argument(
        "--motif-column",
        type=int,
        default=4,
        help="1-based motif column in --bed. Default: 4.",
    )
    parser.add_argument(
        "--out-locus-tsv",
        required=True,
        help="Output per-locus standard error TSV.",
    )
    parser.add_argument(
        "--out-bin-tsv",
        required=True,
        help="Output motif-bin standard error summary TSV.",
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


def numeric_or_none(value: str) -> float | None:
    if value in {"", "NA", "nan", "NaN"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def fmt(value: float | int | str | None) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return value
    return f"{value:.6f}"


def mean_or_none(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def median_or_none(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def sample_sd_or_none(values: list[float]) -> float | None:
    return statistics.stdev(values) if len(values) >= 2 else None


def se_or_none(values: list[float]) -> float | None:
    sd = sample_sd_or_none(values)
    if sd is None:
        return None
    return sd / math.sqrt(len(values))


def min_or_none(values: list[float]) -> float | None:
    return min(values) if values else None


def max_or_none(values: list[float]) -> float | None:
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


def load_summary_loci(
    summary_tsv: str | None,
    motif_by_locus: dict[LocusKey, str],
) -> dict[LocusToolKey, dict[str, str]]:
    if summary_tsv is None:
        return {}

    loci: dict[LocusToolKey, dict[str, str]] = {}
    with open(summary_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"scope", "locus_id", "chrom", "start", "end", "tool"}
        require_columns(summary_tsv, reader.fieldnames, required)
        for row in reader:
            if row["scope"] != "locus":
                continue
            motif = row_motif(row, motif_by_locus)
            if motif_by_locus and motif is None:
                continue
            key = row_locus_tool_key(row)
            loci[key] = {
                "locus_id": row["locus_id"],
                "chrom": key[0],
                "start": key[1],
                "end": key[2],
                "tool": key[3],
                "motif": motif or row["locus_id"],
            }
    return loci


def load_per_read_distances(
    per_read_tsv: str,
    motif_by_locus: dict[LocusKey, str],
) -> tuple[dict[LocusToolKey, list[float]], dict[LocusToolKey, dict[str, str]], int]:
    distances_by_locus: dict[LocusToolKey, list[float]] = defaultdict(list)
    locus_meta: dict[LocusToolKey, dict[str, str]] = {}
    skipped_missing_motif = 0

    with open(per_read_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"locus_id", "chrom", "start", "end", "tool", "distance"}
        require_columns(per_read_tsv, reader.fieldnames, required)
        for row in reader:
            motif = row_motif(row, motif_by_locus)
            if motif_by_locus and motif is None:
                skipped_missing_motif += 1
                continue
            distance = numeric_or_none(row["distance"])
            if distance is None:
                continue
            key = row_locus_tool_key(row)
            distances_by_locus[key].append(distance)
            locus_meta.setdefault(
                key,
                {
                    "locus_id": row["locus_id"],
                    "chrom": key[0],
                    "start": key[1],
                    "end": key[2],
                    "tool": key[3],
                    "motif": motif or row["locus_id"],
                },
            )
    return distances_by_locus, locus_meta, skipped_missing_motif


def tool_sort_key(tool: str) -> tuple[int, str]:
    try:
        return TOOL_ORDER.index(tool), tool
    except ValueError:
        return len(TOOL_ORDER), tool


def locus_sort_key(meta: dict[str, str]) -> tuple[str, int, int, tuple[int, str]]:
    return meta["chrom"], int(meta["start"]), int(meta["end"]), tool_sort_key(meta["tool"])


def build_locus_rows(
    summary_loci: dict[LocusToolKey, dict[str, str]],
    per_read_loci: dict[LocusToolKey, dict[str, str]],
    distances_by_locus: dict[LocusToolKey, list[float]],
) -> list[dict[str, str | int | float | None]]:
    all_keys = set(summary_loci).union(per_read_loci)
    rows: list[dict[str, str | int | float | None]] = []

    for key in all_keys:
        meta = summary_loci.get(key) or per_read_loci[key]
        distances = distances_by_locus.get(key, [])
        rows.append(
            {
                "motif_bin": motif_len_bin(meta["motif"]),
                "motif": meta["motif"],
                "locus_id": meta["locus_id"],
                "chrom": meta["chrom"],
                "start": meta["start"],
                "end": meta["end"],
                "tool": meta["tool"],
                "n_reads": len(distances),
                "mean_distance": mean_or_none(distances),
                "sd_distance": sample_sd_or_none(distances),
                "se_distance": se_or_none(distances),
                "median_distance": median_or_none(distances),
                "min_distance": min_or_none(distances),
                "max_distance": max_or_none(distances),
            }
        )
    return sorted(rows, key=lambda row: locus_sort_key(row))  # type: ignore[arg-type]


def write_locus_rows(rows: list[dict[str, str | int | float | None]], out_tsv: str) -> None:
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "motif_bin",
        "motif",
        "locus_id",
        "chrom",
        "start",
        "end",
        "tool",
        "n_reads",
        "mean_distance",
        "sd_distance",
        "se_distance",
        "median_distance",
        "min_distance",
        "max_distance",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([fmt(row[column]) for column in header])


def write_bin_summary(
    rows: list[dict[str, str | int | float | None]],
    distances_by_locus: dict[LocusToolKey, list[float]],
    out_tsv: str,
) -> None:
    grouped_rows: dict[tuple[str, str], list[dict[str, str | int | float | None]]] = defaultdict(list)
    grouped_read_distances: dict[tuple[str, str], list[float]] = defaultdict(list)

    for row in rows:
        key = (str(row["motif_bin"]), str(row["tool"]))
        grouped_rows[key].append(row)
        locus_key = (
            str(row["chrom"]),
            str(row["start"]),
            str(row["end"]),
            str(row["tool"]),
        )
        grouped_read_distances[key].extend(distances_by_locus.get(locus_key, []))

    tools = sorted({str(row["tool"]) for row in rows}, key=tool_sort_key)
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "motif_bin",
        "tool",
        "n_loci",
        "n_loci_with_reads",
        "n_loci_with_se",
        "n_reads",
        "mean_locus_se_distance",
        "median_locus_se_distance",
        "read_weighted_mean_locus_se_distance",
        "mean_locus_sd_distance",
        "median_locus_sd_distance",
        "read_weighted_mean_locus_sd_distance",
        "mean_locus_mean_distance",
        "median_locus_mean_distance",
        "pooled_read_sd_distance",
        "pooled_read_se_distance",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for motif_bin in BIN_ORDER:
            for tool in tools:
                key = (motif_bin, tool)
                group = grouped_rows.get(key, [])
                pooled_distances = grouped_read_distances.get(key, [])
                se_values = [row["se_distance"] for row in group if isinstance(row["se_distance"], float)]
                sd_values = [row["sd_distance"] for row in group if isinstance(row["sd_distance"], float)]
                mean_values = [row["mean_distance"] for row in group if isinstance(row["mean_distance"], float)]
                n_reads_total = sum(int(row["n_reads"]) for row in group)
                weighted_se_sum = sum(
                    float(row["se_distance"]) * int(row["n_reads"])
                    for row in group
                    if isinstance(row["se_distance"], float) and int(row["n_reads"]) > 0
                )
                weighted_se_n = sum(
                    int(row["n_reads"])
                    for row in group
                    if isinstance(row["se_distance"], float) and int(row["n_reads"]) > 0
                )
                weighted_sd_sum = sum(
                    float(row["sd_distance"]) * int(row["n_reads"])
                    for row in group
                    if isinstance(row["sd_distance"], float) and int(row["n_reads"]) > 0
                )
                weighted_sd_n = sum(
                    int(row["n_reads"])
                    for row in group
                    if isinstance(row["sd_distance"], float) and int(row["n_reads"]) > 0
                )
                writer.writerow(
                    [
                        motif_bin,
                        tool,
                        len(group),
                        sum(1 for row in group if int(row["n_reads"]) > 0),
                        len(se_values),
                        n_reads_total,
                        fmt(mean_or_none(se_values)),
                        fmt(median_or_none(se_values)),
                        fmt(weighted_se_sum / weighted_se_n if weighted_se_n else None),
                        fmt(mean_or_none(sd_values)),
                        fmt(median_or_none(sd_values)),
                        fmt(weighted_sd_sum / weighted_sd_n if weighted_sd_n else None),
                        fmt(mean_or_none(mean_values)),
                        fmt(median_or_none(mean_values)),
                        fmt(sample_sd_or_none(pooled_distances)),
                        fmt(se_or_none(pooled_distances)),
                    ]
                )


def main() -> None:
    args = parse_args()
    motif_by_locus = load_motif_map(args.bed, args.motif_column)
    summary_loci = load_summary_loci(args.summary_tsv, motif_by_locus)
    distances_by_locus, per_read_loci, skipped_missing_motif = load_per_read_distances(
        args.per_read_tsv,
        motif_by_locus,
    )
    locus_rows = build_locus_rows(summary_loci, per_read_loci, distances_by_locus)
    write_locus_rows(locus_rows, args.out_locus_tsv)
    write_bin_summary(locus_rows, distances_by_locus, args.out_bin_tsv)

    if skipped_missing_motif:
        print(
            f"warning: skipped {skipped_missing_motif} per-read rows without BED motif mapping",
            file=sys.stderr,
        )
    print(
        (
            "[read-ref-lev-se] "
            f"locus_rows={len(locus_rows)} "
            f"out_locus={args.out_locus_tsv} "
            f"out_bin={args.out_bin_tsv}"
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
