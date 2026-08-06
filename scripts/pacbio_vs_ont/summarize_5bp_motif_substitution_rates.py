#!/usr/bin/env python3
"""
Summarize mean substitution rate per bp for each 5bp motif.

Input is the per-read TSV from decompose_read_ref_edits_ont_pacbio.py. By
default, locus_id is treated as the motif. If locus_id is not the motif, pass
--bed so BED column 4 can be used as the motif by chrom/start/end.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


TOOL_ORDER = ["ONT", "PacBio", "realigned"]
LocusKey = tuple[str, str, str]
LocusToolKey = tuple[str, str, str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize sub_per_ref_bp for each exact 5bp motif."
    )
    parser.add_argument(
        "--per-read-tsv",
        required=True,
        help="Per-read TSV from decompose_read_ref_edits_ont_pacbio.py.",
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
        "--motif-length",
        type=int,
        default=5,
        help="Exact motif length to summarize. Default: 5.",
    )
    parser.add_argument(
        "--out-tsv",
        required=True,
        help="Output TSV grouped by motif and tool.",
    )
    return parser.parse_args()


def numeric_or_none(value: str) -> float | None:
    if value in {"", "NA", "nan", "NaN"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def fmt(value: float | int | None) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{value:.8f}"


def mean_or_none(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def median_or_none(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


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


def tool_sort_key(tool: str) -> tuple[int, str]:
    try:
        return TOOL_ORDER.index(tool), tool
    except ValueError:
        return len(TOOL_ORDER), tool


def load_substitution_rates(
    per_read_tsv: str,
    motif_by_locus: dict[LocusKey, str],
    motif_length: int,
) -> tuple[
    dict[tuple[str, str], list[float]],
    dict[tuple[str, str], dict[LocusToolKey, list[float]]],
    int,
    int,
]:
    per_read_rates: dict[tuple[str, str], list[float]] = defaultdict(list)
    rates_by_locus: dict[tuple[str, str], dict[LocusToolKey, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    skipped_missing_motif = 0
    skipped_wrong_length = 0

    with open(per_read_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"locus_id", "chrom", "start", "end", "tool", "sub_per_ref_bp"}
        require_columns(per_read_tsv, reader.fieldnames, required)
        for row in reader:
            motif = row_motif(row, motif_by_locus)
            if motif_by_locus and motif is None:
                skipped_missing_motif += 1
                continue
            motif_value = (motif or row["locus_id"]).upper()
            if len(motif_value) != motif_length:
                skipped_wrong_length += 1
                continue
            value = numeric_or_none(row["sub_per_ref_bp"])
            if value is None:
                continue
            key = (motif_value, row["tool"])
            per_read_rates[key].append(value)
            rates_by_locus[key][row_locus_tool_key(row)].append(value)
    return per_read_rates, rates_by_locus, skipped_missing_motif, skipped_wrong_length


def write_output(
    per_read_rates: dict[tuple[str, str], list[float]],
    rates_by_locus: dict[tuple[str, str], dict[LocusToolKey, list[float]]],
    out_tsv: str,
) -> None:
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "motif",
        "motif_length",
        "tool",
        "n_loci",
        "n_reads",
        "mean_sub_per_ref_bp",
        "median_sub_per_ref_bp",
        "mean_locus_mean_sub_per_ref_bp",
        "median_locus_mean_sub_per_ref_bp",
    ]
    keys = sorted(per_read_rates, key=lambda item: (item[0], tool_sort_key(item[1])))
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for motif, tool in keys:
            read_values = per_read_rates[(motif, tool)]
            locus_means = [
                statistics.fmean(values)
                for values in rates_by_locus[(motif, tool)].values()
                if values
            ]
            writer.writerow(
                [
                    motif,
                    len(motif),
                    tool,
                    len(locus_means),
                    len(read_values),
                    fmt(mean_or_none(read_values)),
                    fmt(median_or_none(read_values)),
                    fmt(mean_or_none(locus_means)),
                    fmt(median_or_none(locus_means)),
                ]
            )


def main() -> None:
    args = parse_args()
    if args.motif_length < 1:
        raise SystemExit("--motif-length must be >= 1")
    motif_by_locus = load_motif_map(args.bed, args.motif_column)
    per_read_rates, rates_by_locus, skipped_missing_motif, skipped_wrong_length = load_substitution_rates(
        args.per_read_tsv,
        motif_by_locus,
        args.motif_length,
    )
    write_output(per_read_rates, rates_by_locus, args.out_tsv)
    if skipped_missing_motif:
        print(
            f"warning: skipped {skipped_missing_motif} per-read rows without BED motif mapping",
            file=sys.stderr,
        )
    print(
        (
            "[5bp-motif-substitution-rates] "
            f"motifs={len({key[0] for key in per_read_rates})} "
            f"skipped_wrong_length={skipped_wrong_length} "
            f"out={args.out_tsv}"
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
