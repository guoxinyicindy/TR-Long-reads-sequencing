#!/usr/bin/env python3
"""
Summarize insertion/deletion/substitution rates by motif bins after removing
high-coverage loci.

This combines the broad-bin and exact 1-10bp summaries from
summarize_edit_rates_by_motif_bins.py with the coverage filter used by
summarize_5bp_sub_rate_after_coverage_filter.py.

For each locus/tool pair, n_reads is read from scope=locus rows in the summary
TSV. By default, per-read rows from locus/tool pairs with n_reads greater than
--max-reads-per-locus are excluded before averaging ins_per_ref_bp,
del_per_ref_bp, and sub_per_ref_bp.

By default, locus_id is treated as the motif. If locus_id is not the motif,
pass --bed so BED column 4 can be used as the motif by chrom/start/end.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


BIN_ORDER = ["1bp", "2-10bp", "11-30bp", "31-50bp", "51-100bp", ">100bp"]
TOOL_ORDER = ["ONT", "PacBio", "realigned"]
RATE_COLUMNS = ["ins_per_ref_bp", "del_per_ref_bp", "sub_per_ref_bp"]
RATE_PREFIXES = {
    "ins_per_ref_bp": "ins",
    "del_per_ref_bp": "del",
    "sub_per_ref_bp": "sub",
}
LocusKey = tuple[str, str, str]
LocusToolKey = tuple[str, str, str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize edit operation rates by motif bins and exact motif "
            "lengths after excluding high-coverage locus/tool pairs."
        )
    )
    parser.add_argument(
        "--summary-tsv",
        required=True,
        help="Input edit-op summary TSV with scope=locus rows and n_reads.",
    )
    parser.add_argument(
        "--per-read-tsv",
        required=True,
        help="Per-read TSV from decompose_read_ref_edits_ont_pacbio.py.",
    )
    parser.add_argument(
        "--max-reads-per-locus",
        type=int,
        required=True,
        help="Exclude a locus/tool pair if n_reads is greater than this value.",
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
        "--exclude-mode",
        choices=["locus-tool", "locus"],
        default="locus-tool",
        help=(
            "locus-tool excludes only the high-coverage technology at a locus; "
            "locus excludes the whole locus for every tool if any tool is above "
            "threshold. Default: locus-tool, matching "
            "summarize_5bp_sub_rate_after_coverage_filter.py."
        ),
    )
    parser.add_argument(
        "--out-bin-tsv",
        help=(
            "Output TSV grouped by broad motif length bins. Default is derived "
            "from --per-read-tsv."
        ),
    )
    parser.add_argument(
        "--out-length-1-10-tsv",
        help=(
            "Output TSV grouped by exact motif lengths from 1bp to 10bp. "
            "Default is derived from --per-read-tsv."
        ),
    )
    parser.add_argument(
        "--out-excluded-loci-tsv",
        help="Optional TSV listing excluded locus/tool pairs.",
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


def sorted_locus_tool_keys(keys: set[LocusToolKey]) -> list[LocusToolKey]:
    return sorted(keys, key=lambda item: (item[0], int(item[1]), int(item[2]), tool_sort_key(item[3])))


def default_output_path(per_read_tsv: str, suffix: str, max_reads_per_locus: int) -> str:
    path = Path(per_read_tsv)
    name = path.name
    if name.endswith(".per_read.tsv"):
        prefix = name[: -len(".per_read.tsv")]
    elif name.endswith(".tsv"):
        prefix = name[:-4]
    else:
        prefix = name
    return str(path.with_name(f"{prefix}.max_reads_{max_reads_per_locus}.filtered.{suffix}.tsv"))


def load_filter_sets(
    summary_tsv: str,
    motif_by_locus: dict[LocusKey, str],
    max_reads_per_locus: int,
    exclude_mode: str,
) -> tuple[set[LocusToolKey], set[LocusToolKey], dict[LocusToolKey, dict[str, str | int]], int]:
    all_locus_tools: set[LocusToolKey] = set()
    initially_excluded: set[LocusToolKey] = set()
    meta_by_key: dict[LocusToolKey, dict[str, str | int]] = {}
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
            motif_value = (motif or row["locus_id"]).upper()
            try:
                n_reads = int(row["n_reads"])
            except ValueError as exc:
                raise SystemExit(f"{summary_tsv} has non-integer n_reads: {row['n_reads']}") from exc
            if n_reads < 0:
                raise SystemExit(f"{summary_tsv} has negative n_reads: {n_reads}")

            key = row_locus_tool_key(row)
            all_locus_tools.add(key)
            meta_by_key[key] = {
                "motif": motif_value,
                "locus_id": row["locus_id"],
                "chrom": key[0],
                "start": key[1],
                "end": key[2],
                "tool": key[3],
                "n_reads": n_reads,
            }
            if n_reads > max_reads_per_locus:
                initially_excluded.add(key)

    if exclude_mode == "locus":
        excluded_loci = {key[:3] for key in initially_excluded}
        excluded = {key for key in all_locus_tools if key[:3] in excluded_loci}
    else:
        excluded = initially_excluded

    included = all_locus_tools.difference(excluded)
    return included, excluded, meta_by_key, skipped_missing_motif


def load_filtered_rates(
    per_read_tsv: str,
    motif_by_locus: dict[LocusKey, str],
    included_loci: set[LocusToolKey],
    excluded_loci: set[LocusToolKey],
) -> tuple[
    dict[tuple[str, str], dict[str, list[float]]],
    dict[tuple[int, str], dict[str, list[float]]],
    dict[tuple[str, str], set[LocusToolKey]],
    dict[tuple[int, str], set[LocusToolKey]],
    int,
    int,
    int,
]:
    rates_by_bin_tool: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    rates_by_length_tool: dict[tuple[int, str], dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    loci_by_bin_tool: dict[tuple[str, str], set[LocusToolKey]] = defaultdict(set)
    loci_by_length_tool: dict[tuple[int, str], set[LocusToolKey]] = defaultdict(set)
    skipped_missing_motif = 0
    skipped_excluded_reads = 0
    skipped_not_in_summary = 0

    with open(per_read_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"locus_id", "chrom", "start", "end", "tool", *RATE_COLUMNS}
        require_columns(per_read_tsv, reader.fieldnames, required)
        for row in reader:
            motif = row_motif(row, motif_by_locus)
            if motif_by_locus and motif is None:
                skipped_missing_motif += 1
                continue
            key = row_locus_tool_key(row)
            if key in excluded_loci:
                skipped_excluded_reads += 1
                continue
            if key not in included_loci:
                skipped_not_in_summary += 1
                continue

            motif_value = (motif or row["locus_id"]).upper()
            motif_length = len(motif_value)
            tool = row["tool"]
            bin_key = (motif_len_bin(motif_value), tool)
            length_key = (motif_length, tool)
            loci_by_bin_tool[bin_key].add(key)
            if 1 <= motif_length <= 10:
                loci_by_length_tool[length_key].add(key)
            for column in RATE_COLUMNS:
                value = numeric_or_none(row[column])
                if value is None:
                    continue
                rates_by_bin_tool[bin_key][column].append(value)
                if 1 <= motif_length <= 10:
                    rates_by_length_tool[length_key][column].append(value)

    return (
        rates_by_bin_tool,
        rates_by_length_tool,
        loci_by_bin_tool,
        loci_by_length_tool,
        skipped_missing_motif,
        skipped_excluded_reads,
        skipped_not_in_summary,
    )


def metric_values(rates: dict[str, list[float]], metric: str) -> list[float]:
    return rates.get(metric, [])


def count_loci_by_category(
    meta_by_key: dict[LocusToolKey, dict[str, str | int]],
    keys: set[LocusToolKey],
    category_type: str,
) -> dict[tuple[str | int, str], int]:
    counts: dict[tuple[str | int, str], int] = defaultdict(int)
    for key in keys:
        meta = meta_by_key.get(key)
        if meta is None:
            continue
        motif = str(meta["motif"])
        if category_type == "bin":
            category: str | int = motif_len_bin(motif)
        else:
            motif_length = len(motif)
            if not 1 <= motif_length <= 10:
                continue
            category = motif_length
        counts[(category, key[3])] += 1
    return counts


def tools_from_loci_and_rates(
    included_loci: set[LocusToolKey],
    excluded_loci: set[LocusToolKey],
    *rate_maps: dict[tuple[str | int, str], dict[str, list[float]]],
) -> list[str]:
    tools = {key[3] for key in included_loci.union(excluded_loci)}
    for rate_map in rate_maps:
        tools.update(tool for _category, tool in rate_map)
    return sorted(tools, key=tool_sort_key)


def summary_values(rates: dict[str, list[float]], column: str) -> tuple[list[float], str, str]:
    values = metric_values(rates, column)
    prefix = RATE_PREFIXES[column]
    return values, f"mean_{prefix}_per_ref_bp", f"median_{prefix}_per_ref_bp"


def write_bin_output(
    rates_by_bin_tool: dict[tuple[str, str], dict[str, list[float]]],
    loci_by_bin_tool: dict[tuple[str, str], set[LocusToolKey]],
    included_loci: set[LocusToolKey],
    excluded_loci: set[LocusToolKey],
    meta_by_key: dict[LocusToolKey, dict[str, str | int]],
    max_reads_per_locus: int,
    exclude_mode: str,
    out_tsv: str,
) -> None:
    tools = tools_from_loci_and_rates(included_loci, excluded_loci, rates_by_bin_tool)
    before_counts = count_loci_by_category(meta_by_key, included_loci.union(excluded_loci), "bin")
    excluded_counts = count_loci_by_category(meta_by_key, excluded_loci, "bin")
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "motif_bin",
        "tool",
        "max_reads_per_locus",
        "exclude_mode",
        "n_loci_before_filter",
        "n_loci_after_filter",
        "n_loci_excluded",
        "n_reads_after_filter",
        "mean_ins_per_ref_bp",
        "median_ins_per_ref_bp",
        "mean_del_per_ref_bp",
        "median_del_per_ref_bp",
        "mean_sub_per_ref_bp",
        "median_sub_per_ref_bp",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for motif_bin in BIN_ORDER:
            for tool in tools:
                key = (motif_bin, tool)
                rates = rates_by_bin_tool.get(key, {})
                ins_values, _, _ = summary_values(rates, "ins_per_ref_bp")
                del_values, _, _ = summary_values(rates, "del_per_ref_bp")
                sub_values, _, _ = summary_values(rates, "sub_per_ref_bp")
                writer.writerow(
                    [
                        motif_bin,
                        tool,
                        max_reads_per_locus,
                        exclude_mode,
                        before_counts.get(key, 0),
                        len(loci_by_bin_tool.get(key, set())),
                        excluded_counts.get(key, 0),
                        len(ins_values),
                        fmt(mean_or_none(ins_values)),
                        fmt(median_or_none(ins_values)),
                        fmt(mean_or_none(del_values)),
                        fmt(median_or_none(del_values)),
                        fmt(mean_or_none(sub_values)),
                        fmt(median_or_none(sub_values)),
                    ]
                )


def write_length_output(
    rates_by_length_tool: dict[tuple[int, str], dict[str, list[float]]],
    loci_by_length_tool: dict[tuple[int, str], set[LocusToolKey]],
    included_loci: set[LocusToolKey],
    excluded_loci: set[LocusToolKey],
    meta_by_key: dict[LocusToolKey, dict[str, str | int]],
    max_reads_per_locus: int,
    exclude_mode: str,
    out_tsv: str,
) -> None:
    tools = tools_from_loci_and_rates(included_loci, excluded_loci, rates_by_length_tool)
    before_counts = count_loci_by_category(meta_by_key, included_loci.union(excluded_loci), "length")
    excluded_counts = count_loci_by_category(meta_by_key, excluded_loci, "length")
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "motif_length",
        "motif_length_label",
        "tool",
        "max_reads_per_locus",
        "exclude_mode",
        "n_loci_before_filter",
        "n_loci_after_filter",
        "n_loci_excluded",
        "n_reads_after_filter",
        "mean_ins_per_ref_bp",
        "median_ins_per_ref_bp",
        "mean_del_per_ref_bp",
        "median_del_per_ref_bp",
        "mean_sub_per_ref_bp",
        "median_sub_per_ref_bp",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for motif_length in range(1, 11):
            for tool in tools:
                key = (motif_length, tool)
                rates = rates_by_length_tool.get(key, {})
                ins_values, _, _ = summary_values(rates, "ins_per_ref_bp")
                del_values, _, _ = summary_values(rates, "del_per_ref_bp")
                sub_values, _, _ = summary_values(rates, "sub_per_ref_bp")
                writer.writerow(
                    [
                        motif_length,
                        f"{motif_length}bp",
                        tool,
                        max_reads_per_locus,
                        exclude_mode,
                        before_counts.get(key, 0),
                        len(loci_by_length_tool.get(key, set())),
                        excluded_counts.get(key, 0),
                        len(ins_values),
                        fmt(mean_or_none(ins_values)),
                        fmt(median_or_none(ins_values)),
                        fmt(mean_or_none(del_values)),
                        fmt(median_or_none(del_values)),
                        fmt(mean_or_none(sub_values)),
                        fmt(median_or_none(sub_values)),
                    ]
                )


def write_excluded_loci(
    path: str | None,
    excluded_loci: set[LocusToolKey],
    meta_by_key: dict[LocusToolKey, dict[str, str | int]],
) -> None:
    if path is None:
        return
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = ["motif", "motif_length", "locus_id", "chrom", "start", "end", "tool", "n_reads"]
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for key in sorted_locus_tool_keys(excluded_loci):
            meta = meta_by_key[key]
            motif = str(meta["motif"])
            writer.writerow(
                [
                    motif,
                    len(motif),
                    meta["locus_id"],
                    meta["chrom"],
                    meta["start"],
                    meta["end"],
                    meta["tool"],
                    meta["n_reads"],
                ]
            )


def main() -> None:
    args = parse_args()
    if args.max_reads_per_locus < 0:
        raise SystemExit("--max-reads-per-locus must be >= 0")

    out_bin_tsv = args.out_bin_tsv or default_output_path(
        args.per_read_tsv,
        "rates_by_motif_bins",
        args.max_reads_per_locus,
    )
    out_length_tsv = args.out_length_1_10_tsv or default_output_path(
        args.per_read_tsv,
        "rates_by_motif_length_1_10bp",
        args.max_reads_per_locus,
    )

    motif_by_locus = load_motif_map(args.bed, args.motif_column)
    included_loci, excluded_loci, meta_by_key, skipped_summary_missing_motif = load_filter_sets(
        args.summary_tsv,
        motif_by_locus,
        args.max_reads_per_locus,
        args.exclude_mode,
    )
    (
        rates_by_bin_tool,
        rates_by_length_tool,
        loci_by_bin_tool,
        loci_by_length_tool,
        skipped_per_read_missing_motif,
        skipped_excluded_reads,
        skipped_not_in_summary,
    ) = load_filtered_rates(
        args.per_read_tsv,
        motif_by_locus,
        included_loci,
        excluded_loci,
    )

    write_bin_output(
        rates_by_bin_tool,
        loci_by_bin_tool,
        included_loci,
        excluded_loci,
        meta_by_key,
        args.max_reads_per_locus,
        args.exclude_mode,
        out_bin_tsv,
    )
    write_length_output(
        rates_by_length_tool,
        loci_by_length_tool,
        included_loci,
        excluded_loci,
        meta_by_key,
        args.max_reads_per_locus,
        args.exclude_mode,
        out_length_tsv,
    )
    write_excluded_loci(args.out_excluded_loci_tsv, excluded_loci, meta_by_key)

    if skipped_summary_missing_motif or skipped_per_read_missing_motif:
        print(
            (
                "warning: skipped rows without BED motif mapping "
                f"summary={skipped_summary_missing_motif} "
                f"per_read={skipped_per_read_missing_motif}"
            ),
            file=sys.stderr,
        )
    if skipped_not_in_summary:
        print(
            f"warning: skipped {skipped_not_in_summary} per-read rows without matching summary locus/tool",
            file=sys.stderr,
        )
    print(
        (
            "[filtered-edit-rates-by-motif] "
            f"included_locus_tool={len(included_loci)} "
            f"excluded_locus_tool={len(excluded_loci)} "
            f"skipped_excluded_reads={skipped_excluded_reads} "
            f"out_bin={out_bin_tsv} "
            f"out_length_1_10={out_length_tsv}"
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
