#!/usr/bin/env python3
"""Summarize insertion/deletion/substitution rates by TRF motif-length bins.

The input summary files can be produced by either
``calculate_cigar_edit_rates_ont_pacbio_bed_dir_fasta.py`` or
``decompose_read_ref_edits_ont_pacbio_bed_dir_fasta.py``.  Motif period is
read from the matching original BED file, column 6 by default.  For every
tool and bin, the script averages the per-locus insertion, deletion and
substitution rates.

Three output files are written:

1. exact motif lengths: 1bp, 2bp, ..., 10bp;
2. 10bp bins: 1-10bp, 11-20bp, ..., 91-100bp;
3. 100bp bins: 1-100bp, 101-200bp, ..., 1901-2000bp.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


CoordKey = tuple[str, int, int]


@dataclass(frozen=True)
class Bin:
    label: str
    minimum: int
    maximum: int


@dataclass(frozen=True)
class LocusRate:
    motif_length: int
    tool: str
    insertion_rate: float | None
    deletion_rate: float | None
    substitution_rate: float | None


EXACT_BINS = tuple(Bin(f"{length}bp", length, length) for length in range(1, 11))
TEN_BP_BINS = tuple(
    Bin(f"{start}-{start + 9}bp", start, start + 9)
    for start in range(1, 100, 10)
)
HUNDRED_BP_BINS = tuple(
    Bin(f"{start}-{start + 99}bp", start, start + 99)
    for start in range(1, 2000, 100)
)

CIGAR_FORMAT = ("cigar", "*.cigar_edit_rates.summary.tsv", ".cigar_edit_rates.summary.tsv")
DECOMPOSE_FORMAT = ("decompose", "*.edit_ops.summary.tsv", ".edit_ops.summary.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize insertion/deletion/substitution rates by motif bins."
    )
    parser.add_argument(
        "--summary-dir",
        required=True,
        type=Path,
        help="Directory containing CIGAR or decompose summary TSV files.",
    )
    parser.add_argument(
        "--bed-dir",
        required=True,
        type=Path,
        help="Directory containing the matching original TRF BED files.",
    )
    parser.add_argument(
        "--summary-pattern",
        help="Optional summary glob pattern. If omitted, both supported formats are searched.",
    )
    parser.add_argument(
        "--summary-suffix",
        help="Optional suffix removed to find the matching BED basename.",
    )
    parser.add_argument(
        "--input-format",
        choices=("auto", "cigar", "decompose"),
        default="auto",
        help="Input summary format. Default: auto-detect from columns or filename.",
    )
    parser.add_argument(
        "--motif-column",
        type=int,
        default=6,
        help="1-based motif period column in the BED file. Default: 6.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search summary and BED directories recursively.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for the three output TSV files.",
    )
    return parser.parse_args()


def parse_number(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value.upper() == "NA" or value.lower() == "nan":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int(value: str, path: Path, column: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"{path}: invalid integer in {column}: {value!r}") from exc


def require_columns(path: Path, fieldnames: list[str] | None, required: set[str]) -> None:
    missing = required.difference(fieldnames or [])
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")


def read_bed_motif_lengths(path: Path, motif_column: int) -> dict[CoordKey, int]:
    if motif_column < 1:
        raise ValueError("--motif-column must be >= 1")
    motif_index = motif_column - 1
    result: dict[CoordKey, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if fields[0].lower() in {"track", "browser"}:
                continue
            if len(fields) <= motif_index:
                raise ValueError(f"{path}:{line_no} has fewer than {motif_column} columns")
            try:
                key = (fields[0], int(fields[1]), int(fields[2]))
                motif_length = int(float(fields[motif_index]))
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no} has invalid BED values") from exc
            if motif_length < 1:
                raise ValueError(f"{path}:{line_no} has motif length < 1")
            previous = result.get(key)
            if previous is not None and previous != motif_length:
                raise ValueError(
                    f"{path}:{line_no} has conflicting motif lengths for {key}: "
                    f"{previous} and {motif_length}"
                )
            result[key] = motif_length
    return result


def find_matching_bed(
    summary_path: Path,
    bed_dir: Path,
    summary_suffix: str,
    recursive: bool,
) -> Path | None:
    if summary_suffix and summary_path.name.endswith(summary_suffix):
        base = summary_path.name[: -len(summary_suffix)]
    else:
        base = summary_path.stem
    direct = bed_dir / f"{base}.bed"
    if direct.is_file():
        return direct
    if recursive:
        matches = sorted(bed_dir.rglob(f"{base}.bed"))
        if matches:
            return matches[0]
    return None


def read_summary_rates(
    summary_path: Path,
    bed_path: Path,
    motif_column: int,
    input_format: str,
) -> tuple[list[LocusRate], int]:
    motif_lengths = read_bed_motif_lengths(bed_path, motif_column)
    common_required = {
        "scope",
        "chrom",
        "start",
        "end",
        "tool",
        "n_reads",
    }
    if input_format == "cigar":
        rate_columns = (
            "mean_insertion_rate",
            "mean_deletion_rate",
            "mean_substitution_rate",
        )
    elif input_format == "decompose":
        rate_columns = (
            "mean_ins_per_ref_bp",
            "mean_del_per_ref_bp",
            "mean_sub_per_ref_bp",
        )
    else:
        raise ValueError(f"Unsupported input format: {input_format}")

    required = common_required.union(rate_columns)
    rates: list[LocusRate] = []
    skipped = 0
    with summary_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(summary_path, reader.fieldnames, required)
        for row in reader:
            if row.get("scope") != "locus":
                continue
            n_reads = parse_number(row.get("n_reads"))
            if n_reads is None or n_reads <= 0:
                skipped += 1
                continue
            chrom = row["chrom"]
            start = parse_int(row["start"], summary_path, "start")
            end = parse_int(row["end"], summary_path, "end")
            motif_length = motif_lengths.get((chrom, start, end))
            if motif_length is None:
                raise ValueError(
                    f"{summary_path}: locus {chrom}:{start}-{end} was not found in {bed_path}"
                )
            rates.append(
                LocusRate(
                    motif_length=motif_length,
                    tool=row["tool"],
                    insertion_rate=parse_number(row.get(rate_columns[0])),
                    deletion_rate=parse_number(row.get(rate_columns[1])),
                    substitution_rate=parse_number(row.get(rate_columns[2])),
                )
            )
    return rates, skipped


def detect_input_format(path: Path) -> str | None:
    if path.name.endswith(CIGAR_FORMAT[2]):
        return CIGAR_FORMAT[0]
    if path.name.endswith(DECOMPOSE_FORMAT[2]):
        return DECOMPOSE_FORMAT[0]
    return None


def assign_bin(motif_length: int, bins: tuple[Bin, ...]) -> Bin | None:
    for bin_definition in bins:
        if bin_definition.minimum <= motif_length <= bin_definition.maximum:
            return bin_definition
    return None


def tool_order(tools: set[str]) -> list[str]:
    preferred = ["ONT", "PacBio"]
    return [tool for tool in preferred if tool in tools] + sorted(tools.difference(preferred))


def fmt(value: float | None) -> str:
    return "NA" if value is None else f"{value:.8f}"


def average(values: list[float | None]) -> float | None:
    usable = [value for value in values if value is not None]
    return statistics.fmean(usable) if usable else None


def write_scheme(
    output_path: Path,
    bins: tuple[Bin, ...],
    rates: list[LocusRate],
) -> None:
    grouped: dict[tuple[str, str], list[LocusRate]] = defaultdict(list)
    tools = {item.tool for item in rates}
    for item in rates:
        bin_definition = assign_bin(item.motif_length, bins)
        if bin_definition is not None:
            grouped[(bin_definition.label, item.tool)].append(item)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "motif_bin",
            "motif_min_bp",
            "motif_max_bp",
            "tool",
            "n_loci",
            "n_loci_with_insertion_rate",
            "n_loci_with_deletion_rate",
            "n_loci_with_substitution_rate",
            "mean_insertion_rate",
            "mean_deletion_rate",
            "mean_substitution_rate",
        ])
        for bin_definition in bins:
            for tool in tool_order(tools):
                rows = grouped.get((bin_definition.label, tool), [])
                writer.writerow([
                    bin_definition.label,
                    bin_definition.minimum,
                    bin_definition.maximum,
                    tool,
                    len(rows),
                    sum(item.insertion_rate is not None for item in rows),
                    sum(item.deletion_rate is not None for item in rows),
                    sum(item.substitution_rate is not None for item in rows),
                    fmt(average([item.insertion_rate for item in rows])),
                    fmt(average([item.deletion_rate for item in rows])),
                    fmt(average([item.substitution_rate for item in rows])),
                ])


def main() -> int:
    args = parse_args()
    if not args.summary_dir.is_dir():
        raise SystemExit(f"Summary directory does not exist: {args.summary_dir}")
    if not args.bed_dir.is_dir():
        raise SystemExit(f"BED directory does not exist: {args.bed_dir}")

    globber = args.summary_dir.rglob if args.recursive else args.summary_dir.glob
    if args.summary_pattern:
        patterns = [args.summary_pattern]
    elif args.input_format == "cigar":
        patterns = [CIGAR_FORMAT[1]]
    elif args.input_format == "decompose":
        patterns = [DECOMPOSE_FORMAT[1]]
    else:
        patterns = [CIGAR_FORMAT[1], DECOMPOSE_FORMAT[1]]

    summary_paths = sorted({
        path
        for pattern in patterns
        for path in globber(pattern)
        if path.is_file()
    })
    if not summary_paths:
        raise SystemExit(
            f"No summary files matched {patterns!r} in {args.summary_dir}"
        )

    all_rates: list[LocusRate] = []
    missing_beds: list[Path] = []
    skipped_rows = 0
    for summary_path in summary_paths:
        input_format = args.input_format
        detected_format = detect_input_format(summary_path)
        if input_format == "auto":
            input_format = detected_format or "cigar"
        if args.summary_suffix:
            summary_suffix = args.summary_suffix
        elif input_format == "cigar":
            summary_suffix = CIGAR_FORMAT[2]
        else:
            summary_suffix = DECOMPOSE_FORMAT[2]
        bed_path = find_matching_bed(
            summary_path, args.bed_dir, summary_suffix, args.recursive
        )
        if bed_path is None:
            missing_beds.append(summary_path)
            continue
        rates, skipped = read_summary_rates(
            summary_path, bed_path, args.motif_column, input_format
        )
        all_rates.extend(rates)
        skipped_rows += skipped

    if missing_beds:
        examples = ", ".join(str(path) for path in missing_beds[:3])
        raise SystemExit(
            "Could not find matching BED files for summary files, for example: "
            f"{examples}. Check --bed-dir and --summary-suffix."
        )
    if not all_rates:
        raise SystemExit("No usable locus rows were found in the summary files")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = [
        ("mean_cigar_edit_rates_exact_1_10bp.tsv", EXACT_BINS),
        ("mean_cigar_edit_rates_bins_1_100bp.tsv", TEN_BP_BINS),
        ("mean_cigar_edit_rates_bins_1_2000bp.tsv", HUNDRED_BP_BINS),
    ]
    for filename, bins in outputs:
        write_scheme(args.output_dir / filename, bins, all_rates)
        print(f"wrote={args.output_dir / filename}")
    print(
        f"processed_summary_files={len(summary_paths)} "
        f"usable_loci={len(all_rates)} skipped_loci={skipped_rows}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
