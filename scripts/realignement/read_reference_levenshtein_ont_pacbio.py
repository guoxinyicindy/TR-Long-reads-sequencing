#!/usr/bin/env python3
"""
Calculate read-to-assembly-reference Levenshtein distances at BED loci.

Inputs:
  * maternal and paternal assembly alignments to hg38 (BAM/CRAM)
  * ONT and PacBio read BAM files
  * a BED file of hg38 loci

For each locus, the script extracts the maternal/paternal assembly sequence
over the BED interval plus an optional flank. Each read sequence is extracted
over the same target interval and compared against both haplotypes; the smaller
Levenshtein distance is used for summaries.

The script supports CPU multi-processing across loci via --processes. It does
not currently use GPU acceleration.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQ = 7
CIGAR_X = 8

REF_CONSUMING = {CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X}
QUERY_CONSUMING = {CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_EQ, CIGAR_X}
ALIGNED_QUERY_AND_REF = {CIGAR_M, CIGAR_EQ, CIGAR_X}
REF_GAP = {CIGAR_D, CIGAR_N}


@dataclass(frozen=True)
class Locus:
    chrom: str
    start: int
    end: int
    name: str


@dataclass(frozen=True)
class BamSpec:
    label: str
    path: str


@dataclass(frozen=True)
class ReadDistance:
    read_name: str
    distance: int
    matched_haplotype: str
    read_length: int
    reference_length: int
    left_flank_used: int
    right_flank_used: int


@dataclass(frozen=True)
class LocusToolResult:
    locus: Locus
    tool: str
    n_ref_alleles: int
    read_distances: list[ReadDistance]
    max_observed_read_length: int = 0


_WORKER_STATE = None


def parse_bed(path: str) -> list[Locus]:
    loci: list[Locus] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED line {line_no} has fewer than 3 columns")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                raise ValueError(f"BED line {line_no} has end <= start: {line.rstrip()}")
            name = fields[3] if len(fields) >= 4 and fields[3] else f"{chrom}:{start}-{end}"
            loci.append(Locus(chrom=chrom, start=start, end=end, name=name))
    return loci


def levenshtein_fallback(a: str, b: str) -> int:
    """Memory-efficient Levenshtein distance fallback."""
    if a == b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    if not b:
        return len(a)

    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        current = [i]
        left = i
        diag = i - 1
        for j, cb in enumerate(b, start=1):
            up = previous[j] + 1
            left = left + 1
            subst = diag + (ca != cb)
            value = min(up, left, subst)
            current.append(value)
            diag = previous[j]
            left = value
        previous = current
    return previous[-1]


def get_levenshtein_distance():
    try:
        import Levenshtein as lv  # type: ignore

        return lv.distance
    except ImportError:
        return levenshtein_fallback


def candidate_score(aln) -> tuple[int, int, int, int]:
    try:
        aln_score = int(aln.get_tag("AS"))
    except KeyError:
        aln_score = 0
    primary = 0 if aln.is_supplementary else 1
    span = (aln.reference_end or 0) - (aln.reference_start or 0)
    return primary, int(aln.mapping_quality), aln_score, span


def spans_interval(aln, start: int, end: int) -> bool:
    if aln.reference_start is None or aln.reference_end is None:
        return False
    if start == end:
        return aln.reference_start <= start <= aln.reference_end
    return aln.reference_start <= start and aln.reference_end >= end


def collect_query_sequence(
    aln,
    start: int,
    end: int,
) -> str | None:
    """Extract query bases corresponding to a reference interval.

    Insertions whose anchor position falls inside the half-open interval are
    included.
    """
    if aln.cigartuples is None or aln.query_sequence is None:
        return None
    if not spans_interval(aln, start, end):
        return None

    ref_pos = aln.reference_start
    query_pos = 0
    pieces: list[str] = []

    for op, length in aln.cigartuples:
        if op in ALIGNED_QUERY_AND_REF:
            ref_next = ref_pos + length
            overlap_start = max(start, ref_pos)
            overlap_end = min(end, ref_next)
            if overlap_start < overlap_end:
                q_start = query_pos + (overlap_start - ref_pos)
                q_end = query_pos + (overlap_end - ref_pos)
                pieces.append(aln.query_sequence[q_start:q_end])
            ref_pos = ref_next
            query_pos += length
            continue

        if op == CIGAR_I:
            if start <= ref_pos < end:
                pieces.append(aln.query_sequence[query_pos : query_pos + length])
            query_pos += length
            continue

        if op in REF_GAP:
            ref_pos += length
            continue

        if op == CIGAR_S:
            query_pos += length
            continue

        if op in {CIGAR_H, CIGAR_P}:
            continue

        if op in QUERY_CONSUMING:
            query_pos += length
        if op in REF_CONSUMING:
            ref_pos += length

    return "".join(pieces).upper()


def get_original_read_name(aln) -> str:
    try:
        return aln.get_tag("ZO")
    except KeyError:
        return aln.query_name


def should_skip_read(aln, args: argparse.Namespace) -> bool:
    if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail:
        return True
    if args.primary_only and (aln.is_secondary or aln.is_supplementary):
        return True
    if not args.include_secondary and aln.is_secondary:
        return True
    if aln.mapping_quality < args.min_mapq:
        return True
    return False


def fetch_assembly_haplotype(aln_file, locus: Locus, left_flank: int, right_flank: int) -> str | None:
    start = max(0, locus.start - left_flank)
    end = locus.end + right_flank
    candidates = []
    for aln in aln_file.fetch(locus.chrom, start, end):
        if aln.is_unmapped or aln.is_secondary or aln.is_duplicate or aln.is_qcfail:
            continue
        seq = collect_query_sequence(aln, start, end)
        if seq is None:
            continue
        candidates.append((candidate_score(aln), seq))

    if not candidates:
        return None
    _, best_seq = max(candidates, key=lambda item: item[0])
    return best_seq


def fetch_reference_alleles(mat_aln, pat_aln, locus: Locus, left_flank: int, right_flank: int) -> dict[str, str]:
    alleles: dict[str, str] = {}
    mat_seq = fetch_assembly_haplotype(mat_aln, locus, left_flank, right_flank)
    if mat_seq is not None:
        alleles["maternal"] = mat_seq
    pat_seq = fetch_assembly_haplotype(pat_aln, locus, left_flank, right_flank)
    if pat_seq is not None:
        alleles["paternal"] = pat_seq
    return alleles


def extract_read_sequence(
    aln,
    locus: Locus,
    left_flank: int,
    right_flank: int,
) -> str | None:
    start = max(0, locus.start - left_flank)
    end = locus.end + right_flank
    return collect_query_sequence(aln, start, end)


def iter_distances_for_locus(
    bam,
    spec: BamSpec,
    locus: Locus,
    ref_alleles: dict[str, str],
    flank: int,
    distance_fn,
    args: argparse.Namespace,
) -> tuple[list[ReadDistance], int]:
    if locus.chrom not in bam.references:
        return [], 0

    fetch_start = max(0, locus.start - flank - args.fetch_pad)
    fetch_end = max(locus.end + flank + args.fetch_pad, locus.start + flank + args.fetch_pad + 1)
    best_by_read: dict[str, ReadDistance] = {}
    max_observed_read_length = 0

    for aln in bam.fetch(locus.chrom, fetch_start, fetch_end):
        if should_skip_read(aln, args):
            continue
        left_flank_used = flank
        right_flank_used = flank
        alleles_for_read = ref_alleles
        seq = extract_read_sequence(aln, locus, left_flank_used, right_flank_used)

        if not seq:
            continue
        max_observed_read_length = max(max_observed_read_length, len(seq))

        distances = [
            (hap, distance_fn(seq, ref_seq), len(ref_seq))
            for hap, ref_seq in alleles_for_read.items()
        ]
        matched_hap, distance, ref_len = min(distances, key=lambda item: item[1])
        if abs(len(seq) - ref_len) > args.max_length_difference:
            continue
        read_name = get_original_read_name(aln)
        item = ReadDistance(
            read_name=read_name,
            distance=int(distance),
            matched_haplotype=matched_hap,
            read_length=len(seq),
            reference_length=ref_len,
            left_flank_used=left_flank_used,
            right_flank_used=right_flank_used,
        )
        previous = best_by_read.get(read_name)
        if previous is None or item.distance < previous.distance:
            best_by_read[read_name] = item

    return list(best_by_read.values()), max_observed_read_length


def fmt(value: float | int | None) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{value:.6f}"


def mean_or_none(values: list[int]) -> float | None:
    return statistics.fmean(values) if values else None


def mean_float_or_none(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def median_or_none(values: list[int]) -> float | None:
    return statistics.median(values) if values else None


def write_summary_header(writer: csv.writer) -> None:
    writer.writerow(
        [
            "scope",
            "locus_id",
            "chrom",
            "start",
            "end",
            "tool",
            "n_ref_alleles",
            "n_loci_with_reads",
            "n_reads",
            "mean_distance",
            "median_distance",
            "min_distance",
            "max_distance",
            "mean_read_length",
            "mean_reference_length",
            "mean_left_flank_used",
            "mean_right_flank_used",
            "mean_locus_mean_distance",
        ]
    )


def write_summary_row(
    writer: csv.writer,
    scope: str,
    locus: Locus | None,
    tool: str,
    n_ref_alleles: int,
    distances: list[int],
    read_lengths: list[int],
    ref_lengths: list[int],
    left_flanks: list[int],
    right_flanks: list[int],
    n_loci_with_reads: int,
    mean_locus_mean: float | None = None,
) -> None:
    writer.writerow(
        [
            scope,
            locus.name if locus else "ALL",
            locus.chrom if locus else "ALL",
            locus.start if locus else "NA",
            locus.end if locus else "NA",
            tool,
            n_ref_alleles,
            n_loci_with_reads,
            len(distances),
            fmt(mean_or_none(distances)),
            fmt(median_or_none(distances)),
            fmt(min(distances) if distances else None),
            fmt(max(distances) if distances else None),
            fmt(mean_or_none(read_lengths)),
            fmt(mean_or_none(ref_lengths)),
            fmt(mean_or_none(left_flanks)),
            fmt(mean_or_none(right_flanks)),
            fmt(mean_locus_mean),
        ]
    )


def write_per_read_header(writer: csv.writer) -> None:
    writer.writerow(
        [
            "locus_id",
            "chrom",
            "start",
            "end",
            "tool",
            "read_name",
            "distance",
            "matched_haplotype",
            "read_length",
            "reference_length",
            "left_flank_used",
            "right_flank_used",
        ]
    )


def process_locus(
    locus: Locus,
    args: argparse.Namespace,
    mat_aln,
    pat_aln,
    bam_by_label: dict[str, object],
    bam_specs: list[BamSpec],
    distance_fn,
) -> list[LocusToolResult]:
    if locus.chrom not in mat_aln.references and locus.chrom not in pat_aln.references:
        return []

    ref_alleles = fetch_reference_alleles(mat_aln, pat_aln, locus, args.flank, args.flank)
    n_ref_alleles = len(ref_alleles)
    results: list[LocusToolResult] = []
    for spec in bam_specs:
        read_distances: list[ReadDistance] = []
        max_observed_read_length = 0
        if ref_alleles:
            read_distances, max_observed_read_length = iter_distances_for_locus(
                bam_by_label[spec.label],
                spec,
                locus,
                ref_alleles,
                args.flank,
                distance_fn,
                args,
            )
        results.append(
            LocusToolResult(
                locus=locus,
                tool=spec.label,
                n_ref_alleles=n_ref_alleles,
                read_distances=read_distances,
                max_observed_read_length=max_observed_read_length,
            )
        )
    return results


def record_locus_results(
    writer: csv.writer,
    per_read_writer,
    tool_results: list[LocusToolResult],
    totals: dict[str, list[ReadDistance]],
    locus_means_by_tool: dict[str, list[float]],
    ref_allele_counts_by_tool: dict[str, int],
) -> None:
    for result in tool_results:
        ref_allele_counts_by_tool[result.tool] += result.n_ref_alleles
        distances = [item.distance for item in result.read_distances]
        read_lengths = [item.read_length for item in result.read_distances]
        ref_lengths = [item.reference_length for item in result.read_distances]
        left_flanks = [item.left_flank_used for item in result.read_distances]
        right_flanks = [item.right_flank_used for item in result.read_distances]
        if distances:
            locus_mean = statistics.fmean(distances)
            locus_means_by_tool[result.tool].append(locus_mean)
        else:
            locus_mean = None

        write_summary_row(
            writer,
            "locus",
            result.locus,
            result.tool,
            result.n_ref_alleles,
            distances,
            read_lengths,
            ref_lengths,
            left_flanks,
            right_flanks,
            1 if distances else 0,
            locus_mean,
        )
        totals[result.tool].extend(result.read_distances)

        if per_read_writer is not None:
            for item in result.read_distances:
                per_read_writer.writerow(
                    [
                        result.locus.name,
                        result.locus.chrom,
                        result.locus.start,
                        result.locus.end,
                        result.tool,
                        item.read_name,
                        item.distance,
                        item.matched_haplotype,
                        item.read_length,
                        item.reference_length,
                        item.left_flank_used,
                        item.right_flank_used,
                    ]
                )


def tool_read_counts(tool_results: list[LocusToolResult]) -> dict[str, int]:
    return {result.tool: len(result.read_distances) for result in tool_results}


def tool_max_kept_lengths(tool_results: list[LocusToolResult]) -> dict[str, int]:
    return {
        result.tool: max((item.read_length for item in result.read_distances), default=0)
        for result in tool_results
    }


def tool_max_observed_lengths(tool_results: list[LocusToolResult]) -> dict[str, int]:
    return {result.tool: result.max_observed_read_length for result in tool_results}


def format_read_counts(tool_results: list[LocusToolResult]) -> str:
    counts = tool_read_counts(tool_results)
    total = sum(counts.values())
    return (
        f"reads_total={total} "
        f"ONT={counts.get('ONT', 0)} "
        f"PacBio={counts.get('PacBio', 0)}"
    )


def format_length_diagnostics(locus: Locus, tool_results: list[LocusToolResult], flank: int) -> str:
    max_kept = tool_max_kept_lengths(tool_results)
    max_observed = tool_max_observed_lengths(tool_results)
    locus_len = locus.end - locus.start
    target_len = locus_len + 2 * flank
    return (
        f"locus_len={locus_len} target_len={target_len} "
        f"max_kept_read_len="
        f"ONT:{max_kept.get('ONT', 0)},"
        f"PacBio:{max_kept.get('PacBio', 0)} "
        f"max_observed_read_len="
        f"ONT:{max_observed.get('ONT', 0)},"
        f"PacBio:{max_observed.get('PacBio', 0)}"
    )


def print_locus_progress(
    idx: int,
    total_loci: int,
    locus: Locus,
    flank: int,
    stage: str = "DONE",
    tool_results: list[LocusToolResult] | None = None,
    elapsed: float | None = None,
    completed: int | None = None,
    written: int | None = None,
) -> None:
    parts = [
        f"[read-reference-levenshtein] {stage} locus {idx}/{total_loci}",
        f"{locus.chrom}:{locus.start}-{locus.end}",
    ]
    if tool_results is None:
        locus_len = locus.end - locus.start
        parts.append(f"locus_len={locus_len}")
        parts.append(f"target_len={locus_len + 2 * flank}")
    if tool_results is not None:
        parts.append(format_read_counts(tool_results))
        parts.append(format_length_diagnostics(locus, tool_results, flank))
    if elapsed is not None:
        parts.append(f"elapsed={elapsed:.2f}s")
    if completed is not None and written is not None:
        parts.append(f"completed={completed}/{total_loci}")
        parts.append(f"written={written}/{total_loci}")
    print(" ".join(parts), file=sys.stderr, flush=True)


def init_worker(args: argparse.Namespace) -> None:
    global _WORKER_STATE
    import pysam

    bam_specs = [
        BamSpec("ONT", args.ont_bam),
        BamSpec("PacBio", args.pacbio_bam),
    ]
    mat_aln = pysam.AlignmentFile(args.maternal_aln, "rb")
    pat_aln = pysam.AlignmentFile(args.paternal_aln, "rb")
    ont_bam = pysam.AlignmentFile(args.ont_bam, "rb")
    pacbio_bam = pysam.AlignmentFile(args.pacbio_bam, "rb")
    _WORKER_STATE = {
        "args": args,
        "mat_aln": mat_aln,
        "pat_aln": pat_aln,
        "bam_by_label": {
            "ONT": ont_bam,
            "PacBio": pacbio_bam,
        },
        "bam_specs": bam_specs,
        "distance_fn": get_levenshtein_distance(),
    }


def process_locus_worker(job: tuple[int, Locus]) -> tuple[int, list[LocusToolResult]]:
    index, locus = job
    state = _WORKER_STATE
    if state is None:
        raise RuntimeError("Worker state not initialized")
    results = process_locus(
        locus,
        state["args"],
        state["mat_aln"],
        state["pat_aln"],
        state["bam_by_label"],
        state["bam_specs"],
        state["distance_fn"],
    )
    return index, results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize Levenshtein distances between ONT/PacBio reads and maternal/paternal assembly reference alleles at BED loci."
    )
    parser.add_argument("--bed", required=True, help="BED file with 0-based half-open hg38 loci.")
    parser.add_argument("--maternal-aln", required=True, help="Maternal assembly alignment to hg38, BAM/CRAM.")
    parser.add_argument("--paternal-aln", required=True, help="Paternal assembly alignment to hg38, BAM/CRAM.")
    parser.add_argument("--ont-bam", required=True, help="ONT reads aligned to hg38, BAM/CRAM.")
    parser.add_argument("--pacbio-bam", required=True, help="PacBio reads aligned to hg38, BAM/CRAM.")
    parser.add_argument("-o", "--output-tsv", required=True, help="Output summary TSV.")
    parser.add_argument("--per-read-tsv", help="Optional per-read distance TSV.")
    parser.add_argument(
        "--flank",
        type=int,
        default=0,
        help="Bases of flanking sequence to include on both sides for all comparisons. Default: 0.",
    )
    parser.add_argument(
        "--fetch-pad",
        type=int,
        default=50,
        help="Extra reference bases for fetching reads around each target interval. Default: 50.",
    )
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum read mapping quality. Default: 0.")
    parser.add_argument(
        "--max-length-difference",
        type=int,
        default=100,
        help=(
            "Keep only reads with abs(read_length - matched_reference_length) <= this value. "
            "Default: 100. Use a negative value to disable this filter."
        ),
    )
    parser.add_argument(
        "--primary-only",
        action="store_true",
        help="Skip both secondary and supplementary read alignments.",
    )
    parser.add_argument(
        "--include-secondary",
        action="store_true",
        help="Include secondary read alignments. Supplementary alignments are included unless --primary-only is set.",
    )
    parser.add_argument(
        "--chrom",
        action="append",
        help="Only process this chromosome/contig. Can be provided multiple times.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="Print progress every N loci. Set 0 to disable. Default: 100.",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of CPU worker processes to use across loci. Default: 1.",
    )
    parser.add_argument(
        "--max-locus-length",
        type=int,
        default=1000,
        help=(
            "Skip loci with (end - start) greater than this value before any read processing. "
            "Default: 1000. Use a negative value to disable this filter."
        ),
    )
    return parser


def validate_paths(args: argparse.Namespace) -> None:
    for attr in ("bed", "maternal_aln", "paternal_aln", "ont_bam", "pacbio_bam"):
        path = Path(getattr(args, attr))
        if not path.exists():
            raise SystemExit(f"Input path does not exist: {path}")


def run(args: argparse.Namespace) -> None:
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit("Missing dependency: pysam. Activate/install the project environment first.") from exc

    validate_paths(args)
    loci = parse_bed(args.bed)
    if args.chrom:
        chroms = set(args.chrom)
        loci = [locus for locus in loci if locus.chrom in chroms]
        if not loci:
            raise SystemExit(f"No BED loci matched --chrom: {', '.join(args.chrom)}")
    if args.flank < 0:
        raise SystemExit("--flank must be >= 0")
    if args.max_length_difference < 0:
        args.max_length_difference = sys.maxsize
    if args.max_locus_length < 0:
        args.max_locus_length = sys.maxsize
    if args.processes < 1:
        raise SystemExit("--processes must be >= 1")

    total_loci_before_length_filter = len(loci)
    loci = [locus for locus in loci if (locus.end - locus.start) <= args.max_locus_length]
    skipped_long_loci = total_loci_before_length_filter - len(loci)
    if not loci:
        raise SystemExit(
            f"All loci were filtered out by --max-locus-length {args.max_locus_length}."
        )

    bam_specs = [
        BamSpec("ONT", args.ont_bam),
        BamSpec("PacBio", args.pacbio_bam),
    ]

    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    per_read_handle = None
    per_read_writer = None
    if args.per_read_tsv:
        per_read_path = Path(args.per_read_tsv)
        per_read_path.parent.mkdir(parents=True, exist_ok=True)
        per_read_handle = open(per_read_path, "w", encoding="utf-8", newline="")
        per_read_writer = csv.writer(per_read_handle, delimiter="\t", lineterminator="\n")
        write_per_read_header(per_read_writer)

    totals: dict[str, list[ReadDistance]] = defaultdict(list)
    locus_means_by_tool: dict[str, list[float]] = defaultdict(list)
    ref_allele_counts_by_tool: dict[str, int] = defaultdict(int)

    try:
        with (
            pysam.AlignmentFile(args.maternal_aln, "rb") as mat_aln,
            pysam.AlignmentFile(args.paternal_aln, "rb") as pat_aln,
            pysam.AlignmentFile(args.ont_bam, "rb") as ont_bam,
            pysam.AlignmentFile(args.pacbio_bam, "rb") as pacbio_bam,
            open(output_path, "w", encoding="utf-8", newline="") as out_handle,
        ):
            bam_by_label = {
                "ONT": ont_bam,
                "PacBio": pacbio_bam,
            }
            distance_fn = get_levenshtein_distance()
            writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
            write_summary_header(writer)

            if args.processes == 1:
                for idx, locus in enumerate(loci, start=1):
                    if args.progress_every and (idx == 1 or idx == len(loci) or idx % args.progress_every == 0):
                        print_locus_progress(idx, len(loci), locus, args.flank, stage="START")
                    started_at = time.monotonic()
                    tool_results = process_locus(
                        locus,
                        args,
                        mat_aln,
                        pat_aln,
                        bam_by_label,
                        bam_specs,
                        distance_fn,
                    )
                    elapsed = time.monotonic() - started_at
                    if args.progress_every and (idx == 1 or idx == len(loci) or idx % args.progress_every == 0):
                        print_locus_progress(
                            idx,
                            len(loci),
                            locus,
                            args.flank,
                            stage="DONE",
                            tool_results=tool_results,
                            elapsed=elapsed,
                        )
                    record_locus_results(
                        writer,
                        per_read_writer,
                        tool_results,
                        totals,
                        locus_means_by_tool,
                        ref_allele_counts_by_tool,
                    )
            else:
                pending_results: dict[int, list[LocusToolResult]] = {}
                next_to_write = 0
                completed_jobs = 0
                jobs = list(enumerate(loci, start=1))
                with ProcessPoolExecutor(
                    max_workers=args.processes,
                    initializer=init_worker,
                    initargs=(args,),
                ) as executor:
                    future_map = {
                        executor.submit(process_locus_worker, job): job[0]
                        for job in jobs
                    }
                    for future in as_completed(future_map):
                        idx, tool_results = future.result()
                        completed_jobs += 1
                        pending_results[idx] = tool_results
                        while next_to_write + 1 in pending_results:
                            next_to_write += 1
                            locus = loci[next_to_write - 1]
                            ordered_tool_results = pending_results.pop(next_to_write)
                            if args.progress_every and (
                                next_to_write == 1
                                or next_to_write == len(loci)
                                or next_to_write % args.progress_every == 0
                            ):
                                print_locus_progress(
                                    next_to_write,
                                    len(loci),
                                    locus,
                                    args.flank,
                                    stage="DONE",
                                    tool_results=ordered_tool_results,
                                    completed=completed_jobs,
                                    written=next_to_write,
                                )
                            record_locus_results(
                                writer,
                                per_read_writer,
                                ordered_tool_results,
                                totals,
                                locus_means_by_tool,
                                ref_allele_counts_by_tool,
                            )

            all_tool_distances: list[ReadDistance] = []
            all_tool_locus_means: list[float] = []
            for spec in bam_specs:
                items = totals[spec.label]
                all_tool_distances.extend(items)
                all_tool_locus_means.extend(locus_means_by_tool[spec.label])
                distances = [item.distance for item in items]
                read_lengths = [item.read_length for item in items]
                ref_lengths = [item.reference_length for item in items]
                left_flanks = [item.left_flank_used for item in items]
                right_flanks = [item.right_flank_used for item in items]
                write_summary_row(
                    writer,
                    "tool_total",
                    None,
                    spec.label,
                    ref_allele_counts_by_tool[spec.label],
                    distances,
                    read_lengths,
                    ref_lengths,
                    left_flanks,
                    right_flanks,
                    len(locus_means_by_tool[spec.label]),
                    mean_float_or_none(locus_means_by_tool[spec.label]),
                )

            all_distances = [item.distance for item in all_tool_distances]
            all_read_lengths = [item.read_length for item in all_tool_distances]
            all_ref_lengths = [item.reference_length for item in all_tool_distances]
            all_left_flanks = [item.left_flank_used for item in all_tool_distances]
            all_right_flanks = [item.right_flank_used for item in all_tool_distances]
            write_summary_row(
                writer,
                "all_tools_total",
                None,
                "ALL",
                sum(ref_allele_counts_by_tool.values()),
                all_distances,
                all_read_lengths,
                all_ref_lengths,
                all_left_flanks,
                all_right_flanks,
                len(all_tool_locus_means),
                statistics.fmean(all_tool_locus_means) if all_tool_locus_means else None,
            )
            print(
                (
                    "[read-reference-levenshtein] done "
                    f"processed_loci={len(loci)} "
                    f"skipped_long_loci={skipped_long_loci} "
                    f"max_locus_length={args.max_locus_length}"
                ),
                file=sys.stderr,
                flush=True,
            )
    finally:
        if per_read_handle is not None:
            per_read_handle.close()


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
