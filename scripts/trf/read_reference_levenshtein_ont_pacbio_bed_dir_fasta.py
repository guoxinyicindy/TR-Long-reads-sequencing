#!/usr/bin/env python3
"""Calculate read-to-HG002-FASTA Levenshtein distances for BED files.

This version expects BED coordinates and ONT/PacBio BAM coordinates to be in
the same HG002 FASTA coordinate system. Maternal and paternal reference
sequences are fetched directly from FASTA files; no hg38 liftover BAM is used.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import importlib.machinery
import re
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
HAP_SUFFIX_RE = re.compile(r"_(MATERNAL|PATERNAL)$", re.IGNORECASE)


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
    candidate_read_count: int | None = None
    status: str = "processed"
    skip_reason: str | None = None


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
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval at line {line_no}: {line.rstrip()}")
            name = fields[3] if len(fields) >= 4 and fields[3] else f"{chrom}:{start}-{end}"
            loci.append(Locus(chrom=chrom, start=start, end=end, name=name))
    return loci


def haplotype_from_chrom(chrom: str) -> str | None:
    match = HAP_SUFFIX_RE.search(chrom)
    return match.group(1).lower() if match else None


def base_chrom(chrom: str) -> str:
    return HAP_SUFFIX_RE.sub("", chrom)


def contig_candidates(chrom: str, haplotype: str | None = None) -> list[str]:
    """Return possible FASTA/BAM contig names without guessing a wrong haplotype."""
    explicit_haplotype = haplotype_from_chrom(chrom)
    base = base_chrom(chrom)
    candidates: list[str] = [chrom]

    if explicit_haplotype is not None:
        if haplotype is not None and explicit_haplotype != haplotype:
            return []
        candidates.append(base)
    else:
        candidates.append(base)
        if haplotype is not None:
            candidates.append(f"{base}_{haplotype.upper()}")

    return list(dict.fromkeys(candidates))


def resolve_contig(file_handle, chrom: str, haplotype: str | None = None) -> str | None:
    references = set(file_handle.references)
    for candidate in contig_candidates(chrom, haplotype):
        if candidate in references:
            return candidate
    return None


def get_levenshtein_distance():
    try:
        import Levenshtein as lv  # type: ignore
        backend = importlib.import_module("Levenshtein.levenshtein_cpp")
    except ImportError as exc:
        raise RuntimeError(
            "The compiled Levenshtein extension is required. "
            "Install python-Levenshtein in the same environment used to run this script."
        ) from exc

    extension_suffixes = tuple(importlib.machinery.EXTENSION_SUFFIXES)
    backend_path = getattr(backend, "__file__", "")
    if not backend_path or not backend_path.endswith(extension_suffixes):
        raise RuntimeError(
            "Levenshtein.levenshtein_cpp is not a compiled C/C++ extension. "
            "Install python-Levenshtein in the same environment used to run this script."
        )
    return lv.distance


def spans_interval(aln, start: int, end: int) -> bool:
    if aln.reference_start is None or aln.reference_end is None:
        return False
    return aln.reference_start <= start and aln.reference_end >= end


def collect_query_sequence(aln, start: int, end: int) -> str | None:
    """Extract query bases corresponding to a reference interval."""
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
                q_start = query_pos + overlap_start - ref_pos
                q_end = query_pos + overlap_end - ref_pos
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


def fetch_fasta_sequence(fasta, chrom: str, haplotype: str, start: int, end: int) -> str | None:
    contig = resolve_contig(fasta, chrom, haplotype)
    if contig is None:
        return None
    contig_length = fasta.get_reference_length(contig)
    if start >= contig_length or end <= 0:
        return None
    clipped_start = max(0, start)
    clipped_end = min(contig_length, end)
    if clipped_end <= clipped_start:
        return None
    return fasta.fetch(contig, clipped_start, clipped_end).upper()


def fetch_reference_alleles(
    maternal_fa,
    paternal_fa,
    locus: Locus,
    left_flank: int,
    right_flank: int,
) -> dict[str, str]:
    start = max(0, locus.start - left_flank)
    end = locus.end + right_flank
    locus_haplotype = haplotype_from_chrom(locus.chrom)
    haplotypes = [locus_haplotype] if locus_haplotype else ["maternal", "paternal"]
    fasta_by_haplotype = {
        "maternal": maternal_fa,
        "paternal": paternal_fa,
    }

    alleles: dict[str, str] = {}
    for haplotype in haplotypes:
        sequence = fetch_fasta_sequence(
            fasta_by_haplotype[haplotype],
            locus.chrom,
            haplotype,
            start,
            end,
        )
        if sequence is not None:
            alleles[haplotype] = sequence
    return alleles


def extract_read_sequence(aln, locus: Locus, left_flank: int, right_flank: int) -> str | None:
    start = max(0, locus.start - left_flank)
    end = locus.end + right_flank
    return collect_query_sequence(aln, start, end)


def count_candidate_reads_for_locus(bam, locus: Locus, args: argparse.Namespace) -> int:
    contig = resolve_contig(bam, locus.chrom)
    if contig is None:
        return 0
    fetch_start = max(0, locus.start - args.flank - args.fetch_pad)
    fetch_end = max(
        locus.end + args.flank + args.fetch_pad,
        locus.start + args.flank + args.fetch_pad + 1,
    )
    read_names: set[str] = set()
    for aln in bam.fetch(contig, fetch_start, fetch_end):
        if should_skip_read(aln, args):
            continue
        read_names.add(get_original_read_name(aln))
    return len(read_names)


def iter_distances_for_locus(
    bam,
    locus: Locus,
    ref_alleles: dict[str, str],
    flank: int,
    distance_fn,
    args: argparse.Namespace,
) -> tuple[list[ReadDistance], int]:
    contig = resolve_contig(bam, locus.chrom)
    if contig is None:
        return [], 0

    fetch_start = max(0, locus.start - flank - args.fetch_pad)
    fetch_end = max(locus.end + flank + args.fetch_pad, locus.start + flank + args.fetch_pad + 1)
    best_by_read: dict[str, ReadDistance] = {}
    max_observed_read_length = 0

    for aln in bam.fetch(contig, fetch_start, fetch_end):
        if should_skip_read(aln, args):
            continue
        sequence = extract_read_sequence(aln, locus, flank, flank)
        if not sequence:
            continue
        max_observed_read_length = max(max_observed_read_length, len(sequence))

        distances = [
            (haplotype, distance_fn(sequence, reference), len(reference))
            for haplotype, reference in ref_alleles.items()
        ]
        matched_haplotype, distance, reference_length = min(distances, key=lambda item: item[1])
        read_name = get_original_read_name(aln)
        item = ReadDistance(
            read_name=read_name,
            distance=int(distance),
            matched_haplotype=matched_haplotype,
            read_length=len(sequence),
            reference_length=reference_length,
            left_flank_used=flank,
            right_flank_used=flank,
        )
        previous = best_by_read.get(read_name)
        if previous is None or item.distance < previous.distance:
            best_by_read[read_name] = item

    return list(best_by_read.values()), max_observed_read_length


def process_locus(
    locus: Locus,
    args: argparse.Namespace,
    maternal_fa,
    paternal_fa,
    bam_by_label: dict[str, object],
    bam_specs: list[BamSpec],
    distance_fn,
) -> list[LocusToolResult]:
    ref_alleles = fetch_reference_alleles(
        maternal_fa,
        paternal_fa,
        locus,
        args.flank,
        args.flank,
    )
    n_ref_alleles = len(ref_alleles)
    if not ref_alleles:
        return [
            LocusToolResult(
                locus=locus,
                tool=spec.label,
                n_ref_alleles=0,
                read_distances=[],
                status="no_reference",
                skip_reason="contig_not_found_in_expected_fasta",
            )
            for spec in bam_specs
        ]

    candidate_counts: dict[str, int] = {}
    if args.max_reads_per_locus >= 0:
        for spec in bam_specs:
            candidate_counts[spec.label] = count_candidate_reads_for_locus(
                bam_by_label[spec.label], locus, args
            )
        over_limit = [
            spec.label
            for spec in bam_specs
            if candidate_counts[spec.label] > args.max_reads_per_locus
        ]
        if over_limit:
            reason = ";".join(
                f"{tool}_candidate_reads>{args.max_reads_per_locus}" for tool in over_limit
            )
            return [
                LocusToolResult(
                    locus=locus,
                    tool=spec.label,
                    n_ref_alleles=n_ref_alleles,
                    read_distances=[],
                    candidate_read_count=candidate_counts[spec.label],
                    status="skipped",
                    skip_reason=reason,
                )
                for spec in bam_specs
            ]

    results: list[LocusToolResult] = []
    for spec in bam_specs:
        distances, max_observed = iter_distances_for_locus(
            bam_by_label[spec.label],
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
                read_distances=distances,
                max_observed_read_length=max_observed,
                candidate_read_count=candidate_counts.get(spec.label),
            )
        )
    return results


def fmt(value: float | int | None) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{value:.6f}"


def mean_int(values: list[int]) -> float | None:
    return statistics.fmean(values) if values else None


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
            "n_candidate_reads",
            "locus_status",
            "skip_reason",
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
    candidate_read_count: int | None = None,
    locus_status: str = "aggregate",
    skip_reason: str | None = None,
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
            fmt(mean_int(distances)),
            fmt(statistics.median(distances) if distances else None),
            fmt(min(distances) if distances else None),
            fmt(max(distances) if distances else None),
            fmt(mean_int(read_lengths)),
            fmt(mean_int(ref_lengths)),
            fmt(mean_int(left_flanks)),
            fmt(mean_int(right_flanks)),
            fmt(mean_locus_mean),
            fmt(candidate_read_count),
            locus_status,
            skip_reason or "NA",
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


def record_locus_results(
    writer: csv.writer,
    per_read_writer,
    tool_results: list[LocusToolResult],
    totals: dict[str, list[ReadDistance]],
    locus_means: dict[str, list[float]],
    ref_counts: dict[str, int],
) -> None:
    for result in tool_results:
        distances = [item.distance for item in result.read_distances]
        read_lengths = [item.read_length for item in result.read_distances]
        ref_lengths = [item.reference_length for item in result.read_distances]
        left_flanks = [item.left_flank_used for item in result.read_distances]
        right_flanks = [item.right_flank_used for item in result.read_distances]
        locus_mean = statistics.fmean(distances) if distances else None
        if locus_mean is not None:
            locus_means[result.tool].append(locus_mean)
        ref_counts[result.tool] += result.n_ref_alleles

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
            result.candidate_read_count,
            result.status,
            result.skip_reason,
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


def print_progress(
    idx: int,
    total: int,
    locus: Locus,
    results: list[LocusToolResult] | None = None,
    completed: int | None = None,
    written: int | None = None,
) -> None:
    message = f"[read-reference-levenshtein-fasta] locus {idx}/{total} {locus.chrom}:{locus.start}-{locus.end}"
    if results is not None:
        counts = ",".join(f"{r.tool}:{len(r.read_distances)}" for r in results)
        message += f" reads={counts}"
        skipped = [f"{r.tool}:{r.candidate_read_count}" for r in results if r.status == "skipped"]
        if skipped:
            message += " skipped_candidate_reads=" + ",".join(skipped)
    if completed is not None and written is not None:
        message += f" completed={completed}/{total} written={written}/{total}"
    print(message, file=sys.stderr, flush=True)


def init_worker(args: argparse.Namespace) -> None:
    global _WORKER_STATE
    import pysam

    _WORKER_STATE = {
        "args": args,
        "maternal_fa": pysam.FastaFile(args.maternal_fa),
        "paternal_fa": pysam.FastaFile(args.paternal_fa),
        "bam_by_label": {
            "ONT": pysam.AlignmentFile(args.ont_bam, "rb"),
            "PacBio": pysam.AlignmentFile(args.pacbio_bam, "rb"),
        },
        "bam_specs": [
            BamSpec("ONT", args.ont_bam),
            BamSpec("PacBio", args.pacbio_bam),
        ],
        "distance_fn": get_levenshtein_distance(),
    }


def process_locus_worker(job: tuple[int, Locus]) -> tuple[int, list[LocusToolResult]]:
    index, locus = job
    if _WORKER_STATE is None:
        raise RuntimeError("Worker state not initialized")
    results = process_locus(
        locus,
        _WORKER_STATE["args"],
        _WORKER_STATE["maternal_fa"],
        _WORKER_STATE["paternal_fa"],
        _WORKER_STATE["bam_by_label"],
        _WORKER_STATE["bam_specs"],
        _WORKER_STATE["distance_fn"],
    )
    return index, results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate Levenshtein distances between ONT/PacBio reads and "
            "maternal/paternal HG002 FASTA sequences for every BED file."
        )
    )
    parser.add_argument("--bed-dir", required=True, type=Path, help="Directory containing BED files.")
    parser.add_argument("--bed-pattern", default="*.bed", help="BED glob pattern. Default: *.bed")
    parser.add_argument("--recursive", action="store_true", help="Search BED files recursively.")
    parser.add_argument("--maternal-fa", required=True, help="Maternal HG002 reference FASTA.")
    parser.add_argument("--paternal-fa", required=True, help="Paternal HG002 reference FASTA.")
    parser.add_argument("--ont-bam", required=True, help="ONT reads aligned to the HG002 FASTA.")
    parser.add_argument("--pacbio-bam", required=True, help="PacBio reads aligned to the HG002 FASTA.")
    parser.add_argument("-o", "--output-dir", required=True, type=Path, help="Directory for summary TSV files.")
    parser.add_argument(
        "--summary-suffix",
        default=".read_ref_lev.summary.tsv",
        help="Summary filename suffix. Default: .read_ref_lev.summary.tsv",
    )
    parser.add_argument("--write-per-read", action="store_true", help="Also write per-read distance TSV files.")
    parser.add_argument(
        "--per-read-suffix",
        default=".read_ref_lev.per_read.tsv",
        help="Per-read filename suffix. Default: .read_ref_lev.per_read.tsv",
    )
    parser.add_argument("--skip-existing", action="store_true", help="Skip BED files with existing summary output.")
    parser.add_argument("--continue-on-error", action="store_true", help="Continue if one BED file fails.")
    parser.add_argument("--flank", type=int, default=0, help="Reference/read flank bases. Default: 0.")
    parser.add_argument("--fetch-pad", type=int, default=50, help="Extra BAM fetch padding. Default: 50.")
    parser.add_argument(
        "--max-reads-per-locus",
        type=int,
        default=-1,
        help="Skip a locus if either platform has more than this many candidate reads. Default: disabled.",
    )
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum read MAPQ. Default: 0.")
    parser.add_argument("--primary-only", action="store_true", help="Skip secondary and supplementary alignments.")
    parser.add_argument(
        "--include-secondary",
        action="store_true",
        help="Include secondary alignments. Supplementary are included by default.",
    )
    parser.add_argument("--chrom", action="append", help="Process only this contig; may be repeated.")
    parser.add_argument("--progress-every", type=int, default=100, help="Progress interval. 0 disables it.")
    parser.add_argument("--processes", type=int, default=1, help="Worker processes per BED. Default: 1.")
    return parser


def validate_paths(args: argparse.Namespace) -> None:
    for path in (args.maternal_fa, args.paternal_fa, args.ont_bam, args.pacbio_bam):
        if not Path(path).exists():
            raise SystemExit(f"Input path does not exist: {path}")
    for fasta in (args.maternal_fa, args.paternal_fa):
        fai = Path(fasta + ".fai")
        if not fai.exists():
            raise SystemExit(f"Missing FASTA index: {fai}; run samtools faidx first")


def run_one_bed(args: argparse.Namespace) -> None:
    import pysam

    validate_paths(args)
    loci = parse_bed(args.bed)
    if args.chrom:
        requested = set(args.chrom)
        loci = [locus for locus in loci if locus.chrom in requested or base_chrom(locus.chrom) in requested]
    if not loci:
        raise SystemExit(f"No loci to process in {args.bed}")
    if args.flank < 0 or args.fetch_pad < 0:
        raise SystemExit("--flank and --fetch-pad must be >= 0")
    if args.max_reads_per_locus < -1:
        raise SystemExit("--max-reads-per-locus must be >= 0 or -1 to disable it")
    if args.processes < 1:
        raise SystemExit("--processes must be >= 1")
    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    per_read_handle = None
    try:
        per_read_writer = None
        if args.per_read_tsv:
            per_read_path = Path(args.per_read_tsv)
            per_read_path.parent.mkdir(parents=True, exist_ok=True)
            per_read_handle = open(per_read_path, "w", encoding="utf-8", newline="")
            per_read_writer = csv.writer(per_read_handle, delimiter="\t", lineterminator="\n")
            write_per_read_header(per_read_writer)

        totals: dict[str, list[ReadDistance]] = defaultdict(list)
        locus_means: dict[str, list[float]] = defaultdict(list)
        ref_counts: dict[str, int] = defaultdict(int)
        bam_specs = [BamSpec("ONT", args.ont_bam), BamSpec("PacBio", args.pacbio_bam)]

        with (
            pysam.FastaFile(args.maternal_fa) as maternal_fa,
            pysam.FastaFile(args.paternal_fa) as paternal_fa,
            pysam.AlignmentFile(args.ont_bam, "rb") as ont_bam,
            pysam.AlignmentFile(args.pacbio_bam, "rb") as pacbio_bam,
            open(output_path, "w", encoding="utf-8", newline="") as out_handle,
        ):
            writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
            write_summary_header(writer)
            bam_by_label = {"ONT": ont_bam, "PacBio": pacbio_bam}
            distance_fn = get_levenshtein_distance()

            if args.processes == 1:
                for idx, locus in enumerate(loci, 1):
                    if args.progress_every and (idx == 1 or idx == len(loci) or idx % args.progress_every == 0):
                        print_progress(idx, len(loci), locus)
                    started = time.monotonic()
                    results = process_locus(
                        locus, args, maternal_fa, paternal_fa, bam_by_label, bam_specs, distance_fn
                    )
                    record_locus_results(writer, per_read_writer, results, totals, locus_means, ref_counts)
                    if args.progress_every and (idx == 1 or idx == len(loci) or idx % args.progress_every == 0):
                        print_progress(idx, len(loci), locus, results)
                    _ = time.monotonic() - started
            else:
                pending: dict[int, list[LocusToolResult]] = {}
                next_to_write = 0
                completed = 0
                with ProcessPoolExecutor(
                    max_workers=args.processes,
                    initializer=init_worker,
                    initargs=(args,),
                ) as executor:
                    futures = {
                        executor.submit(process_locus_worker, (idx, locus)): idx
                        for idx, locus in enumerate(loci, 1)
                    }
                    for future in as_completed(futures):
                        idx, results = future.result()
                        completed += 1
                        pending[idx] = results
                        while next_to_write + 1 in pending:
                            next_to_write += 1
                            ordered = pending.pop(next_to_write)
                            locus = loci[next_to_write - 1]
                            record_locus_results(
                                writer, per_read_writer, ordered, totals, locus_means, ref_counts
                            )
                            if args.progress_every and (
                                next_to_write == 1
                                or next_to_write == len(loci)
                                or next_to_write % args.progress_every == 0
                            ):
                                print_progress(
                                    next_to_write,
                                    len(loci),
                                    locus,
                                    ordered,
                                    completed,
                                    next_to_write,
                                )

            all_locus_means: list[float] = []
            for spec in bam_specs:
                items = totals[spec.label]
                distances = [item.distance for item in items]
                read_lengths = [item.read_length for item in items]
                ref_lengths = [item.reference_length for item in items]
                left_flanks = [item.left_flank_used for item in items]
                right_flanks = [item.right_flank_used for item in items]
                all_locus_means.extend(locus_means[spec.label])
                write_summary_row(
                    writer,
                    "tool_total",
                    None,
                    spec.label,
                    ref_counts[spec.label],
                    distances,
                    read_lengths,
                    ref_lengths,
                    left_flanks,
                    right_flanks,
                    len(locus_means[spec.label]),
                    statistics.fmean(locus_means[spec.label]) if locus_means[spec.label] else None,
                )

            all_items = totals["ONT"] + totals["PacBio"]
            write_summary_row(
                writer,
                "all_tools_total",
                None,
                "ALL",
                sum(ref_counts.values()),
                [item.distance for item in all_items],
                [item.read_length for item in all_items],
                [item.reference_length for item in all_items],
                [item.left_flank_used for item in all_items],
                [item.right_flank_used for item in all_items],
                len(all_locus_means),
                statistics.fmean(all_locus_means) if all_locus_means else None,
            )
            print(
                f"[read-reference-levenshtein-fasta] done bed={args.bed} "
                f"processed_loci={len(loci)}",
                file=sys.stderr,
                flush=True,
            )
    finally:
        if per_read_handle is not None:
            per_read_handle.close()


def iter_bed_files(bed_dir: Path, pattern: str, recursive: bool) -> list[Path]:
    globber = bed_dir.rglob if recursive else bed_dir.glob
    return sorted(path for path in globber(pattern) if path.is_file())


def build_single_bed_args(args: argparse.Namespace, bed_path: Path) -> argparse.Namespace:
    stem = bed_path.name[:-4] if bed_path.name.endswith(".bed") else bed_path.stem
    per_read_tsv = None
    if args.write_per_read:
        per_read_tsv = str(args.output_dir / f"{stem}{args.per_read_suffix}")
    return argparse.Namespace(
        bed=str(bed_path),
        maternal_fa=args.maternal_fa,
        paternal_fa=args.paternal_fa,
        ont_bam=args.ont_bam,
        pacbio_bam=args.pacbio_bam,
        output_tsv=str(args.output_dir / f"{stem}{args.summary_suffix}"),
        per_read_tsv=per_read_tsv,
        flank=args.flank,
        fetch_pad=args.fetch_pad,
        max_reads_per_locus=args.max_reads_per_locus,
        min_mapq=args.min_mapq,
        primary_only=args.primary_only,
        include_secondary=args.include_secondary,
        chrom=args.chrom,
        progress_every=args.progress_every,
        processes=args.processes,
    )


def main() -> int:
    args = build_parser().parse_args()
    if not args.bed_dir.is_dir():
        raise SystemExit(f"BED directory does not exist: {args.bed_dir}")
    bed_files = iter_bed_files(args.bed_dir, args.bed_pattern, args.recursive)
    if not bed_files:
        raise SystemExit(f"No BED files matched {args.bed_pattern!r} in {args.bed_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    failures: list[tuple[Path, BaseException]] = []

    for index, bed_path in enumerate(bed_files, 1):
        single_args = build_single_bed_args(args, bed_path)
        output_path = Path(single_args.output_tsv)
        if args.skip_existing and output_path.exists():
            print(
                f"[batch-read-reference-levenshtein-fasta] SKIP {index}/{len(bed_files)} {bed_path}",
                file=sys.stderr,
                flush=True,
            )
            continue
        print(
            f"[batch-read-reference-levenshtein-fasta] START {index}/{len(bed_files)} {bed_path}",
            file=sys.stderr,
            flush=True,
        )
        try:
            run_one_bed(single_args)
        except (Exception, SystemExit) as exc:
            failures.append((bed_path, exc))
            print(
                f"[batch-read-reference-levenshtein-fasta] ERROR {bed_path}: {exc}",
                file=sys.stderr,
                flush=True,
            )
            if not args.continue_on_error:
                raise

    if failures:
        print(
            f"[batch-read-reference-levenshtein-fasta] failed_beds={len(failures)}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
