#!/usr/bin/env python3
"""
Decompose read-to-reference Levenshtein distance into insertions/deletions/substitutions.

This script uses the same core filtering and sequence extraction rules as
read_reference_levenshtein_ont_pacbio.py, but writes one row per kept read with
edit operation counts against the best-matching maternal/paternal allele.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
import time
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
TOOLS = ["ONT", "PacBio"]


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
class EditCounts:
    distance: int
    n_ins: int
    n_del: int
    n_sub: int
    n_match: int


@dataclass(frozen=True)
class ReadEdit:
    locus: Locus
    tool: str
    read_name: str
    matched_haplotype: str
    read_length: int
    reference_length: int
    left_flank_used: int
    right_flank_used: int
    edits: EditCounts


def parse_bed(path: str) -> list[Locus]:
    loci: list[Locus] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise SystemExit(f"{path}:{line_no} has fewer than 3 BED columns")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                raise SystemExit(f"{path}:{line_no} has end <= start")
            name = fields[3] if len(fields) >= 4 and fields[3] else f"{chrom}:{start}-{end}"
            loci.append(Locus(chrom=chrom, start=start, end=end, name=name))
    return loci


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


def collect_query_sequence(aln, start: int, end: int) -> str | None:
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


def extract_read_sequence(aln, locus: Locus, left_flank: int, right_flank: int) -> str | None:
    start = max(0, locus.start - left_flank)
    end = locus.end + right_flank
    return collect_query_sequence(aln, start, end)


def decompose_levenshtein(query: str, reference: str) -> EditCounts:
    """Return one minimal edit decomposition from query to reference.

    n_ins is bases present in query but not reference. n_del is bases missing
    from query relative to reference.
    """
    n = len(query)
    m = len(reference)
    previous = list(range(m + 1))
    trace: list[bytearray] = []

    for i, qbase in enumerate(query, start=1):
        current = [i]
        trace_row = bytearray(m + 1)
        trace_row[0] = ord("I")
        for j, rbase in enumerate(reference, start=1):
            sub_cost = 0 if qbase == rbase else 1
            diag = previous[j - 1] + sub_cost
            delete_from_query = previous[j] + 1
            insert_to_query = current[j - 1] + 1
            best = min(diag, delete_from_query, insert_to_query)
            current.append(best)
            if best == diag:
                trace_row[j] = ord("M") if sub_cost == 0 else ord("S")
            elif best == delete_from_query:
                trace_row[j] = ord("I")
            else:
                trace_row[j] = ord("D")
        previous = current
        trace.append(trace_row)

    i = n
    j = m
    n_ins = 0
    n_del = 0
    n_sub = 0
    n_match = 0
    while i > 0 or j > 0:
        if i == 0:
            n_del += j
            break
        if j == 0:
            n_ins += i
            break
        op = chr(trace[i - 1][j])
        if op == "M":
            n_match += 1
            i -= 1
            j -= 1
        elif op == "S":
            n_sub += 1
            i -= 1
            j -= 1
        elif op == "I":
            n_ins += 1
            i -= 1
        elif op == "D":
            n_del += 1
            j -= 1
        else:
            raise RuntimeError(f"unexpected traceback op: {op!r}")

    return EditCounts(
        distance=n_ins + n_del + n_sub,
        n_ins=n_ins,
        n_del=n_del,
        n_sub=n_sub,
        n_match=n_match,
    )


def iter_edits_for_locus(
    bam,
    spec: BamSpec,
    locus: Locus,
    ref_alleles: dict[str, str],
    flank: int,
    args: argparse.Namespace,
) -> list[ReadEdit]:
    if locus.chrom not in bam.references:
        return []

    fetch_start = max(0, locus.start - flank - args.fetch_pad)
    fetch_end = max(locus.end + flank + args.fetch_pad, locus.start + flank + args.fetch_pad + 1)
    best_by_read: dict[str, ReadEdit] = {}

    for aln in bam.fetch(locus.chrom, fetch_start, fetch_end):
        if should_skip_read(aln, args):
            continue
        seq = extract_read_sequence(aln, locus, flank, flank)
        if not seq:
            continue

        best_hap = None
        best_ref_seq = None
        best_edits = None
        for hap, ref_seq in ref_alleles.items():
            edits = decompose_levenshtein(seq, ref_seq)
            if best_edits is None or edits.distance < best_edits.distance:
                best_hap = hap
                best_ref_seq = ref_seq
                best_edits = edits
        if best_edits is None or best_hap is None or best_ref_seq is None:
            continue
        if abs(len(seq) - len(best_ref_seq)) > args.max_length_difference:
            continue

        read_name = get_original_read_name(aln)
        item = ReadEdit(
            locus=locus,
            tool=spec.label,
            read_name=read_name,
            matched_haplotype=best_hap,
            read_length=len(seq),
            reference_length=len(best_ref_seq),
            left_flank_used=flank,
            right_flank_used=flank,
            edits=best_edits,
        )
        previous = best_by_read.get(read_name)
        if previous is None or item.edits.distance < previous.edits.distance:
            best_by_read[read_name] = item

    return list(best_by_read.values())


def fmt(value: float | int | None) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{value:.6f}"


def mean_or_none(values: list[int]) -> float | None:
    return statistics.fmean(values) if values else None


def write_per_read_header(writer: csv.writer) -> None:
    writer.writerow(
        [
            "locus_id",
            "chrom",
            "start",
            "end",
            "tool",
            "read_name",
            "matched_haplotype",
            "distance",
            "n_ins",
            "n_del",
            "n_sub",
            "n_match",
            "read_length",
            "reference_length",
            "length_diff",
            "abs_length_diff",
            "ins_per_ref_bp",
            "del_per_ref_bp",
            "sub_per_ref_bp",
            "distance_per_ref_bp",
            "left_flank_used",
            "right_flank_used",
        ]
    )


def write_summary_header(writer: csv.writer) -> None:
    writer.writerow(
        [
            "scope",
            "locus_id",
            "chrom",
            "start",
            "end",
            "tool",
            "n_reads",
            "mean_distance",
            "mean_n_ins",
            "mean_n_del",
            "mean_n_sub",
            "mean_read_length",
            "mean_reference_length",
            "mean_length_diff",
            "mean_abs_length_diff",
            "mean_ins_per_ref_bp",
            "mean_del_per_ref_bp",
            "mean_sub_per_ref_bp",
            "mean_distance_per_ref_bp",
        ]
    )


def write_per_read_row(writer: csv.writer, item: ReadEdit) -> None:
    length_diff = item.read_length - item.reference_length
    ref_len = item.reference_length
    writer.writerow(
        [
            item.locus.name,
            item.locus.chrom,
            item.locus.start,
            item.locus.end,
            item.tool,
            item.read_name,
            item.matched_haplotype,
            item.edits.distance,
            item.edits.n_ins,
            item.edits.n_del,
            item.edits.n_sub,
            item.edits.n_match,
            item.read_length,
            item.reference_length,
            length_diff,
            abs(length_diff),
            fmt(item.edits.n_ins / ref_len if ref_len else None),
            fmt(item.edits.n_del / ref_len if ref_len else None),
            fmt(item.edits.n_sub / ref_len if ref_len else None),
            fmt(item.edits.distance / ref_len if ref_len else None),
            item.left_flank_used,
            item.right_flank_used,
        ]
    )


def write_summary_row(
    writer: csv.writer,
    scope: str,
    locus: Locus | None,
    tool: str,
    items: list[ReadEdit],
) -> None:
    distances = [item.edits.distance for item in items]
    n_ins = [item.edits.n_ins for item in items]
    n_del = [item.edits.n_del for item in items]
    n_sub = [item.edits.n_sub for item in items]
    read_lengths = [item.read_length for item in items]
    ref_lengths = [item.reference_length for item in items]
    length_diffs = [item.read_length - item.reference_length for item in items]
    abs_length_diffs = [abs(value) for value in length_diffs]
    ins_rates = [item.edits.n_ins / item.reference_length for item in items if item.reference_length]
    del_rates = [item.edits.n_del / item.reference_length for item in items if item.reference_length]
    sub_rates = [item.edits.n_sub / item.reference_length for item in items if item.reference_length]
    distance_rates = [
        item.edits.distance / item.reference_length for item in items if item.reference_length
    ]
    writer.writerow(
        [
            scope,
            locus.name if locus else "ALL",
            locus.chrom if locus else "ALL",
            locus.start if locus else "NA",
            locus.end if locus else "NA",
            tool,
            len(items),
            fmt(mean_or_none(distances)),
            fmt(mean_or_none(n_ins)),
            fmt(mean_or_none(n_del)),
            fmt(mean_or_none(n_sub)),
            fmt(mean_or_none(read_lengths)),
            fmt(mean_or_none(ref_lengths)),
            fmt(mean_or_none(length_diffs)),
            fmt(mean_or_none(abs_length_diffs)),
            fmt(statistics.fmean(ins_rates) if ins_rates else None),
            fmt(statistics.fmean(del_rates) if del_rates else None),
            fmt(statistics.fmean(sub_rates) if sub_rates else None),
            fmt(statistics.fmean(distance_rates) if distance_rates else None),
        ]
    )


def process_locus(
    locus: Locus,
    mat_aln,
    pat_aln,
    bam_by_label: dict[str, object],
    bam_specs: list[BamSpec],
    args: argparse.Namespace,
) -> dict[str, list[ReadEdit]]:
    if locus.chrom not in mat_aln.references and locus.chrom not in pat_aln.references:
        return {spec.label: [] for spec in bam_specs}

    ref_alleles = fetch_reference_alleles(mat_aln, pat_aln, locus, args.flank, args.flank)
    results: dict[str, list[ReadEdit]] = {}
    for spec in bam_specs:
        if ref_alleles:
            results[spec.label] = iter_edits_for_locus(
                bam_by_label[spec.label],
                spec,
                locus,
                ref_alleles,
                args.flank,
                args,
            )
        else:
            results[spec.label] = []
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Decompose ONT/PacBio read-reference Levenshtein distance into ins/del/sub counts."
    )
    parser.add_argument("--bed", required=True, help="BED file with 0-based half-open hg38 loci.")
    parser.add_argument("--maternal-aln", required=True, help="Maternal assembly alignment to hg38, BAM/CRAM.")
    parser.add_argument("--paternal-aln", required=True, help="Paternal assembly alignment to hg38, BAM/CRAM.")
    parser.add_argument("--ont-bam", required=True, help="ONT reads aligned to hg38, BAM/CRAM.")
    parser.add_argument("--pacbio-bam", required=True, help="PacBio reads aligned to hg38, BAM/CRAM.")
    parser.add_argument("--out-per-read-tsv", required=True, help="Output per-read edit decomposition TSV.")
    parser.add_argument("--out-summary-tsv", required=True, help="Output summary edit decomposition TSV.")
    parser.add_argument("--flank", type=int, default=0, help="Bases of flank to include on both sides. Default: 0.")
    parser.add_argument("--fetch-pad", type=int, default=50, help="Extra reference bases for fetching reads. Default: 50.")
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum read mapping quality. Default: 0.")
    parser.add_argument(
        "--max-length-difference",
        type=int,
        default=100,
        help="Keep reads with abs(read_length - matched_reference_length) <= this value. Default: 100. Negative disables.",
    )
    parser.add_argument("--primary-only", action="store_true", help="Skip secondary and supplementary alignments.")
    parser.add_argument("--include-secondary", action="store_true", help="Include secondary alignments.")
    parser.add_argument("--chrom", action="append", help="Only process this chromosome/contig. Can be repeated.")
    parser.add_argument("--progress-every", type=int, default=100, help="Print progress every N loci. Default: 100.")
    parser.add_argument(
        "--max-locus-length",
        type=int,
        default=1000,
        help="Skip loci with BED length greater than this value. Default: 1000. Negative disables.",
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
    if args.flank < 0:
        raise SystemExit("--flank must be >= 0")
    if args.max_length_difference < 0:
        args.max_length_difference = sys.maxsize
    if args.max_locus_length < 0:
        args.max_locus_length = sys.maxsize

    loci = parse_bed(args.bed)
    if args.chrom:
        chroms = set(args.chrom)
        loci = [locus for locus in loci if locus.chrom in chroms]
        if not loci:
            raise SystemExit(f"No BED loci matched --chrom: {', '.join(args.chrom)}")
    total_loci_before_length_filter = len(loci)
    loci = [locus for locus in loci if (locus.end - locus.start) <= args.max_locus_length]
    skipped_long_loci = total_loci_before_length_filter - len(loci)
    if not loci:
        raise SystemExit(f"All loci were filtered out by --max-locus-length {args.max_locus_length}.")

    per_read_path = Path(args.out_per_read_tsv)
    summary_path = Path(args.out_summary_tsv)
    per_read_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    bam_specs = [BamSpec("ONT", args.ont_bam), BamSpec("PacBio", args.pacbio_bam)]
    totals: dict[str, list[ReadEdit]] = defaultdict(list)

    with (
        pysam.AlignmentFile(args.maternal_aln, "rb") as mat_aln,
        pysam.AlignmentFile(args.paternal_aln, "rb") as pat_aln,
        pysam.AlignmentFile(args.ont_bam, "rb") as ont_bam,
        pysam.AlignmentFile(args.pacbio_bam, "rb") as pacbio_bam,
        open(per_read_path, "w", encoding="utf-8", newline="") as per_read_handle,
        open(summary_path, "w", encoding="utf-8", newline="") as summary_handle,
    ):
        bam_by_label = {"ONT": ont_bam, "PacBio": pacbio_bam}
        per_read_writer = csv.writer(per_read_handle, delimiter="\t", lineterminator="\n")
        summary_writer = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
        write_per_read_header(per_read_writer)
        write_summary_header(summary_writer)

        for idx, locus in enumerate(loci, start=1):
            do_progress = args.progress_every and (
                idx == 1 or idx == len(loci) or idx % args.progress_every == 0
            )
            if do_progress:
                print(
                    f"[decompose-read-ref-edits] START locus {idx}/{len(loci)} "
                    f"{locus.chrom}:{locus.start}-{locus.end}",
                    file=sys.stderr,
                    flush=True,
                )
            started_at = time.monotonic()
            results = process_locus(locus, mat_aln, pat_aln, bam_by_label, bam_specs, args)
            for spec in bam_specs:
                items = results[spec.label]
                for item in items:
                    write_per_read_row(per_read_writer, item)
                write_summary_row(summary_writer, "locus", locus, spec.label, items)
                totals[spec.label].extend(items)
            if do_progress:
                elapsed = time.monotonic() - started_at
                print(
                    f"[decompose-read-ref-edits] DONE locus {idx}/{len(loci)} "
                    f"{locus.chrom}:{locus.start}-{locus.end} "
                    f"ONT={len(results.get('ONT', []))} PacBio={len(results.get('PacBio', []))} "
                    f"elapsed={elapsed:.2f}s",
                    file=sys.stderr,
                    flush=True,
                )

        all_items: list[ReadEdit] = []
        for spec in bam_specs:
            items = totals[spec.label]
            all_items.extend(items)
            write_summary_row(summary_writer, "tool_total", None, spec.label, items)
        write_summary_row(summary_writer, "all_tools_total", None, "ALL", all_items)

    print(
        (
            "[decompose-read-ref-edits] done "
            f"processed_loci={len(loci)} skipped_long_loci={skipped_long_loci} "
            f"out_per_read={per_read_path} out_summary={summary_path}"
        ),
        file=sys.stderr,
        flush=True,
    )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
