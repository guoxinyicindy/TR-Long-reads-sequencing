#!/usr/bin/env python3
"""Lift TRF BED intervals from haplotype assembly coordinates to hg38.

The input BED files are expected to use assembly contig names such as
chr1_MATERNAL and chr1_PATERNAL. The script uses maternal/paternal
assembly-to-hg38 BAM/CRAM alignments to project those query intervals onto
hg38 reference coordinates.

Output BED files use the hg38 reference contig name in column 1, e.g. chr1.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


HAP_SUFFIX_RE = re.compile(r"_(MATERNAL|PATERNAL)$")
HAP_IN_NAME_RE = re.compile(r"(?P<chrom>.+)_(?P<hap>MATERNAL|PATERNAL)(?:\.|_|$)")

CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQ = 7
CIGAR_X = 8

ALIGNED_QUERY_AND_REF = {CIGAR_M, CIGAR_EQ, CIGAR_X}
REF_GAP = {CIGAR_D, CIGAR_N}
QUERY_ONLY = {CIGAR_I, CIGAR_S}
QUERY_COORD_CONSUMING = {CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_H, CIGAR_EQ, CIGAR_X}


@dataclass(frozen=True)
class BedRecord:
    raw_chrom: str
    start: int
    end: int
    fields: tuple[str, ...]
    line_no: int
    source_haplotype: str = "unknown"

    @property
    def haplotype(self) -> str:
        match = HAP_SUFFIX_RE.search(self.raw_chrom)
        if match:
            return match.group(1).lower()
        return self.source_haplotype

    @property
    def normalized_chrom(self) -> str:
        return normalize_haplotype_chrom(self.raw_chrom)

    @property
    def locus_id(self) -> str:
        return f"{self.raw_chrom}:{self.start}-{self.end}"


@dataclass(frozen=True)
class LiftedRecord:
    record: BedRecord
    chrom: str
    start: int
    end: int
    alignment_query: str
    strand: str
    mapq: int
    aligned_query_bases: int
    query_part_start: int
    query_part_end: int
    liftover_kind: str


@dataclass(frozen=True)
class FailedRecord:
    record: BedRecord
    reason: str


@dataclass(frozen=True)
class QueryAlignment:
    query_name: str
    reference_name: str
    reference_start: int
    reference_end: int
    mapping_quality: int
    is_reverse: bool
    is_supplementary: bool
    cigartuples: tuple[tuple[int, int], ...]
    query_length: int
    query_span_start: int
    query_span_end: int
    alignment_score: int


@dataclass(frozen=True)
class LiftoverJob:
    bed_path: Path
    output_bed_path: Path
    failed_path: Path


_WORKER_MATERNAL_INDEX: dict[str, list[QueryAlignment]] | None = None
_WORKER_PATERNAL_INDEX: dict[str, list[QueryAlignment]] | None = None
_WORKER_ARGS: argparse.Namespace | None = None


def normalize_haplotype_chrom(chrom: str) -> str:
    return HAP_SUFFIX_RE.sub("", chrom)


def parse_haplotype_from_filename(path: Path) -> tuple[str | None, str]:
    match = HAP_IN_NAME_RE.search(path.name)
    if not match:
        return None, "unknown"
    return match.group("chrom"), match.group("hap").lower()


def iter_bed_files(input_dir: Path, pattern: str, recursive: bool) -> list[Path]:
    globber = input_dir.rglob if recursive else input_dir.glob
    return sorted(path for path in globber(pattern) if path.is_file())


def output_stem_for_bed(bed_path: Path, input_dir: Path, recursive: bool) -> str:
    if not recursive:
        return bed_path.name[:-4] if bed_path.name.endswith(".bed") else bed_path.stem
    rel = bed_path.relative_to(input_dir)
    name = str(rel)
    if name.endswith(".bed"):
        name = name[:-4]
    return name.replace("/", "__")


def read_bed(path: Path) -> list[BedRecord]:
    records: list[BedRecord] = []
    filename_chrom, filename_haplotype = parse_haplotype_from_filename(path)
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = tuple(line.split())
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no} has fewer than 3 BED columns")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no} has non-integer start/end") from exc
            if end <= start:
                raise ValueError(f"{path}:{line_no} has end <= start")
            records.append(
                BedRecord(
                    raw_chrom=fields[0],
                    start=start,
                    end=end,
                    fields=fields,
                    line_no=line_no,
                    source_haplotype=filename_haplotype,
                )
            )
    return records


def candidate_score(aln, aligned_query_bases: int) -> tuple[int, int, int, int, int]:
    aln_score = int(getattr(aln, "alignment_score", 0))
    primary = 0 if aln.is_supplementary else 1
    span = (aln.reference_end or 0) - (aln.reference_start or 0)
    return primary, int(aln.mapping_quality), aln_score, aligned_query_bases, span


def should_skip_alignment(aln, args: argparse.Namespace) -> bool:
    if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail:
        return True
    if aln.is_secondary and not args.include_secondary:
        return True
    if aln.is_supplementary and args.primary_only:
        return True
    if aln.mapping_quality < args.min_mapq:
        return True
    return False


def build_query_alignment_index(aln_path: str, args: argparse.Namespace) -> dict[str, list[QueryAlignment]]:
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit("Missing dependency: pysam. Activate/install the project environment first.") from exc

    index: dict[str, list[QueryAlignment]] = defaultdict(list)
    with pysam.AlignmentFile(aln_path, "rb") as aln_file:
        for aln in aln_file.fetch(until_eof=True):
            if should_skip_alignment(aln, args):
                continue
            if aln.query_name is None:
                continue
            if aln.cigartuples is None or aln.reference_start is None:
                continue
            if aln.reference_name is None or aln.reference_end is None:
                continue
            query_length = get_query_length(aln)
            if query_length is None:
                continue
            cigartuples = tuple((int(op), int(length)) for op, length in aln.cigartuples)
            query_span = alignment_query_span_from_parts(
                cigartuples,
                bool(aln.is_reverse),
                query_length,
            )
            if query_span is None:
                continue
            try:
                alignment_score = int(aln.get_tag("AS"))
            except KeyError:
                alignment_score = 0
            index[aln.query_name].append(
                QueryAlignment(
                    query_name=aln.query_name,
                    reference_name=aln.reference_name,
                    reference_start=int(aln.reference_start),
                    reference_end=int(aln.reference_end),
                    mapping_quality=int(aln.mapping_quality),
                    is_reverse=bool(aln.is_reverse),
                    is_supplementary=bool(aln.is_supplementary),
                    cigartuples=cigartuples,
                    query_length=query_length,
                    query_span_start=query_span[0],
                    query_span_end=query_span[1],
                    alignment_score=alignment_score,
                )
            )

    for alignments in index.values():
        alignments.sort(
            key=lambda aln: (
                *(alignment_query_span(aln) or (sys.maxsize, sys.maxsize)),
                aln.reference_name,
                aln.reference_start,
            )
        )
    return index


def candidate_query_names(record: BedRecord) -> tuple[str, ...]:
    names = [record.raw_chrom, record.normalized_chrom]
    if record.haplotype in {"maternal", "paternal"}:
        names.append(f"{record.normalized_chrom}_{record.haplotype.upper()}")
    seen = set()
    unique = []
    for name in names:
        if name not in seen:
            seen.add(name)
            unique.append(name)
    return tuple(unique)


def get_query_length(aln) -> int | None:
    if getattr(aln, "cigartuples", None) is not None:
        query_length_from_cigar = sum(
            int(length)
            for op, length in aln.cigartuples
            if op in QUERY_COORD_CONSUMING
        )
        if query_length_from_cigar > 0:
            return query_length_from_cigar
    if hasattr(aln, "infer_query_length"):
        try:
            inferred_length = aln.infer_query_length(always=True)
            if inferred_length is not None:
                return int(inferred_length)
        except (TypeError, ValueError):
            pass
    query_length = getattr(aln, "query_length", None)
    if query_length is not None:
        return int(query_length)
    return None


def query_block_for_cigar_op_parts(
    is_reverse: bool,
    query_length: int,
    query_cursor: int,
    length: int,
) -> tuple[int, int]:
    if not is_reverse:
        return query_cursor, query_cursor + length
    return query_length - query_cursor - length, query_length - query_cursor


def query_block_for_cigar_op(aln, query_cursor: int, length: int) -> tuple[int, int] | None:
    query_length = get_query_length(aln)
    if query_length is None:
        return None
    return query_block_for_cigar_op_parts(aln.is_reverse, query_length, query_cursor, length)


def alignment_query_span_from_parts(
    cigartuples: tuple[tuple[int, int], ...],
    is_reverse: bool,
    query_length: int,
) -> tuple[int, int] | None:
    query_cursor = 0
    ranges: list[tuple[int, int]] = []
    for op, length in cigartuples:
        if op in ALIGNED_QUERY_AND_REF or op == CIGAR_I:
            ranges.append(
                query_block_for_cigar_op_parts(
                    is_reverse,
                    query_length,
                    query_cursor,
                    length,
                )
            )
            query_cursor += length
            continue
        if op in {CIGAR_S, CIGAR_H}:
            query_cursor += length
            continue
        if op in REF_GAP or op == CIGAR_P:
            continue
    if not ranges:
        return None
    return min(start for start, _ in ranges), max(end for _, end in ranges)


def alignment_query_span(aln) -> tuple[int, int] | None:
    if hasattr(aln, "query_span_start") and hasattr(aln, "query_span_end"):
        return int(aln.query_span_start), int(aln.query_span_end)
    if aln.cigartuples is None:
        return None
    query_length = get_query_length(aln)
    if query_length is None:
        return None
    return alignment_query_span_from_parts(
        tuple((int(op), int(length)) for op, length in aln.cigartuples),
        bool(aln.is_reverse),
        query_length,
    )


def add_ref_projection_for_query_overlap(
    aln,
    query_cursor: int,
    ref_cursor: int,
    length: int,
    start: int,
    end: int,
    ref_ranges: list[tuple[int, int]],
) -> int:
    query_block = query_block_for_cigar_op(aln, query_cursor, length)
    if query_block is None:
        return 0
    query_block_start, query_block_end = query_block
    overlap_start = max(start, query_block_start)
    overlap_end = min(end, query_block_end)
    if overlap_start >= overlap_end:
        return 0

    if aln.is_reverse:
        ref_start = ref_cursor + (query_block_end - overlap_end)
        ref_end = ref_cursor + (query_block_end - overlap_start)
    else:
        ref_start = ref_cursor + (overlap_start - query_block_start)
        ref_end = ref_cursor + (overlap_end - query_block_start)
    ref_ranges.append((ref_start, ref_end))
    return overlap_end - overlap_start


def map_query_interval(
    aln,
    start: int,
    end: int,
    require_full_coverage: bool,
) -> tuple[str, int, int, int] | None:
    query_span = alignment_query_span(aln)
    if query_span is None:
        return None
    query_alignment_start, query_alignment_end = query_span
    if require_full_coverage and (
        start < query_alignment_start or end > query_alignment_end
    ):
        return None
    if end <= query_alignment_start or start >= query_alignment_end:
        return None

    if aln.cigartuples is None or aln.reference_start is None:
        return None

    query_cursor = 0
    ref_cursor = aln.reference_start
    ref_ranges: list[tuple[int, int]] = []
    aligned_query_bases = 0

    for op, length in aln.cigartuples:
        if op in ALIGNED_QUERY_AND_REF:
            aligned_query_bases += add_ref_projection_for_query_overlap(
                aln,
                query_cursor,
                ref_cursor,
                length,
                start,
                end,
                ref_ranges,
            )
            query_cursor += length
            ref_cursor += length
            continue

        if op in QUERY_ONLY:
            query_cursor += length
            continue

        if op in REF_GAP:
            ref_cursor += length
            continue

        if op == CIGAR_H:
            query_cursor += length
            continue

        if op == CIGAR_P:
            continue

    if not ref_ranges:
        return None

    ref_chrom = normalize_haplotype_chrom(aln.reference_name)
    ref_start = min(start for start, _ in ref_ranges)
    ref_end = max(end for _, end in ref_ranges)
    return ref_chrom, ref_start, ref_end, aligned_query_bases


def choose_alignment(
    record: BedRecord,
    alignments_by_query: dict[str, list[QueryAlignment]],
    require_full_coverage: bool,
) -> tuple[QueryAlignment, str, int, int, int] | None:
    candidates = []
    for query_name in candidate_query_names(record):
        for aln in alignments_by_query.get(query_name, []):
            mapped = map_query_interval(
                aln,
                record.start,
                record.end,
                require_full_coverage=require_full_coverage,
            )
            if mapped is None:
                continue
            ref_chrom, ref_start, ref_end, aligned_query_bases = mapped
            candidates.append(
                (
                    candidate_score(aln, aligned_query_bases),
                    aln,
                    ref_chrom,
                    ref_start,
                    ref_end,
                    aligned_query_bases,
                )
            )
    if not candidates:
        return None
    _, aln, ref_chrom, ref_start, ref_end, aligned_query_bases = max(
        candidates, key=lambda item: item[0]
    )
    return aln, ref_chrom, ref_start, ref_end, aligned_query_bases


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    merged: list[tuple[int, int]] = []
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
            continue
        merged.append((current_start, current_end))
        current_start, current_end = start, end
    merged.append((current_start, current_end))
    return merged


def intervals_cover_span(
    intervals: list[tuple[int, int]],
    start: int,
    end: int,
) -> bool:
    cursor = start
    for interval_start, interval_end in merge_intervals(intervals):
        if interval_end <= cursor:
            continue
        if interval_start > cursor:
            return False
        cursor = max(cursor, interval_end)
        if cursor >= end:
            return True
    return cursor >= end


def record_alignments(
    record: BedRecord,
    indexes: list[dict[str, list[QueryAlignment]]],
) -> list[QueryAlignment]:
    alignments: list[QueryAlignment] = []
    seen = set()
    for index in indexes:
        for query_name in candidate_query_names(record):
            for aln in index.get(query_name, []):
                key = (
                    aln.query_name,
                    aln.reference_name,
                    aln.reference_start,
                    aln.reference_end,
                    aln.is_reverse,
                    aln.is_supplementary,
                    aln.cigartuples,
                )
                if key in seen:
                    continue
                seen.add(key)
                alignments.append(aln)
    return alignments


def split_liftover_record(
    record: BedRecord,
    indexes: list[dict[str, list[QueryAlignment]]],
    args: argparse.Namespace,
) -> list[LiftedRecord] | FailedRecord:
    alignments = []
    covered_query_parts = []
    for aln in record_alignments(record, indexes):
        query_span = alignment_query_span(aln)
        if query_span is None:
            continue
        query_start, query_end = query_span
        overlap_start = max(record.start, query_start)
        overlap_end = min(record.end, query_end)
        if overlap_start >= overlap_end:
            continue
        alignments.append((aln, overlap_start, overlap_end))
        covered_query_parts.append((overlap_start, overlap_end))

    if not alignments:
        mode = "full" if not args.allow_partial else "partial"
        return FailedRecord(record, f"no_{mode}_alignment_for_query_interval")

    if not args.allow_partial and not intervals_cover_span(
        covered_query_parts,
        record.start,
        record.end,
    ):
        return FailedRecord(record, "no_full_alignment_for_query_interval")

    lifted_parts = []
    for aln, overlap_start, overlap_end in alignments:
        mapped = map_query_interval(
            aln,
            record.start,
            record.end,
            require_full_coverage=False,
        )
        if mapped is None:
            continue
        ref_chrom, ref_start, ref_end, aligned_query_bases = mapped
        lifted_parts.append(
            LiftedRecord(
                record=record,
                chrom=ref_chrom,
                start=ref_start,
                end=ref_end,
                alignment_query=aln.query_name,
                strand="-" if aln.is_reverse else "+",
                mapq=int(aln.mapping_quality),
                aligned_query_bases=aligned_query_bases,
                query_part_start=overlap_start,
                query_part_end=overlap_end,
                liftover_kind="split",
            )
        )

    if not lifted_parts:
        return FailedRecord(record, "no_reference_projection_for_query_interval")

    lifted_parts.sort(
        key=lambda item: (
            item.query_part_start,
            item.query_part_end,
            item.chrom,
            item.start,
            item.end,
        )
    )
    return lifted_parts


def liftover_record(
    record: BedRecord,
    maternal_index: dict[str, list[QueryAlignment]],
    paternal_index: dict[str, list[QueryAlignment]],
    args: argparse.Namespace,
) -> list[LiftedRecord] | FailedRecord:
    if record.haplotype == "maternal":
        indexes = [maternal_index]
    elif record.haplotype == "paternal":
        indexes = [paternal_index]
    else:
        indexes = [maternal_index, paternal_index]

    candidates = []
    for index in indexes:
        chosen = choose_alignment(
            record,
            index,
            require_full_coverage=not args.allow_partial,
        )
        if chosen is not None:
            aln, ref_chrom, ref_start, ref_end, aligned_query_bases = chosen
            candidates.append(
                (
                    candidate_score(aln, aligned_query_bases),
                    aln,
                    ref_chrom,
                    ref_start,
                    ref_end,
                    aligned_query_bases,
                )
            )

    if not candidates:
        if args.allow_split_liftover:
            return split_liftover_record(record, indexes, args)
        mode = "full" if not args.allow_partial else "partial"
        return FailedRecord(record, f"no_{mode}_alignment_for_query_interval")

    _, aln, ref_chrom, ref_start, ref_end, aligned_query_bases = max(
        candidates, key=lambda item: item[0]
    )
    return [
        LiftedRecord(
            record=record,
            chrom=ref_chrom,
            start=ref_start,
            end=ref_end,
            alignment_query=aln.query_name,
            strand="-" if aln.is_reverse else "+",
            mapq=int(aln.mapping_quality),
            aligned_query_bases=aligned_query_bases,
            query_part_start=record.start,
            query_part_end=record.end,
            liftover_kind="single",
        )
    ]


def lifted_fields(item: LiftedRecord) -> list[str | int]:
    record = item.record
    return [
        item.chrom,
        item.start,
        item.end,
        record.locus_id,
        record.raw_chrom,
        record.start,
        record.end,
        record.haplotype.upper(),
        item.alignment_query,
        item.strand,
        item.mapq,
        item.aligned_query_bases,
        *record.fields[3:],
    ]


def write_lifted_bed(path: Path, records: Iterable[LiftedRecord]) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for record in records:
            writer.writerow(lifted_fields(record))
            count += 1
    return count


def write_failed_tsv(path: Path, records: Iterable[FailedRecord]) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "input_chrom",
                "input_start",
                "input_end",
                "input_line_no",
                "haplotype",
                "reason",
            ]
        )
        for item in records:
            writer.writerow(
                [
                    item.record.raw_chrom,
                    item.record.start,
                    item.record.end,
                    item.record.line_no,
                    item.record.haplotype,
                    item.reason,
                ]
            )
            count += 1
    return count


def process_bed_file(
    bed_path: Path,
    output_bed_path: Path,
    failed_path: Path,
    maternal_index: dict[str, list[QueryAlignment]],
    paternal_index: dict[str, list[QueryAlignment]],
    args: argparse.Namespace,
) -> tuple[int, int]:
    lifted: list[LiftedRecord] = []
    failed: list[FailedRecord] = []
    for record in read_bed(bed_path):
        result = liftover_record(record, maternal_index, paternal_index, args)
        if isinstance(result, FailedRecord):
            failed.append(result)
        else:
            lifted.extend(result)

    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    failed_path.parent.mkdir(parents=True, exist_ok=True)
    lifted_count = write_lifted_bed(output_bed_path, lifted)
    failed_count = write_failed_tsv(failed_path, failed)
    return lifted_count, failed_count


def init_worker(
    maternal_index: dict[str, list[QueryAlignment]],
    paternal_index: dict[str, list[QueryAlignment]],
    args: argparse.Namespace,
) -> None:
    global _WORKER_MATERNAL_INDEX, _WORKER_PATERNAL_INDEX, _WORKER_ARGS
    _WORKER_MATERNAL_INDEX = maternal_index
    _WORKER_PATERNAL_INDEX = paternal_index
    _WORKER_ARGS = args


def process_bed_file_worker(job: LiftoverJob) -> tuple[Path, int, int]:
    if _WORKER_MATERNAL_INDEX is None or _WORKER_PATERNAL_INDEX is None or _WORKER_ARGS is None:
        raise RuntimeError("Worker state was not initialized")
    lifted_count, failed_count = process_bed_file(
        job.bed_path,
        job.output_bed_path,
        job.failed_path,
        _WORKER_MATERNAL_INDEX,
        _WORKER_PATERNAL_INDEX,
        _WORKER_ARGS,
    )
    return job.bed_path, lifted_count, failed_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Lift TRF BED intervals from haplotype assembly coordinates to hg38 BED coordinates."
    )
    parser.add_argument("--bed-dir", required=True, type=Path, help="Directory containing TRF BED files.")
    parser.add_argument(
        "--bed-pattern",
        default="*.bed",
        help="Glob pattern for BED files inside --bed-dir. Default: *.bed",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search for BED files recursively under --bed-dir.",
    )
    parser.add_argument("--maternal-aln", required=True, help="Maternal assembly-to-hg38 BAM/CRAM.")
    parser.add_argument("--paternal-aln", required=True, help="Paternal assembly-to-hg38 BAM/CRAM.")
    parser.add_argument("-o", "--output-dir", required=True, type=Path, help="Directory for lifted BED files.")
    parser.add_argument(
        "--output-suffix",
        default=".hg38.bed",
        help="Suffix appended to each input BED basename. Default: .hg38.bed",
    )
    parser.add_argument(
        "--failed-suffix",
        default=".unlifted.tsv",
        help="Suffix appended to each input BED basename for failed records. Default: .unlifted.tsv",
    )
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Allow partially aligned TRF intervals to be projected using their aligned query bases.",
    )
    parser.add_argument(
        "--allow-split-liftover",
        action="store_true",
        help=(
            "If no single alignment can liftover an interval, allow multiple "
            "alignments to cover it and write one hg38 BED row per projected "
            "piece. Query insertion bases inside covered alignments are allowed "
            "but only M/=/X bases contribute hg38 coordinates."
        ),
    )
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum assembly alignment MAPQ. Default: 0.")
    parser.add_argument(
        "--include-secondary",
        action="store_true",
        help="Include secondary assembly alignments. Default: skip secondary alignments.",
    )
    parser.add_argument(
        "--primary-only",
        action="store_true",
        help="Skip supplementary assembly alignments in addition to secondary alignments.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip an input BED if its lifted BED already exists.",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of BED files to process in parallel. Default: 1.",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.bed_dir.is_dir():
        raise SystemExit(f"BED directory does not exist: {args.bed_dir}")
    if args.processes < 1:
        raise SystemExit("--processes must be >= 1")
    for label, path in (
        ("maternal alignment", args.maternal_aln),
        ("paternal alignment", args.paternal_aln),
    ):
        if not Path(path).exists():
            raise SystemExit(f"{label} does not exist: {path}")


def main() -> int:
    args = build_parser().parse_args()
    validate_args(args)

    bed_files = iter_bed_files(args.bed_dir, args.bed_pattern, args.recursive)
    if not bed_files:
        raise SystemExit(f"No BED files matched {args.bed_pattern!r} in {args.bed_dir}")

    print("[liftover-trf-bed] indexing maternal assembly alignment", file=sys.stderr, flush=True)
    maternal_index = build_query_alignment_index(args.maternal_aln, args)
    print(
        f"[liftover-trf-bed] maternal query contigs indexed={len(maternal_index)}",
        file=sys.stderr,
        flush=True,
    )
    print("[liftover-trf-bed] indexing paternal assembly alignment", file=sys.stderr, flush=True)
    paternal_index = build_query_alignment_index(args.paternal_aln, args)
    print(
        f"[liftover-trf-bed] paternal query contigs indexed={len(paternal_index)}",
        file=sys.stderr,
        flush=True,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    total_lifted = 0
    total_failed = 0
    skipped = 0
    jobs: list[LiftoverJob] = []

    for idx, bed_path in enumerate(bed_files, 1):
        stem = output_stem_for_bed(bed_path, args.bed_dir, args.recursive)
        output_bed_path = args.output_dir / f"{stem}{args.output_suffix}"
        failed_path = args.output_dir / f"{stem}{args.failed_suffix}"
        if args.skip_existing and output_bed_path.exists():
            skipped += 1
            print(
                f"[liftover-trf-bed] SKIP {idx}/{len(bed_files)} {bed_path} output_exists={output_bed_path}",
                file=sys.stderr,
                flush=True,
            )
            continue
        jobs.append(LiftoverJob(bed_path, output_bed_path, failed_path))

    if args.processes == 1:
        for done_idx, job in enumerate(jobs, 1):
            lifted_count, failed_count = process_bed_file(
                job.bed_path,
                job.output_bed_path,
                job.failed_path,
                maternal_index,
                paternal_index,
                args,
            )
            total_lifted += lifted_count
            total_failed += failed_count
            print(
                (
                    f"[liftover-trf-bed] DONE {done_idx}/{len(jobs)} {job.bed_path} "
                    f"lifted={lifted_count} unlifted={failed_count} output={job.output_bed_path}"
                ),
                file=sys.stderr,
                flush=True,
            )
    else:
        job_by_path = {job.bed_path: job for job in jobs}
        with ProcessPoolExecutor(
            max_workers=args.processes,
            initializer=init_worker,
            initargs=(maternal_index, paternal_index, args),
        ) as executor:
            future_map = {
                executor.submit(process_bed_file_worker, job): job
                for job in jobs
            }
            for done_idx, future in enumerate(as_completed(future_map), 1):
                bed_path, lifted_count, failed_count = future.result()
                job = job_by_path[bed_path]
                total_lifted += lifted_count
                total_failed += failed_count
                print(
                    (
                        f"[liftover-trf-bed] DONE {done_idx}/{len(jobs)} {bed_path} "
                        f"lifted={lifted_count} unlifted={failed_count} output={job.output_bed_path}"
                    ),
                    file=sys.stderr,
                    flush=True,
                )

    print(
        (
            "[liftover-trf-bed] finished "
            f"bed_files={len(bed_files)} skipped={skipped} "
            f"lifted_records={total_lifted} unlifted_records={total_failed}"
        ),
        file=sys.stderr,
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
