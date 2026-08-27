#!/usr/bin/env python3
"""Calculate CIGAR-based insertion/deletion/substitution rates for BED loci.

For every BED locus, this script examines ONT and PacBio alignments to the
same HG002 haplotype FASTA coordinate system.  The reference haplotype is
chosen from the BAM target contig when it has a ``_MATERNAL`` or
``_PATERNAL`` suffix; otherwise the haplotype encoded by the BED contig is
used, and an unqualified locus is compared with both FASTA files.

Insertion and deletion counts come directly from CIGAR ``I`` and ``D``
operations.  Substitutions are counted by comparing bases in CIGAR ``M``,
``=`` and ``X`` blocks with the selected haplotype FASTA sequence.  ``N``
(reference skip) is reported separately and is not counted as a deletion.

When one read has multiple eligible alignments, including primary and
supplementary alignments, the alignment with the smallest compiled
Levenshtein distance to the selected reference sequence is retained.  Ties
prefer primary, then supplementary, then secondary alignments.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import re
import statistics
import sys
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

ALIGNED_QUERY_AND_REF = {CIGAR_M, CIGAR_EQ, CIGAR_X}
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
class EditMetrics:
    selection_distance: int
    insertion_events: int
    insertion_bases: int
    deletion_events: int
    deletion_bases: int
    substitution_bases: int
    match_bases: int
    reference_skip_bases: int
    query_sequence: str
    reference_sequence: str


@dataclass(frozen=True)
class ReadMetrics:
    locus: Locus
    tool: str
    read_name: str
    reference_name: str
    matched_haplotype: str
    alignment_type: str
    mapq: int
    is_reverse: bool
    metrics: EditMetrics


_WORKER_STATE = None


def parse_bed(path: Path) -> list[Locus]:
    loci: list[Locus] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no} has fewer than 3 BED columns")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no} has invalid interval")
            name = fields[3] if len(fields) >= 4 and fields[3] else f"{chrom}:{start}-{end}"
            loci.append(Locus(chrom, start, end, name))
    return loci


def haplotype_from_chrom(chrom: str) -> str | None:
    match = HAP_SUFFIX_RE.search(chrom)
    return match.group(1).lower() if match else None


def base_chrom(chrom: str) -> str:
    return HAP_SUFFIX_RE.sub("", chrom)


def unique(items: list[str]) -> list[str]:
    return list(dict.fromkeys(items))


def bam_contigs_for_locus(bam, chrom: str) -> list[str]:
    """Return BAM contigs that can represent a BED chromosome."""
    references = set(bam.references)
    explicit_hap = haplotype_from_chrom(chrom)
    base = base_chrom(chrom)
    candidates = [chrom, base]
    if explicit_hap:
        candidates.append(f"{base}_{explicit_hap.upper()}")
    else:
        candidates.extend([f"{base}_MATERNAL", f"{base}_PATERNAL"])
    return [candidate for candidate in unique(candidates) if candidate in references]


def fasta_contig_candidates(reference_name: str, locus_chrom: str, haplotype: str) -> list[str]:
    reference_base = base_chrom(reference_name)
    locus_base = base_chrom(locus_chrom)
    return unique(
        [
            reference_name,
            reference_base,
            locus_chrom,
            locus_base,
            f"{reference_base}_{haplotype.upper()}",
            f"{locus_base}_{haplotype.upper()}",
        ]
    )


def fetch_reference(
    maternal_fa,
    paternal_fa,
    reference_name: str,
    locus_chrom: str,
    start: int,
    end: int,
) -> list[tuple[str, str]]:
    """Fetch candidate haplotype sequences for one BAM target contig."""
    alignment_hap = haplotype_from_chrom(reference_name)
    locus_hap = haplotype_from_chrom(locus_chrom)
    haplotypes = [alignment_hap or locus_hap] if (alignment_hap or locus_hap) else ["maternal", "paternal"]
    fasta_by_hap = {"maternal": maternal_fa, "paternal": paternal_fa}
    result: list[tuple[str, str]] = []
    for haplotype in haplotypes:
        fasta = fasta_by_hap[haplotype]
        contig = next(
            (candidate for candidate in fasta_contig_candidates(reference_name, locus_chrom, haplotype)
             if candidate in fasta.references),
            None,
        )
        if contig is None:
            continue
        contig_length = fasta.get_reference_length(contig)
        fetch_start = max(0, start)
        fetch_end = min(contig_length, end)
        if fetch_end <= fetch_start:
            continue
        result.append((haplotype, fasta.fetch(contig, fetch_start, fetch_end).upper()))
    return result


def spans_interval(aln, start: int, end: int) -> bool:
    return (
        aln.reference_start is not None
        and aln.reference_end is not None
        and aln.reference_start <= start
        and aln.reference_end >= end
    )


def collect_query_sequence(aln, start: int, end: int) -> str | None:
    """Extract query bases corresponding to a fully covered reference interval."""
    if aln.cigartuples is None or aln.query_sequence is None or not spans_interval(aln, start, end):
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
        elif op == CIGAR_I:
            if start <= ref_pos < end:
                pieces.append(aln.query_sequence[query_pos : query_pos + length])
            query_pos += length
        elif op in {CIGAR_D, CIGAR_N}:
            ref_pos += length
        elif op == CIGAR_S:
            query_pos += length
        elif op in {CIGAR_H, CIGAR_P}:
            continue
        else:
            query_pos += length
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
    return aln.mapping_quality < args.min_mapq


def alignment_type(aln) -> str:
    if aln.is_secondary:
        return "secondary"
    if aln.is_supplementary:
        return "supplementary"
    return "primary"


def cigar_metrics(aln, start: int, end: int, reference: str, query: str, selection_distance: int) -> EditMetrics:
    """Count I/D from CIGAR and mismatch bases in M/=/X blocks."""
    ref_pos = aln.reference_start
    query_pos = 0
    insertion_events = insertion_bases = 0
    deletion_events = deletion_bases = 0
    substitution_bases = match_bases = 0
    reference_skip_bases = 0

    for op, length in aln.cigartuples or ():
        if op in ALIGNED_QUERY_AND_REF:
            ref_next = ref_pos + length
            overlap_start = max(start, ref_pos)
            overlap_end = min(end, ref_next)
            if overlap_start < overlap_end:
                ref_offset = overlap_start - start
                q_offset = query_pos + overlap_start - ref_pos
                n = overlap_end - overlap_start
                ref_slice = reference[ref_offset : ref_offset + n]
                query_slice = aln.query_sequence[q_offset : q_offset + n]
                for ref_base, query_base in zip(ref_slice.upper(), query_slice.upper()):
                    if ref_base == query_base:
                        match_bases += 1
                    else:
                        substitution_bases += 1
            ref_pos = ref_next
            query_pos += length
            continue

        if op == CIGAR_I:
            if start <= ref_pos < end:
                insertion_events += 1
                insertion_bases += length
            query_pos += length
            continue

        if op == CIGAR_D:
            overlap_start = max(start, ref_pos)
            overlap_end = min(end, ref_pos + length)
            if overlap_start < overlap_end:
                deletion_events += 1
                deletion_bases += overlap_end - overlap_start
            ref_pos += length
            continue

        if op == CIGAR_N:
            overlap_start = max(start, ref_pos)
            overlap_end = min(end, ref_pos + length)
            if overlap_start < overlap_end:
                reference_skip_bases += overlap_end - overlap_start
            ref_pos += length
            continue

        if op == CIGAR_S:
            query_pos += length
        elif op in {CIGAR_H, CIGAR_P}:
            continue
        else:
            query_pos += length
            ref_pos += length

    return EditMetrics(
        selection_distance=selection_distance,
        insertion_events=insertion_events,
        insertion_bases=insertion_bases,
        deletion_events=deletion_events,
        deletion_bases=deletion_bases,
        substitution_bases=substitution_bases,
        match_bases=match_bases,
        reference_skip_bases=reference_skip_bases,
        query_sequence=query,
        reference_sequence=reference,
    )


def get_levenshtein_distance():
    try:
        import Levenshtein
        backend = importlib.import_module("Levenshtein.levenshtein_cpp")
    except (ImportError, ModuleNotFoundError) as exc:
        raise RuntimeError(
            "This script requires the compiled python-Levenshtein extension. "
            "Install python-Levenshtein in the environment running the script."
        ) from exc

    backend_path = str(getattr(backend, "__file__", ""))
    if not backend_path.endswith((".so", ".pyd", ".dll")):
        raise RuntimeError(
            "Levenshtein.levenshtein_cpp was not loaded from a compiled extension: "
            f"{backend_path or 'unknown path'}"
        )

    distance = getattr(Levenshtein, "distance", None)
    if not callable(distance) or distance("A", "G") != 1:
        raise RuntimeError("The imported Levenshtein implementation is unusable")
    return distance


def candidate_rank(item: ReadMetrics) -> tuple[int, int, int, int]:
    type_rank = {"primary": 0, "supplementary": 1, "secondary": 2}[item.alignment_type]
    return (item.metrics.selection_distance, type_rank, -item.mapq, -len(item.metrics.query_sequence))


def iter_read_metrics_for_locus(
    bam,
    spec: BamSpec,
    locus: Locus,
    maternal_fa,
    paternal_fa,
    args: argparse.Namespace,
    distance_fn,
) -> list[ReadMetrics]:
    start = max(0, locus.start - args.flank)
    end = locus.end + args.flank
    fetch_start = max(0, start - args.fetch_pad)
    fetch_end = end + args.fetch_pad
    best_by_read: dict[str, ReadMetrics] = {}

    for contig in bam_contigs_for_locus(bam, locus.chrom):
        for aln in bam.fetch(contig, fetch_start, fetch_end):
            if should_skip_read(aln, args) or not spans_interval(aln, start, end):
                continue
            query = collect_query_sequence(aln, start, end)
            if not query:
                continue
            references = fetch_reference(maternal_fa, paternal_fa, aln.reference_name, locus.chrom, start, end)
            if not references:
                continue

            best_candidate: ReadMetrics | None = None
            for haplotype, reference in references:
                selection_distance = distance_fn(query, reference)
                metrics = cigar_metrics(aln, start, end, reference, query, selection_distance)
                candidate = ReadMetrics(
                    locus=locus,
                    tool=spec.label,
                    read_name=get_original_read_name(aln),
                    reference_name=aln.reference_name,
                    matched_haplotype=haplotype,
                    alignment_type=alignment_type(aln),
                    mapq=aln.mapping_quality,
                    is_reverse=bool(aln.is_reverse),
                    metrics=metrics,
                )
                if best_candidate is None or candidate_rank(candidate) < candidate_rank(best_candidate):
                    best_candidate = candidate

            if best_candidate is None:
                continue
            previous = best_by_read.get(best_candidate.read_name)
            if previous is None or candidate_rank(best_candidate) < candidate_rank(previous):
                best_by_read[best_candidate.read_name] = best_candidate

    return list(best_by_read.values())


def process_locus(locus: Locus, args, maternal_fa, paternal_fa, bam_by_label, bam_specs, distance_fn):
    return {
        spec.label: iter_read_metrics_for_locus(
            bam_by_label[spec.label], spec, locus, maternal_fa, paternal_fa, args, distance_fn
        )
        for spec in bam_specs
    }


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
        "bam_specs": [BamSpec("ONT", args.ont_bam), BamSpec("PacBio", args.pacbio_bam)],
        "distance_fn": get_levenshtein_distance(),
    }


def process_locus_worker(job):
    index, locus = job
    return index, process_locus(
        locus,
        _WORKER_STATE["args"],
        _WORKER_STATE["maternal_fa"],
        _WORKER_STATE["paternal_fa"],
        _WORKER_STATE["bam_by_label"],
        _WORKER_STATE["bam_specs"],
        _WORKER_STATE["distance_fn"],
    )


def fmt(value: float | int | None) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{value:.8f}"


def mean(values: list[int | float]) -> float | None:
    return statistics.fmean(values) if values else None


def rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def write_per_read_header(writer: csv.writer) -> None:
    writer.writerow([
        "locus_id", "chrom", "start", "end", "tool", "read_name",
        "reference_name", "matched_haplotype", "alignment_type", "mapq",
        "is_reverse", "selection_edit_distance", "insertion_events",
        "insertion_bases", "deletion_events", "deletion_bases",
        "substitution_bases", "match_bases", "reference_skip_bases",
        "reference_length", "insertion_rate", "deletion_rate",
        "substitution_rate", "total_edit_rate", "query_length",
    ])


def write_per_read_row(writer: csv.writer, item: ReadMetrics) -> None:
    m = item.metrics
    denominator = len(m.reference_sequence)
    writer.writerow([
        item.locus.name, item.locus.chrom, item.locus.start, item.locus.end,
        item.tool, item.read_name, item.reference_name, item.matched_haplotype,
        item.alignment_type, item.mapq, int(item.is_reverse), m.selection_distance,
        m.insertion_events, m.insertion_bases, m.deletion_events, m.deletion_bases,
        m.substitution_bases, m.match_bases, m.reference_skip_bases, denominator,
        fmt(rate(m.insertion_bases, denominator)), fmt(rate(m.deletion_bases, denominator)),
        fmt(rate(m.substitution_bases, denominator)),
        fmt(rate(m.insertion_bases + m.deletion_bases + m.substitution_bases, denominator)),
        len(m.query_sequence),
    ])


def write_summary_header(writer: csv.writer) -> None:
    writer.writerow([
        "scope", "locus_id", "chrom", "start", "end", "tool", "n_reads",
        "n_loci_with_reads", "mean_selection_edit_distance", "mean_insertion_events",
        "mean_insertion_bases", "mean_deletion_events", "mean_deletion_bases",
        "mean_substitution_bases", "mean_match_bases", "mean_reference_skip_bases",
        "mean_reference_length", "mean_insertion_rate", "mean_deletion_rate",
        "mean_substitution_rate", "mean_total_edit_rate",
    ])


def write_summary_row(writer, scope: str, locus: Locus | None, tool: str, items: list[ReadMetrics]) -> None:
    metrics = [item.metrics for item in items]
    lengths = [len(item.metrics.reference_sequence) for item in items]
    insertion_rates = [rate(m.insertion_bases, len(m.reference_sequence)) for m in metrics]
    deletion_rates = [rate(m.deletion_bases, len(m.reference_sequence)) for m in metrics]
    substitution_rates = [rate(m.substitution_bases, len(m.reference_sequence)) for m in metrics]
    total_rates = [rate(m.insertion_bases + m.deletion_bases + m.substitution_bases, len(m.reference_sequence)) for m in metrics]
    writer.writerow([
        scope,
        locus.name if locus else "ALL",
        locus.chrom if locus else "ALL",
        locus.start if locus else "NA",
        locus.end if locus else "NA",
        tool,
        len(items),
        1 if items else 0,
        fmt(mean([m.selection_distance for m in metrics])),
        fmt(mean([m.insertion_events for m in metrics])),
        fmt(mean([m.insertion_bases for m in metrics])),
        fmt(mean([m.deletion_events for m in metrics])),
        fmt(mean([m.deletion_bases for m in metrics])),
        fmt(mean([m.substitution_bases for m in metrics])),
        fmt(mean([m.match_bases for m in metrics])),
        fmt(mean([m.reference_skip_bases for m in metrics])),
        fmt(mean(lengths)),
        fmt(mean([x for x in insertion_rates if x is not None])),
        fmt(mean([x for x in deletion_rates if x is not None])),
        fmt(mean([x for x in substitution_rates if x is not None])),
        fmt(mean([x for x in total_rates if x is not None])),
    ])


def validate_paths(args: argparse.Namespace) -> None:
    for path in (args.maternal_fa, args.paternal_fa, args.ont_bam, args.pacbio_bam):
        if not Path(path).exists():
            raise SystemExit(f"Input path does not exist: {path}")
    for fasta in (args.maternal_fa, args.paternal_fa):
        if not Path(str(fasta) + ".fai").exists():
            raise SystemExit(f"Missing FASTA index: {fasta}.fai; run samtools faidx first")
    for bam in (args.ont_bam, args.pacbio_bam):
        path = Path(bam)
        candidates = [Path(str(path) + suffix) for suffix in (".bai", ".csi")]
        candidates.extend([path.with_suffix(suffix) for suffix in (".bai", ".csi")])
        if not any(candidate.exists() for candidate in candidates):
            raise SystemExit(f"Missing BAM index for {bam}; run samtools index first")


def print_progress(idx: int, total: int, locus: Locus, completed: int | None = None) -> None:
    message = f"[cigar-edit-rates] locus {idx}/{total} {locus.chrom}:{locus.start}-{locus.end}"
    if completed is not None:
        message += f" completed={completed}/{total}"
    print(message, file=sys.stderr, flush=True)


def run_one_bed(args: argparse.Namespace) -> None:
    import pysam

    validate_paths(args)
    loci = parse_bed(Path(args.bed))
    if args.chrom:
        requested = set(args.chrom)
        loci = [locus for locus in loci if locus.chrom in requested or base_chrom(locus.chrom) in requested]
    if not loci:
        raise SystemExit(f"No loci to process in {args.bed}")

    output_path = Path(args.output_tsv)
    per_read_path = Path(args.per_read_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    per_read_path.parent.mkdir(parents=True, exist_ok=True)
    bam_specs = [BamSpec("ONT", args.ont_bam), BamSpec("PacBio", args.pacbio_bam)]
    totals: dict[str, list[ReadMetrics]] = defaultdict(list)

    with (
        pysam.FastaFile(args.maternal_fa) as maternal_fa,
        pysam.FastaFile(args.paternal_fa) as paternal_fa,
        pysam.AlignmentFile(args.ont_bam, "rb") as ont_bam,
        pysam.AlignmentFile(args.pacbio_bam, "rb") as pacbio_bam,
        open(output_path, "w", encoding="utf-8", newline="") as summary_handle,
        open(per_read_path, "w", encoding="utf-8", newline="") as per_read_handle,
    ):
        summary_writer = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
        per_read_writer = csv.writer(per_read_handle, delimiter="\t", lineterminator="\n")
        write_summary_header(summary_writer)
        write_per_read_header(per_read_writer)
        bam_by_label = {"ONT": ont_bam, "PacBio": pacbio_bam}
        distance_fn = get_levenshtein_distance()

        def record(idx: int, result: dict[str, list[ReadMetrics]]) -> None:
            locus = loci[idx - 1]
            for tool, items in result.items():
                for item in items:
                    write_per_read_row(per_read_writer, item)
                write_summary_row(summary_writer, "locus", locus, tool, items)
                totals[tool].extend(items)

        if args.processes == 1:
            for idx, locus in enumerate(loci, 1):
                result = process_locus(locus, args, maternal_fa, paternal_fa, bam_by_label, bam_specs, distance_fn)
                record(idx, result)
                if args.progress_every and (idx == 1 or idx == len(loci) or idx % args.progress_every == 0):
                    print_progress(idx, len(loci), locus)
        else:
            pending: dict[int, dict[str, list[ReadMetrics]]] = {}
            next_to_write = 1
            completed = 0
            with ProcessPoolExecutor(max_workers=args.processes, initializer=init_worker, initargs=(args,)) as executor:
                futures = {executor.submit(process_locus_worker, (idx, locus)): idx for idx, locus in enumerate(loci, 1)}
                for future in as_completed(futures):
                    idx, result = future.result()
                    pending[idx] = result
                    completed += 1
                    while next_to_write in pending:
                        record(next_to_write, pending.pop(next_to_write))
                        locus = loci[next_to_write - 1]
                        if args.progress_every and (next_to_write == 1 or next_to_write == len(loci) or next_to_write % args.progress_every == 0):
                            print_progress(next_to_write, len(loci), locus, completed)
                        next_to_write += 1

        for spec in bam_specs:
            write_summary_row(summary_writer, "tool_total", None, spec.label, totals[spec.label])
        all_items = totals["ONT"] + totals["PacBio"]
        write_summary_row(summary_writer, "all_tools_total", None, "ALL", all_items)

    print(f"[cigar-edit-rates] done bed={args.bed} processed_loci={len(loci)}", file=sys.stderr, flush=True)


def iter_bed_files(bed_dir: Path, pattern: str, recursive: bool) -> list[Path]:
    globber = bed_dir.rglob if recursive else bed_dir.glob
    return sorted(path for path in globber(pattern) if path.is_file())


def build_single_bed_args(args: argparse.Namespace, bed_path: Path) -> argparse.Namespace:
    stem = bed_path.name[:-4] if bed_path.name.endswith(".bed") else bed_path.stem
    return argparse.Namespace(
        bed=str(bed_path), maternal_fa=args.maternal_fa, paternal_fa=args.paternal_fa,
        ont_bam=args.ont_bam, pacbio_bam=args.pacbio_bam,
        output_tsv=str(args.output_dir / f"{stem}{args.summary_suffix}"),
        per_read_tsv=str(args.output_dir / f"{stem}{args.per_read_suffix}"),
        flank=args.flank, fetch_pad=args.fetch_pad, min_mapq=args.min_mapq,
        primary_only=args.primary_only, include_secondary=args.include_secondary,
        chrom=args.chrom, progress_every=args.progress_every, processes=args.processes,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Calculate CIGAR-based ONT/PacBio edit rates for BED loci.")
    parser.add_argument("--bed-dir", required=True, type=Path)
    parser.add_argument("--bed-pattern", default="*.bed")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--maternal-fa", required=True)
    parser.add_argument("--paternal-fa", required=True)
    parser.add_argument("--ont-bam", required=True)
    parser.add_argument("--pacbio-bam", required=True)
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    parser.add_argument("--summary-suffix", default=".cigar_edit_rates.summary.tsv")
    parser.add_argument("--per-read-suffix", default=".cigar_edit_rates.per_read.tsv")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--flank", type=int, default=0)
    parser.add_argument("--fetch-pad", type=int, default=50)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--primary-only", action="store_true")
    parser.add_argument("--include-secondary", action="store_true")
    parser.add_argument("--chrom", action="append")
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--processes", type=int, default=1)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.bed_dir.is_dir():
        raise SystemExit(f"BED directory does not exist: {args.bed_dir}")
    if args.flank < 0 or args.fetch_pad < 0 or args.min_mapq < 0 or args.processes < 1 or args.progress_every < 0:
        raise SystemExit("--flank, --fetch-pad, --min-mapq and --processes must be valid non-negative values")
    get_levenshtein_distance()
    bed_files = iter_bed_files(args.bed_dir, args.bed_pattern, args.recursive)
    if not bed_files:
        raise SystemExit(f"No BED files matched {args.bed_pattern!r} in {args.bed_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    failures: list[tuple[Path, BaseException]] = []
    for index, bed_path in enumerate(bed_files, 1):
        single_args = build_single_bed_args(args, bed_path)
        summary_path = Path(single_args.output_tsv)
        if args.skip_existing and summary_path.exists():
            print(f"[batch-cigar-edit-rates] SKIP {index}/{len(bed_files)} {bed_path}", file=sys.stderr, flush=True)
            continue
        print(f"[batch-cigar-edit-rates] START {index}/{len(bed_files)} {bed_path}", file=sys.stderr, flush=True)
        try:
            run_one_bed(single_args)
        except (Exception, SystemExit) as exc:
            failures.append((bed_path, exc))
            print(f"[batch-cigar-edit-rates] ERROR {bed_path}: {exc}", file=sys.stderr, flush=True)
            if not args.continue_on_error:
                raise
    if failures:
        print(f"[batch-cigar-edit-rates] failed_beds={len(failures)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
