#!/usr/bin/env python3
"""Merge windowed TRF .dat files into chromosome/haplotype BED files.

Input .dat filenames are expected to look like:

    chr1_MATERNAL_1_500000.fa.2.7.7.80.10.50.2000.dat

TRF coordinates inside each .dat file are 1-based relative to the window
sequence. This script converts them back to chromosome/haplotype coordinates
and writes tab-delimited BED-like files grouped by chr/haplotype.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DAT_NAME_RE = re.compile(
    r"^(?P<chrom>.+)_(?P<hap>MATERNAL|PATERNAL)_(?P<start>\d+)_(?P<end>\d+)"
    r"\.(?:fa|fasta|fna)(?:\..*)?\.dat$"
)


@dataclass(frozen=True)
class WindowDat:
    path: Path
    chrom: str
    haplotype: str
    window_start1: int
    window_end1: int

    @property
    def chrom_haplotype(self) -> str:
        return f"{self.chrom}_{self.haplotype}"


@dataclass(frozen=True)
class BedRow:
    chrom: str
    start0: int
    end0: int
    fields: tuple[str, ...]

    @property
    def sort_key(self) -> tuple[str, int, int, tuple[str, ...]]:
        return self.chrom, self.start0, self.end0, self.fields

    @property
    def region_key(self) -> tuple[str, int, int]:
        return self.chrom, self.start0, self.end0

    @property
    def trf_score(self) -> float:
        if len(self.fields) <= 10:
            return 0.0
        try:
            return float(self.fields[10])
        except ValueError:
            return 0.0

    @property
    def region_dedup_sort_key(self) -> tuple[str, int, int, float, tuple[str, ...]]:
        return self.chrom, self.start0, self.end0, -self.trf_score, self.fields


def parse_dat_name(path: Path) -> WindowDat | None:
    match = DAT_NAME_RE.match(path.name)
    if not match:
        return None
    return WindowDat(
        path=path,
        chrom=match.group("chrom"),
        haplotype=match.group("hap"),
        window_start1=int(match.group("start")),
        window_end1=int(match.group("end")),
    )


def iter_dat_files(dat_dir: Path, pattern: str, recursive: bool) -> list[Path]:
    globber = dat_dir.rglob if recursive else dat_dir.glob
    return sorted(path for path in globber(pattern) if path.is_file())


def is_trf_data_line(line: str) -> bool:
    return bool(re.match(r"^\d+\s", line))


def read_trf_dat(dat: WindowDat) -> Iterable[BedRow]:
    with dat.path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or not is_trf_data_line(line):
                continue
            trf_fields = line.split()
            if len(trf_fields) < 3:
                raise ValueError(f"{dat.path}:{line_no} has fewer than 3 TRF columns")
            try:
                rel_start1 = int(trf_fields[0])
                rel_end1 = int(trf_fields[1])
            except ValueError as exc:
                raise ValueError(f"{dat.path}:{line_no} has non-integer TRF start/end") from exc
            if rel_end1 < rel_start1:
                raise ValueError(f"{dat.path}:{line_no} has TRF end < start")

            global_start1 = dat.window_start1 + rel_start1 - 1
            global_end1 = dat.window_start1 + rel_end1 - 1
            start0 = global_start1 - 1
            end0 = global_end1

            yield BedRow(
                chrom=dat.chrom,
                start0=start0,
                end0=end0,
                fields=(
                    dat.chrom,
                    str(start0),
                    str(end0),
                    str(global_start1),
                    str(global_end1),
                    *trf_fields[2:],
                    dat.path.name,
                ),
            )


def dedup_rows(rows: Iterable[BedRow]) -> list[BedRow]:
    seen: set[tuple[str, ...]] = set()
    deduped: list[BedRow] = []
    for row in rows:
        if row.fields in seen:
            continue
        seen.add(row.fields)
        deduped.append(row)
    return deduped


def dedup_rows_by_region(rows: Iterable[BedRow]) -> list[BedRow]:
    seen: set[tuple[str, int, int]] = set()
    deduped: list[BedRow] = []
    for row in sorted(rows, key=lambda item: item.region_dedup_sort_key):
        if row.region_key in seen:
            continue
        seen.add(row.region_key)
        deduped.append(row)
    return deduped


def write_bed(path: Path, rows: Iterable[BedRow]) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow(row.fields)
            count += 1
    return count


def write_summary(path: Path, summary_rows: list[tuple[str, int, int]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom_haplotype", "dat_files", "trf_records"])
        writer.writerows(summary_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Merge windowed TRF .dat files into tab-delimited BED files grouped by chromosome/haplotype."
    )
    parser.add_argument("--dat-dir", required=True, type=Path, help="Directory containing TRF .dat files.")
    parser.add_argument("-o", "--out-dir", required=True, type=Path, help="Output directory for grouped BED files.")
    parser.add_argument(
        "--dat-pattern",
        default="*.dat",
        help="Glob pattern for .dat files inside --dat-dir. Default: *.dat",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search for .dat files recursively under --dat-dir.",
    )
    parser.add_argument(
        "--dedup-exact",
        action="store_true",
        help="Remove exact duplicate BED rows after coordinate conversion.",
    )
    parser.add_argument(
        "--dedup-region",
        action="store_true",
        help=(
            "For duplicate chrom/start/end regions, keep one row with the highest "
            "TRF score in BED column 11, matching sort -k1,1 -k2,2n -k3,3n -k11,11nr."
        ),
    )
    parser.add_argument(
        "--summary-name",
        default="merge_window_trf_dat_to_bed.summary.tsv",
        help="Summary TSV filename written inside --out-dir.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.dat_dir.is_dir():
        raise SystemExit(f"DAT directory does not exist: {args.dat_dir}")

    dat_paths = iter_dat_files(args.dat_dir, args.dat_pattern, args.recursive)
    if not dat_paths:
        raise SystemExit(f"No .dat files matched {args.dat_pattern!r} in {args.dat_dir}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    grouped_rows: dict[str, list[BedRow]] = {}
    grouped_dat_counts: dict[str, int] = {}
    skipped_names: list[Path] = []

    for dat_path in dat_paths:
        dat = parse_dat_name(dat_path)
        if dat is None:
            skipped_names.append(dat_path)
            continue
        key = dat.chrom_haplotype
        grouped_dat_counts[key] = grouped_dat_counts.get(key, 0) + 1
        grouped_rows.setdefault(key, []).extend(read_trf_dat(dat))

    summary_rows: list[tuple[str, int, int]] = []
    for key in sorted(grouped_rows):
        rows = sorted(grouped_rows[key], key=lambda row: row.sort_key)
        if args.dedup_region:
            rows = dedup_rows_by_region(rows)
        if args.dedup_exact:
            rows = dedup_rows(rows)
        rows = sorted(rows, key=lambda row: row.sort_key)
        out_path = args.out_dir / f"{key}.trf.raw.bed"
        record_count = write_bed(out_path, rows)
        summary_rows.append((key, grouped_dat_counts.get(key, 0), record_count))
        print(
            f"[merge-window-trf] wrote {out_path} dat_files={grouped_dat_counts.get(key, 0)} records={record_count}",
            flush=True,
        )

    write_summary(args.out_dir / args.summary_name, summary_rows)
    print(
        (
            f"[merge-window-trf] finished dat_files={len(dat_paths)} "
            f"groups={len(grouped_rows)} skipped_unrecognized_names={len(skipped_names)}"
        ),
        flush=True,
    )
    if skipped_names:
        print("[merge-window-trf] first skipped filenames:", flush=True)
        for path in skipped_names[:10]:
            print(f"[merge-window-trf]   {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
