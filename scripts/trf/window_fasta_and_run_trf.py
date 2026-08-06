#!/usr/bin/env python3
"""Split FASTA files into overlapping windows and run TRF in parallel."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


FASTA_EXTENSIONS = (".fa", ".fasta", ".fna")
COORD_SUFFIX_RE = re.compile(r"^(?P<name>.+)_(?P<start>\d+)_(?P<end>\d+)$")


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str


@dataclass(frozen=True)
class WindowFile:
    path: Path
    source_fasta: Path
    chrom: str
    start1: int
    end1: int


def iter_fasta_files(input_dir: Path, recursive: bool) -> list[Path]:
    pattern = "**/*" if recursive else "*"
    return sorted(
        path
        for path in input_dir.glob(pattern)
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    name: str | None = None
    pieces: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append(FastaRecord(name=name, sequence="".join(pieces)))
                name = line[1:].split()[0]
                pieces = []
            else:
                pieces.append(line)
    if name is not None:
        records.append(FastaRecord(name=name, sequence="".join(pieces)))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def fasta_base(path: Path) -> str:
    name = path.name
    for suffix in FASTA_EXTENSIONS:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def source_name_and_offset(path: Path, record: FastaRecord) -> tuple[str, int]:
    base = fasta_base(path)
    match = COORD_SUFFIX_RE.match(base)
    if match:
        return match.group("name"), int(match.group("start"))
    return record.name, 1


def wrap_sequence(seq: str, width: int = 80) -> Iterable[str]:
    for idx in range(0, len(seq), width):
        yield seq[idx : idx + width]


def write_window_fasta(
    path: Path,
    chrom: str,
    start1: int,
    end1: int,
    sequence: str,
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{chrom}_{start1}_{end1}\n")
        for line in wrap_sequence(sequence):
            handle.write(line + "\n")


def split_fasta_file(
    fasta_path: Path,
    windows_dir: Path,
    window_size: int,
    overlap: int,
    force: bool,
) -> list[WindowFile]:
    if overlap >= window_size:
        raise ValueError("--overlap must be smaller than --window-size")
    step = window_size - overlap
    windows: list[WindowFile] = []
    records = read_fasta(fasta_path)

    for record in records:
        chrom, source_start1 = source_name_and_offset(fasta_path, record)
        sequence = record.sequence.upper()
        out_dir = windows_dir / f"{fasta_base(fasta_path)}_subwindows_overlap"
        out_dir.mkdir(parents=True, exist_ok=True)

        for start0 in range(0, len(sequence), step):
            end0 = min(start0 + window_size, len(sequence))
            if end0 <= start0:
                continue
            start1 = source_start1 + start0
            end1 = source_start1 + end0 - 1
            out_path = out_dir / f"{chrom}_{start1}_{end1}.fa"
            if force or not out_path.exists():
                write_window_fasta(out_path, chrom, start1, end1, sequence[start0:end0])
            windows.append(
                WindowFile(
                    path=out_path,
                    source_fasta=fasta_path,
                    chrom=chrom,
                    start1=start1,
                    end1=end1,
                )
            )
            if end0 == len(sequence):
                break
    return windows


def trf_dat_name(window_path: Path, trf_params: tuple[str, ...]) -> str:
    return f"{window_path.name}.{'.'.join(trf_params)}.dat"


def run_trf(
    window: WindowFile,
    trf_bin: str,
    trf_params: tuple[str, ...],
    trf_flags: tuple[str, ...],
    trf_out_dir: Path,
    skip_existing: bool,
) -> tuple[WindowFile, int, str]:
    trf_out_dir.mkdir(parents=True, exist_ok=True)
    window_path = window.path.resolve()
    dat_path = trf_out_dir / trf_dat_name(window_path, trf_params)
    if skip_existing and dat_path.exists():
        return window, 0, "skipped"

    cmd = [trf_bin, str(window_path), *trf_params, *trf_flags]
    result = subprocess.run(
        cmd,
        cwd=trf_out_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    combined_output = "\n".join(part for part in (result.stdout, result.stderr) if part)
    if result.returncode != 0 and dat_path.exists() and "Done." in combined_output:
        return window, 0, f"done_returncode_{result.returncode}"
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        if not message:
            message = "TRF exited with no stdout/stderr"
        return window, result.returncode, message
    return window, 0, "done"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Split all FASTA files in a directory into overlapping windows, then run TRF in parallel."
    )
    parser.add_argument("--input-dir", required=True, type=Path, help="Directory containing .fa/.fasta/.fna files.")
    parser.add_argument("--windows-dir", required=True, type=Path, help="Directory for window FASTA files.")
    parser.add_argument("--trf-out-dir", required=True, type=Path, help="Directory for TRF output files.")
    parser.add_argument("--window-size", type=int, default=500_000, help="Window size in bp. Default: 500000.")
    parser.add_argument("--overlap", type=int, default=50_000, help="Window overlap in bp. Default: 50000.")
    parser.add_argument("--recursive", action="store_true", help="Search input FASTA files recursively.")
    parser.add_argument(
        "--skip-windowing",
        action="store_true",
        help="Skip splitting input FASTA files and run TRF on existing FASTA windows in --windows-dir.",
    )
    parser.add_argument("--jobs", type=int, default=4, help="Number of parallel TRF jobs. Default: 4.")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="Print progress every N finished windows. Default: 100.",
    )
    parser.add_argument("--trf-bin", default="trf", help="TRF executable name/path. Default: trf.")
    parser.add_argument(
        "--trf-params",
        nargs=7,
        default=("2", "7", "7", "80", "10", "50", "2000"),
        metavar=("MATCH", "MISMATCH", "DELTA", "PM", "PI", "MINSCORE", "MAXPERIOD"),
        help="Seven TRF scoring/search parameters. Default: 2 7 7 80 10 50 2000.",
    )
    parser.add_argument(
        "--trf-flags",
        nargs="*",
        default=("-d", "-h"),
        help="Extra TRF flags. Default: -d -h.",
    )
    parser.add_argument("--skip-existing", action="store_true", help="Skip TRF if expected .dat output exists.")
    parser.add_argument("--force-windows", action="store_true", help="Overwrite existing window FASTA files.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.input_dir.is_dir() and not args.skip_windowing:
        raise SystemExit(f"Input directory does not exist: {args.input_dir}")
    if args.window_size <= 0:
        raise SystemExit("--window-size must be > 0")
    if args.overlap < 0 or args.overlap >= args.window_size:
        raise SystemExit("--overlap must be >= 0 and smaller than --window-size")
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")
    if args.progress_every < 1:
        raise SystemExit("--progress-every must be >= 1")

    args.input_dir = args.input_dir.resolve()
    args.windows_dir = args.windows_dir.resolve()
    args.trf_out_dir = args.trf_out_dir.resolve()

    all_windows: list[WindowFile] = []
    if args.skip_windowing:
        if not args.windows_dir.is_dir():
            raise SystemExit(f"Windows directory does not exist: {args.windows_dir}")
        window_files = iter_fasta_files(args.windows_dir, recursive=True)
        if not window_files:
            raise SystemExit(f"No existing FASTA windows found in {args.windows_dir}")
        all_windows = [
            WindowFile(
                path=path,
                source_fasta=path,
                chrom=fasta_base(path),
                start1=0,
                end1=0,
            )
            for path in window_files
        ]
        fasta_files = []
        print(
            f"[window-trf] using existing windows from {args.windows_dir} windows={len(all_windows)}",
            file=sys.stderr,
            flush=True,
        )
    else:
        fasta_files = iter_fasta_files(args.input_dir, args.recursive)
        if not fasta_files:
            raise SystemExit(f"No FASTA files found in {args.input_dir}")
        for fasta_path in fasta_files:
            windows = split_fasta_file(
                fasta_path,
                args.windows_dir,
                args.window_size,
                args.overlap,
                args.force_windows,
            )
            all_windows.extend(windows)
            print(
                f"[window-trf] split {fasta_path} windows={len(windows)}",
                file=sys.stderr,
                flush=True,
            )

    failures: list[tuple[WindowFile, int, str]] = []
    completed = 0
    skipped = 0
    trf_params = tuple(args.trf_params)
    trf_flags = tuple(args.trf_flags)

    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_map = {
            executor.submit(
                run_trf,
                window,
                args.trf_bin,
                trf_params,
                trf_flags,
                args.trf_out_dir,
                args.skip_existing,
            ): window
            for window in all_windows
        }
        for future in as_completed(future_map):
            window, returncode, message = future.result()
            if returncode == 0 and message == "skipped":
                skipped += 1
            elif returncode == 0:
                completed += 1
            else:
                failures.append((window, returncode, message))
                if len(failures) <= 5:
                    print(
                        (
                            f"[window-trf] failed returncode={returncode} "
                            f"window={window.path} message={message}"
                        ),
                        file=sys.stderr,
                        flush=True,
                    )
            done = completed + skipped + len(failures)
            if done == 1 or done == len(all_windows) or done % args.progress_every == 0:
                print(
                    (
                        f"[window-trf] progress {done}/{len(all_windows)} "
                        f"completed={completed} skipped={skipped} failed={len(failures)}"
                    ),
                    file=sys.stderr,
                    flush=True,
                )

    print(
        (
            f"[window-trf] finished fasta_files={len(fasta_files)} "
            f"windows={len(all_windows)} completed={completed} skipped={skipped} "
            f"failed={len(failures)}"
        ),
        file=sys.stderr,
        flush=True,
    )
    if failures:
        for window, returncode, message in failures[:20]:
            print(
                f"[window-trf] failed returncode={returncode} window={window.path} message={message}",
                file=sys.stderr,
            )
        raise SystemExit("Some TRF jobs failed; see stderr for details.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
