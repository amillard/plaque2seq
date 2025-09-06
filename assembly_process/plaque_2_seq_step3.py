#!/usr/bin/env python3

#assemble genomes multiple at the same time

import os
import glob
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

def assemble_reads(input_file: str, threads: int, output_base: str) -> str:
    """
    Assemble single-end normalized reads with SPAdes (--only-assembler, single-cell mode).

    Parameters:
        input_file:   Path to a *_norm.fq.gz file.
        threads:      Number of threads to pass to SPAdes.
        output_base:  Base directory to store output folders.

    Returns:
        Path to the SPAdes output directory.
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    fname = os.path.basename(input_file)
    if "_short_" not in fname or "_norm" not in fname:
        raise ValueError(f"Cannot parse sample name from {fname!r}")

    sample = fname.split("_short_")[0]
    out_dir = os.path.join(output_base, f"ass_{sample}")
    os.makedirs(output_base, exist_ok=True)

    if os.path.exists(out_dir):
        print(f"    Skipping assembly for {input_file!r}: {out_dir!r} already exists.")
        return out_dir

    cmd = [
        "spades.py",
        "-s", input_file,
        "--only-assembler",
        "-t", str(threads),
        "-o", out_dir,
    ]

    print(f" Assembling {input_file!r} → {out_dir!r} with {threads} threads")
    subprocess.run(cmd, check=True)
    print(f"Assembly complete: {out_dir!r}")
    return out_dir

def main():
    parser = argparse.ArgumentParser(description="Parallel SPAdes assembly of *_norm.fq.gz files.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing *_norm.fq.gz files")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Threads per SPAdes job (default: 4) More is less -or quciker ")
    parser.add_argument("-o", "--output_base", default=".", help="Base output directory (default: current dir)")
    parser.add_argument("-w", "--workers", type=int, default=None, help="Max number of concurrent assemblies (default: all cores)")

    args = parser.parse_args()

    pattern = os.path.join(args.directory, "*_norm.fq.gz")
    files = sorted(glob.glob(pattern))

    if not files:
        print(" No *_norm.fq.gz files found.")
        return

    print(f" Found {len(files)} files. Starting parallel assembly...")

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_file = {
            executor.submit(assemble_reads, f, args.threads, args.output_base): f for f in files
        }

        for future in as_completed(future_to_file):
            f = future_to_file[future]
            try:
                result = future.result()
            except Exception as e:
                print(f" Failed to assemble {f}: {e}")

if __name__ == "__main__":
    main()
