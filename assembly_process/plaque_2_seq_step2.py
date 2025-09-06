#!/usr/bin/env python3
#normalise reads for faster assembly with spades

import os
import subprocess
import argparse
import glob

def normalize_reads(input_file: str, target: int = 150, output_dir: str = None) -> str:
    """
    Normalize read coverage using bbnorm.sh, forcing qin=33 and ignoring bad quality scores.
    """
    if not input_file.endswith(".fq.gz"):
        raise ValueError(f"Unsupported extension: {input_file!r}")

    base = os.path.basename(input_file)[:-len(".fq.gz")]
    out_dir = output_dir or os.path.dirname(input_file) or "."
    os.makedirs(out_dir, exist_ok=True)
    output_file = os.path.join(out_dir, f"{base}_norm.fq.gz")

    if os.path.exists(output_file):
        print(f"Skipping {input_file!r}: {output_file!r} already exists.")
        return output_file

    cmd = (
        f"bbnorm.sh in={input_file} out={output_file} "
        f"target={target} qin=33 ignorebadquality"
    )

    print(f"Normalizing {input_file!r} → {output_file!r} (target={target})")
    try:
        subprocess.run(["bash", "-lc", cmd], check=True)
        print(f"Normalization complete: {output_file!r}")
    except subprocess.CalledProcessError as e:
        print(f"Normalization failed for {input_file!r}\nError: {e}")
        raise

    return output_file

def main():
    parser = argparse.ArgumentParser(description="Normalize all *short*_*.fq.gz reads in a directory using bbnorm.sh.")
    parser.add_argument("input_dir", help="Directory containing FASTQ files")
    parser.add_argument("-t", "--target", type=int, default=150, help="Target coverage (default: 150)")
    parser.add_argument("-o", "--output_dir", default=None, help="Optional output directory (default: input_dir)")

    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir or input_dir
    pattern = os.path.join(input_dir, "*short*_*.fq.gz")
    files = sorted(glob.glob(pattern))

    if not files:
        print(f"!No files found matching pattern: {pattern}")
        return

    for fq in files:
        if "norm" in os.path.basename(fq):
            print(f"  Skipping {fq} (already normalized)")
            continue
        normalize_reads(fq, target=args.target, output_dir=output_dir)

if __name__ == "__main__":
    main()
