#!/usr/bin/env python3

#merge multiple phage genomes into a single file for mapping data. Part of plaque-2-seq pipeline 


import os
import argparse
from collections import defaultdict
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Merge FASTA files by sample base into a new directory.")
    parser.add_argument("--input_dir", default=".", help="Directory containing FASTA files (default: current dir)")
    parser.add_argument("--output_dir", default="merged", help="Directory to write merged files (default: merged/)")
    return parser.parse_args()

def get_sample_name(filename):
    """Extract sample name from filename: join first two underscore-separated parts"""
    parts = filename.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else None

def find_fasta_files(input_dir):
    """Find .fa or .fasta files in the input directory"""
    fasta_files = []
    for ext in ("*.fa", "*.fasta"):
        fasta_files.extend(Path(input_dir).glob(ext))
    return fasta_files

def merge_fastas_by_sample(fasta_files, output_dir):
    """Group and merge FASTA files by sample name"""
    sample_to_files = defaultdict(list)

    for f in fasta_files:
        sample = get_sample_name(f.name)
        if sample:
            sample_to_files[sample].append(f)

    os.makedirs(output_dir, exist_ok=True)

    for sample, files in sample_to_files.items():
        output_path = Path(output_dir) / f"{sample}_catenated.fa"
        with open(output_path, "w") as out_fh:
            for f in sorted(files):  # optional sort by filename
                with open(f, "r") as in_fh:
                    out_fh.write(in_fh.read())
        print(f"Merged {len(files)} file(s) into {output_path}")

def main():
    args = parse_args()
    fasta_files = find_fasta_files(args.input_dir)

    if not fasta_files:
        print(f"o FASTA files found in {args.input_dir}")
        return

    merge_fastas_by_sample(fasta_files, args.output_dir)

if __name__ == "__main__":
    main()
