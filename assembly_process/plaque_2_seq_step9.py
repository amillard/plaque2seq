#!/usr/bin/env python3

""""Runs pilon on a set of bam files, matching based on bc_2sample file as usual
corrects any issues.In testing didnt find errors when > 30x was used """

import os
import glob
import argparse
import subprocess

def read_sample_names(dict_file):
    """Read unique sample names from the 2nd column of a tab-delimited file."""
    samples = set()
    with open(dict_file) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    samples.add(parts[1])
    return samples

def find_matching_pairs(sample, clean_dir, bam_dir):
    """Match FASTA files and BAM files based on sample name and shared prefixes."""
    fasta_candidates = glob.glob(os.path.join(clean_dir, "*.fa")) + glob.glob(os.path.join(clean_dir, "*.fasta"))
    bam_candidates = glob.glob(os.path.join(bam_dir, "*.bam"))
    pairs = []

    for fasta in fasta_candidates:
        fasta_base = os.path.basename(fasta)
        if sample not in fasta_base:
            continue
        fasta_core = os.path.splitext(fasta_base)[0]

        # Match BAM based on FASTA core
        matching_bam = next((bam for bam in bam_candidates if fasta_core in os.path.basename(bam)), None)
        if matching_bam:
            pairs.append((sample, fasta, matching_bam))

    return pairs

def run_pilon(sample, fasta, bam, outdir):
    """Run Pilon unless output already exists."""
    os.makedirs(outdir, exist_ok=True)

    fasta_core = os.path.splitext(os.path.basename(fasta))[0]
    output_prefix = fasta_core
    expected_output = os.path.join(outdir, f"{output_prefix}.fasta")
    changes_file = os.path.join(outdir, f"{sample}.changes.txt")

    if os.path.exists(expected_output):
        print(f"Skipping {output_prefix}: output already exists.")
        return

    cmd = [
        "pilon",
        "--genome", fasta,
        "--unpaired", bam,
        "--output", output_prefix,
        "--outdir", outdir,
        "--changes"
    ]

    try:
        print(f"Running Pilon: sample={sample}, genome={os.path.basename(fasta)}, bam={os.path.basename(bam)}")
        subprocess.run(cmd, check=True)
        print(f"Finished: {output_prefix}, changes → {changes_file}")
    except subprocess.CalledProcessError as e:
        print(f" Pilon failed for {output_prefix}:\n{e}")

def main():
    parser = argparse.ArgumentParser(description="Run Pilon on genome/BAM pairs per sample.")
    parser.add_argument("--dict", required=True, help="Tab-delimited file, sample names in column 2")
    parser.add_argument("--clean_dir", required=True, help="Directory containing *.fa or *.fasta files")
    parser.add_argument("--bam_dir", required=True, help="Directory containing matching BAM files")
    parser.add_argument("--outdir", default="pilon", help="Output directory for Pilon results")

    args = parser.parse_args()

    samples = read_sample_names(args.dict)

    for sample in sorted(samples):
        pairs = find_matching_pairs(sample, args.clean_dir, args.bam_dir)
        if not pairs:
            print(f" No FASTA/BAM pair found for {sample}")
            continue

        for _, fasta, bam in pairs:
            run_pilon(sample, fasta, bam, args.outdir)

if __name__ == "__main__":
    main()
