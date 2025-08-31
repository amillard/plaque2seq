#!/usr/bin/env python3
#Map Reads to Viral Contigs and Compute Coverage

import os
import glob
import argparse
import subprocess
import pysam
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

def read_sample_names(dict_file):
    sample_names = set()
    with open(dict_file) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    sample_names.add(parts[1])
    return sample_names

def find_fasta_matches(clean_dir, sample_names):
    fasta_files = glob.glob(os.path.join(clean_dir, "*.fa")) + glob.glob(os.path.join(clean_dir, "*.fasta"))
    matches = []
    for sample in sample_names:
        for fasta in fasta_files:
            if sample in os.path.basename(fasta):
                matches.append((sample, fasta))
    return matches

def run_minimap2(sample, fasta_file, reads_dir, bam_dir):
    reads_file = os.path.join(reads_dir, f"{sample}_short_300.fq.gz")
    if not os.path.isfile(reads_file):
        return sample, fasta_file, None, f"❌ Reads file not found: {reads_file}"

    os.makedirs(bam_dir, exist_ok=True)
    #bam_name = f"{sample}_{os.path.basename(fasta_file).replace('.fasta','').replace('.fa','')}.bam"
    bam_name = f"{os.path.basename(fasta_file).replace('.fasta','').replace('.fa','')}.bam"
    bam_path = os.path.join(bam_dir, bam_name)

    if os.path.exists(bam_path) and os.path.exists(bam_path + ".bai"):
        return sample, fasta_file, bam_path, None

    cmd = (
        f"minimap2 -t 4 -a {fasta_file} {reads_file} | "
        f"samtools sort -@ 2 -o {bam_path} && "
        f"samtools index {bam_path}"
    )

    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        return sample, fasta_file, bam_path, None
    except subprocess.CalledProcessError as e:
        return sample, fasta_file, None, f"❌ Mapping failed: {e}"

def compute_bam_stats(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    total_reads = 0
    mapped_reads = 0
    total_covered_bases = 0
    total_ref_length = 0
    all_depths = []

    for read in bam.fetch(until_eof=True):
        total_reads += 1
        if not read.is_unmapped:
            mapped_reads += 1
    unmapped_reads = total_reads - mapped_reads

    for ref in bam.references:
        ref_len = bam.get_reference_length(ref)
        total_ref_length += ref_len
        coverage = bam.count_coverage(contig=ref)
        per_base_depth = np.sum(np.array(coverage), axis=0) // 4
        total_covered_bases += per_base_depth.sum()
        all_depths.extend(per_base_depth.tolist())

    bam.close()

    pct_mapped = (mapped_reads / total_reads * 100) if total_reads else 0
    mean_cov = (total_covered_bases / total_ref_length) if total_ref_length else 0
    median_cov = np.median(all_depths) if all_depths else 0

    return total_reads, mapped_reads, unmapped_reads, pct_mapped, mean_cov, median_cov

def worker(args):
    sample, fasta_file, reads_dir, bam_dir = args
    sample, fasta_file, bam_path, error = run_minimap2(sample, fasta_file, reads_dir, bam_dir)
    if error or not bam_path:
        return (sample, fasta_file, None, None, None, None, None, None, error)
    try:
        total, mapped, unmapped, pct, mean, median = compute_bam_stats(bam_path)
        return (sample, fasta_file, os.path.basename(bam_path), total, mapped, unmapped, pct, mean, median, None)
    except Exception as e:
        return (sample, fasta_file, os.path.basename(bam_path), None, None, None, None, None, None, str(e))

def main():
    parser = argparse.ArgumentParser(description="Map reads to contigs and compute stats (parallelized).")
    parser.add_argument("--dict", required=True, help="2-column tab file with sample names in column 2")
    parser.add_argument("--clean_dir", required=True, help="Directory with *.fa / *.fasta")
    parser.add_argument("--reads_dir", required=True, help="Directory with *_short_300.fq.gz files")
    parser.add_argument("--bam_dir", required=True, help="Output directory for BAMs")
    parser.add_argument("--out_tsv", default="mapping_summary.tsv", help="Summary output TSV")
    parser.add_argument("-w", "--workers", type=int, default=4, help="Number of parallel workers")

    args = parser.parse_args()

    sample_names = read_sample_names(args.dict)
    matches = find_fasta_matches(args.clean_dir, sample_names)

    if not matches:
        print("❌ No matching FASTA files found.")
        return

    print(f"🔍 {len(matches)} matching FASTA files found. Starting parallel processing...\n")

    tasks = [(sample, fasta, args.reads_dir, args.bam_dir) for sample, fasta in matches]

    with ProcessPoolExecutor(max_workers=args.workers) as executor, open(args.out_tsv, "w") as out:
        out.write("Sample\tContigFile\tBAM\tTotalReads\tMappedReads\tUnmappedReads\tPctMapped\tMeanCoverage\tMedianCoverage\n")
        futures = {executor.submit(worker, t): t for t in tasks}
        for future in as_completed(futures):
            sample, fasta_file, bam, total, mapped, unmapped, pct, mean, median, error = future.result()
            if error:
                print(f"❌ {sample} → Error: {error}")
                continue
            out.write(f"{sample}\t{os.path.basename(fasta_file)}\t{bam}\t{total}\t{mapped}\t{unmapped}\t{pct:.2f}\t{mean:.2f}\t{median:.2f}\n")
            print(f"✅ {sample}: mapped={mapped}/{total}, unmapped={unmapped}, mean={mean:.2f}x, median={median:.2f}x")

if __name__ == "__main__":
    main()
