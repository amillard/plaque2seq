#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import glob
import tempfile
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import pysam
import gzip
import csv

def parse_args():
    p = argparse.ArgumentParser(
        description="Map multiple long-read FASTQ files against multiple references and report %% mapped per reference"
    )
    p.add_argument("-i", "--input", required=True, help="Input directory of FASTQ files")
    p.add_argument("--pattern", default="*.fq.gz", help="Glob pattern to match FASTQ files (default: *.fq.gz)")
    grp = p.add_mutually_exclusive_group(required=False)
    grp.add_argument("-r", "--references", nargs="+", help="One or more reference FASTA files")
    grp.add_argument("-d", "--ref-dir", help="Directory containing reference FASTA files (*.fna)")
    p.add_argument("-a", "--accession", help="NCBI nucleotide accession to download")
    p.add_argument("-t", "--threads", type=int, default=1, help="Threads per minimap2 run (default: 1)")
    p.add_argument("-w", "--workers", type=int, default=1, help="Number of concurrent workers (default: 1)")
    p.add_argument("-p", "--preset", choices=["map-ont", "map-pb"], default="map-ont", help="minimap2 preset (default: map-ont)")
    p.add_argument("--output", default="Genome_contamination_mapping_summary.tsv", help="TSV output file (default: mapping_summary.tsv)")
    p.add_argument("--dry-run", action="store_true", help="Preview files and references without mapping")
    return p.parse_args()

def count_total_reads(fastq_path):
    opener = gzip.open if fastq_path.endswith(".gz") else open
    mode = "rt" if fastq_path.endswith(".gz") else "r"
    with opener(fastq_path, mode) as handle:
        return sum(1 for _ in SeqIO.parse(handle, "fastq"))

def download_reference(accession, output_dir="downloaded_refs"):
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{accession}.fna")
    if os.path.exists(output_path):
        print(f"Reference {accession} already exists. Skipping download.")
        return output_path
    print(f"Downloading reference {accession}...")
    cmd = f"esearch -db nucleotide -query {accession} | efetch -format fasta > {output_path}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        sys.exit(f" Failed to download accession {accession}")
    return output_path

def map_and_count_mapped(fastq_path, ref_path, preset, threads):
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_file = os.path.join(tmpdir, "out.bam")
        cmd = [
            "minimap2", "-a", "--secondary=no", "-x", preset,
            "-t", str(threads), ref_path, fastq_path
        ]
        cmd_str = f"{' '.join(cmd)} | samtools view -Sb - | samtools sort -o {bam_file} -"
        try:
            subprocess.run(cmd_str, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            sys.exit(f" minimap2 or samtools failed on {ref_path}:\n{e.stderr.decode()}")
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            return sum(1 for r in bam.fetch(until_eof=True) if not r.is_unmapped)

def process_one(fastq, ref_path, preset, threads):
    total_reads = count_total_reads(fastq)
    if total_reads == 0:
        return (os.path.basename(fastq), os.path.basename(ref_path), 0, 0, 0)
    mapped = map_and_count_mapped(fastq, ref_path, preset, threads)
    pct = (mapped / total_reads) * 100
    return (os.path.basename(fastq), os.path.basename(ref_path), total_reads, mapped, pct)

def main():
    args = parse_args()

    fastqs = sorted(glob.glob(os.path.join(args.input, args.pattern)))
    if not fastqs:
        sys.exit(f" No files matching {args.pattern} found in {args.input}")

    # Get references
    if args.accession:
        refs = [download_reference(args.accession)]
    elif args.references:
        refs = args.references
    elif args.ref_dir:
        refs = sorted(glob.glob(os.path.join(args.ref_dir, "*.fna")))
        if not refs:
            sys.exit(f"No .fna files found in {args.ref_dir}")
    else:
        sys.exit("Must provide --references, --ref-dir, or --accession")

    if args.dry_run:
        print("Dry run: Files and references to be processed:")
        for fq in fastqs:
            for ref in refs:
                print(f"[DRY-RUN] {os.path.basename(fq)} vs {os.path.basename(ref)}")
        sys.exit(0)

    print(f"# Mapping {len(fastqs)} FASTQ files against {len(refs)} references")
    print(f"{'Sample':<25} {'Reference':<25} {'Total':>8} {'Mapped':>10} {'% Mapped':>9}")
    print("-" * 80)

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(process_one, fq, ref, args.preset, args.threads)
            for fq in fastqs for ref in refs
        ]
        for fut in as_completed(futures):
            sample, refname, total, mapped, pct = fut.result()
            print(f"{sample:<25} {refname:<25} {total:>8} {mapped:>10} {pct:9.2f}")
            results.append((sample, refname, total, mapped, pct))

    # Write TSV
    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["Sample", "Reference", "Total_Reads", "Mapped_Reads", "Percent_Mapped"])
        for row in results:
            writer.writerow(row)

    print(f"\nResults written to {args.output}")

if __name__ == "__main__":
    main()
