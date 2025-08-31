# extract viral/phage contigs from assemblies
import os
import pandas as pd
from Bio import SeqIO



"""
Module: was extract_viral.py

Reads a combined quality summary TSV (no header) and extracts contigs
where viral_genes > host_genes from each assembly’s contigs.fasta.
Writes each selected contig to its own FASTA in an output directory.
"""


def extract_contigs(
    summary_path: str,
    output_base_dir: str = "extracted_contigs",
    dry_run: bool = False
) -> None:
    """
    Parse the summary TSV and for each row where viral_genes > host_genes,
    extract the matching contig from dir/contigs.fasta.

    summary_path:       path to combined_quality_summary.tsv
    output_base_dir:    where to write per-contig FASTAs
    dry_run:            if True, only print what would be done
    """
    # Define the expected columns (no header in TSV)
    columns = [
        "dir", "contig_id", "contig_length", "provirus", "proviral_length",
        "gene_count", "viral_genes", "host_genes", "checkv_quality",
        "miuvig_quality", "completeness", "completeness_method",
        "contamination", "kmer_freq", "warnings"
    ]

    # Load summary
    df = pd.read_csv(summary_path, sep="\t", header=None, names=columns, dtype=str)
    # Clean and convert
    df["dir"] = df["dir"].str.strip()
    df["contig_id"] = df["contig_id"].str.strip()
    df["viral_genes"] = pd.to_numeric(df["viral_genes"], errors="coerce").fillna(0).astype(int)
    df["host_genes"] = pd.to_numeric(df["host_genes"], errors="coerce").fillna(0).astype(int)

    # Filter for more viral genes than host genes
    df = df[df["viral_genes"] > df["host_genes"]].copy()
    if df.empty:
        print("[INFO] No contigs pass viral > host filter.")
        return

    # Assign group IDs per (dir, checkv_quality)
    df["group_id"] = df.groupby(["dir", "checkv_quality"]).cumcount() + 1

    # Prepare output directory
    if not dry_run:
        os.makedirs(output_base_dir, exist_ok=True)

    # Iterate and extract
    for idx, row in df.iterrows():
        dir_path = row["dir"]
        contig_id = row["contig_id"]
        quality = row["checkv_quality"]
        group_id = row["group_id"]

        sample_name = os.path.basename(os.path.normpath(dir_path))
        out_fname = f"{sample_name}_{quality}_{group_id}.fasta"
        out_path = os.path.join(output_base_dir, out_fname)

        print(f"[INFO] Row {idx}: sample={sample_name}, contig={contig_id}, "
              f"quality={quality}, group={group_id}")

        fasta_in = os.path.join(dir_path, "contigs.fasta")
        if not os.path.isfile(fasta_in):
            print(f"[ERROR] Missing contigs.fasta in {dir_path}; skipping.")
            continue

        if dry_run:
            print(f"  [DRY RUN] Would extract {contig_id} → {out_path}")
            continue

        # Perform extraction
        found = False
        try:
            with open(fasta_in) as hin:
                for rec in SeqIO.parse(hin, "fasta"):
                    if rec.id == contig_id:
                        # rename record to the output basename (no suffix)
                        new_id = os.path.splitext(out_fname)[0]
                        rec.id = rec.description = new_id
                        with open(out_path, "w") as hout:
                            SeqIO.write(rec, hout, "fasta")
                        print(f"  [OK] Wrote {contig_id} → {out_path}")
                        found = True
                        break
            if not found:
                print(f"  [WARN] Contig {contig_id} not found in {fasta_in}")
        except Exception as e:
            print(f"  [ERROR] Failed to extract from {fasta_in}: {e}")

if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Extract viral contigs (viral_genes > host_genes) from assemblies"
    )
    p.add_argument("summary", help="Path to combined_quality_summary.tsv")
    p.add_argument(
        "--outdir", "-o", default="extracted_contigs",
        help="Directory to write extracted FASTAs"
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="Only print actions, do not write files"
    )
    args = p.parse_args()

    extract_contigs(
        summary_path=args.summary,
        output_base_dir=args.outdir,
        dry_run=args.dry_run
    )
