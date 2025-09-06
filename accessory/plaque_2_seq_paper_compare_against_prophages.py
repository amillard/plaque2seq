#!/usr/bin/env python3

import argparse
import subprocess
import os
from glob import glob
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

def run_minimap2(query_fasta, target_fasta):
    """Run minimap2 and parse output."""
    try:
        result = subprocess.run(
            ["minimap2", "-c", "--score-N=0", target_fasta, query_fasta],
            check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        print("Error running minimap2:", e.stderr)
        return []

    alignments = []
    for line in result.stdout.strip().split("\n"):
        if line.startswith("@") or not line.strip():
            continue
        fields = line.split("\t")
        qname, tname = fields[0], fields[5]
        aln_len = int(fields[10])
        n_matches = int(fields[9])
        pct_id = 100 * n_matches / aln_len if aln_len > 0 else 0
        alignments.append((qname, tname, pct_id, aln_len))
    return alignments

def get_lengths(fasta_path):
    """Return a dict of sequence lengths."""
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}

def process_file(input_fasta, prophage_fasta, prophage_ids, prophage_lengths):
    query_lengths = get_lengths(input_fasta)
    query_ids = list(query_lengths.keys())

    data = defaultdict(lambda: {"Query_length": 0})
    query_hits = defaultdict(list)

    for qid in query_ids:
        data[qid]["Query_length"] = query_lengths[qid]
        for pid in prophage_ids:
            data[qid][f"%_sim_to_{pid}"] = 0.0
            data[qid][f"aln_len_to_{pid}"] = 0

    alignments = run_minimap2(input_fasta, prophage_fasta)
    for qname, tname, pct_id, aln_len in alignments:
        data[qname][f"%_sim_to_{tname}"] = round(pct_id, 2)
        data[qname][f"aln_len_to_{tname}"] = aln_len
        query_hits[qname].append((pct_id, aln_len))

    rows = []
    for qid in query_ids:
        query_len = data[qid]["Query_length"]
        row = {
            "File": os.path.basename(input_fasta),
            "Query": qid,
            "Query_length": query_len
        }

        for pid in prophage_ids:
            row[f"%_sim_to_{pid}"] = data[qid][f"%_sim_to_{pid}"]
            row[f"aln_len_to_{pid}"] = data[qid][f"aln_len_to_{pid}"]

        # Updated classification logic: alignment must be > 70% of query length
        row["Classification"] = "Phage"
        for pct_id, aln_len in query_hits[qid]:
            if pct_id >= 98.0 and aln_len > 0.7 * query_len:
                row["Classification"] = "Host"
                break

        rows.append(row)

    return rows




def main():
    parser = argparse.ArgumentParser(description="Compare all FASTA files in a directory to prophage reference")
    parser.add_argument("-d", "--indir", required=True, help="Directory containing FASTA files")
    parser.add_argument("-p", "--prophages", default="LEI_0034_prophages.fna", help="Prophage reference file")
    parser.add_argument("-o", "--out", default="prophage_alignment_summary.tsv", help="Output TSV file name")
    args = parser.parse_args()

    # Load prophage data
    prophage_lengths = get_lengths(args.prophages)
    prophage_ids = list(prophage_lengths.keys())

    all_rows = []
    fasta_files = glob(os.path.join(args.indir, "*.fasta")) + glob(os.path.join(args.indir, "*.fa")) + glob(os.path.join(args.indir, "*.fna"))
    if not fasta_files:
        print(f"No FASTA files found in {args.indir}")
        return

    for fasta in sorted(fasta_files):
        rows = process_file(fasta, args.prophages, prophage_ids, prophage_lengths)
        all_rows.extend(rows)

    # Final output
    all_columns = ["File", "Query", "Query_length"]
    for pid in prophage_ids:
        all_columns += [f"%_sim_to_{pid}", f"aln_len_to_{pid}"]
    all_columns.append("Classification")

    df = pd.DataFrame(all_rows)[all_columns]

    # Save and print
    df.to_csv(args.out, sep="\t", index=False)
    print(df.to_csv(sep="\t", index=False))

if __name__ == "__main__":
    main()
