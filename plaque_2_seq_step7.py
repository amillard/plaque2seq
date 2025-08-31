#!/usr/bin/env python3
"""
Script: extend_prefix_match.py

Detects and trims circular overlap in assembled contigs by matching
prefixes to suffixes. Produces trimmed FASTA and logs results.
SPAdes is great but it puts repeat sequences at the end of contigs that are 
circular, this is by design. Needs to be trimmed prior to re-ordering with 
dnaapler

Usage:
    python extend_prefix_match.py --indir assemblies/ --outdir trimmed/
"""

import sys
import os
import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def process_file(fasta_file, outdir, log_handle):
    """Process one FASTA: find longest prefix-suffix match, trim, write outputs, log result."""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        print(f"{fasta_file}: no records found, skipping.")
        return

    record = records[0]
    seq_str = str(record.seq)
    n = len(seq_str)

    if n < 50:
        print(f"{fasta_file}: sequence <50 bp, skipping.")
        return

    # --- Trim ass_ from filename ---
    base_name = os.path.basename(os.path.splitext(fasta_file)[0])
    if base_name.startswith("ass_"):
        base_name = base_name[4:]

    trimmed_file = os.path.join(outdir, base_name + '.fa')

    if os.path.exists(trimmed_file):
        print(f"{fasta_file}: already processed, skipping.")
        return

    # Search for longest prefix match in last 200 bp
    last200 = seq_str[-200:] if n >= 200 else seq_str
    max_len = 0
    match_start = match_end = None

    for k in range(50, n + 1):
        prefix = seq_str[:k]
        idx = last200.find(prefix)
        if idx != -1:
            max_len = k
            offset = n - len(last200)
            match_start = offset + idx + 1
            match_end = match_start + k - 1
        else:
            break

    trimmed_length = max_len if max_len >= 50 else 0
    trimmed_seq_str = seq_str[:-trimmed_length] if trimmed_length > 0 else seq_str

    # --- Trim ass_ from header ---
    header_id = record.id
    if header_id.startswith("ass_"):
        header_id = header_id[4:]

    # Write trimmed FASTA
    trimmed_record = SeqRecord(
        Seq(trimmed_seq_str),
        id=header_id,
        description=header_id
    )
    SeqIO.write(trimmed_record, trimmed_file, "fasta")

    # Log entry
    match_pos = f"{match_start}-{match_end}" if trimmed_length > 0 else ''
    log_handle.write(f"{fasta_file}\t{('Yes' if trimmed_length > 0 else 'No')}\t{trimmed_length}\t{match_pos}\n")
    log_handle.flush()

    print(f"Processed {fasta_file}: trimmed_file={trimmed_file}")

def collect_fasta_files(indir):
    """Return all *.fasta and *.fa files in a directory."""
    fasta_patterns = [os.path.join(indir, "*.fasta"), os.path.join(indir, "*.fa")]
    files = []
    for pattern in fasta_patterns:
        files.extend(glob.glob(pattern))
    return files

def main():
    parser = argparse.ArgumentParser(
        description="Detect and trim prefix-suffix matches in assembled circular contigs."
    )
    parser.add_argument("fastas", nargs="*", help="FASTA files to process (optional if using --indir)")
    parser.add_argument("--indir", help="Directory containing FASTA files (*.fasta or *.fa)")
    parser.add_argument("--outdir", "-o", default="trimmed", help="Directory to write trimmed FASTAs")
    parser.add_argument("--log", default="combined.log", help="Path to log file")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Collect all input FASTAs
    input_files = set(args.fastas)
    if args.indir:
        input_files.update(collect_fasta_files(args.indir))

    input_files = sorted(input_files)
    if not input_files:
        print("❌ No input FASTA files provided or found.")
        sys.exit(1)

    # Open log file once for writing
    with open(args.log, 'w') as log_handle:
        log_handle.write("Input_file\tTrimmed\tLength_trimmed\tMatch_positions\n")

        for fasta in input_files:
            process_file(fasta, args.outdir, log_handle)

if __name__ == "__main__":
    main()
