#!/usr/bin/env python3
"""
Create short reads from long on a single file
Response to reviewer who requested a script to do so 
Split a single FASTQ/FASTA file into fixed-length fragments.
By default uses 300 bp size, and **keeps** the final short fragment (< chunk_size).
"""

import os
import gzip
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _test_gzip(path: str) -> bool:
    """Return True if `path` is a valid gzip file, else False."""
    try:
        subprocess.run(
            ["gzip", "-t", path],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except subprocess.CalledProcessError:
        return False


def chunk_file(
    input_file: str,
    sample_name: str,
    chunk_size: int = 300,
    output_dir: str = ".",
) -> str:
    """
    Split input FASTQ/FASTA into chunks of size `chunk_size`.
    Keeps the final short fragment (< chunk_size).
    Writes gzipped output to:
        <output_dir>/<sample_name>_short_<chunk_size>.<fq.gz|fa.gz>
    """
    if input_file.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
        fmt, ext = "fastq", "fq.gz"
    elif input_file.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
        fmt, ext = "fasta", "fa.gz"
    else:
        raise ValueError(f"Unrecognized extension: {input_file!r}")

    os.makedirs(output_dir, exist_ok=True)
    out_name = os.path.join(output_dir, f"{sample_name}_short_{chunk_size}.{ext}")

    # If output exists and is valid, skip; if corrupt, remove and regenerate
    if os.path.exists(out_name):
        if _test_gzip(out_name):
            print(f"Skipping: {out_name} already exists and is valid.")
            return out_name
        print(f"Corrupted output found, removing: {out_name}")
        os.remove(out_name)

    # Write output (retry once if gzip test fails)
    for attempt in (1, 2):
        with gzip.open(out_name, "wt") as out_fh:
            opener = gzip.open if input_file.endswith(".gz") else open
            with opener(input_file, "rt") as in_fh:
                for rec in SeqIO.parse(in_fh, fmt):
                    seq = rec.seq
                    quals = rec.letter_annotations.get("phred_quality")

                    for i in range(0, len(seq), chunk_size):
                        chunk = seq[i : i + chunk_size]
                        if len(chunk) == 0:
                            continue  # defensive; should never happen

                        if fmt == "fastq":
                            if not quals:
                                continue
                            # Keep short fragments: slice qualities to match chunk length
                            chunk_quals = quals[i : i + len(chunk)]
                            if len(chunk_quals) != len(chunk):
                                continue  # defensive

                        new_id = f"{rec.id}_chunk{i}"
                        new_rec = SeqRecord(chunk, id=new_id, description="")

                        if fmt == "fastq":
                            new_rec.letter_annotations["phred_quality"] = chunk_quals

                        SeqIO.write(new_rec, out_fh, fmt)

        if _test_gzip(out_name):
            print(f"✓ Created {out_name}")
            return out_name

        print(f"Attempt {attempt} failed gzip test, retrying...")
        os.remove(out_name)

    raise IOError(f"Failed to write valid gzip for {out_name}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Split a single FASTQ/FASTA file into fixed-length fragments (keeps short tail). Allows single files to be split, rather than part of the entire sequence workflow which processes a directory and automates the naming proces "
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input FASTQ/FASTA file (optionally gzipped)."
    )
    parser.add_argument(
        "-s", "--sample", required=True,
        help="Will be the name of the output file."
    )
    parser.add_argument(
        "-c", "--chunk_size", type=int, default=300,
        help="Chunk size (default: 300)."
    )
    parser.add_argument(
        "-o", "--output", default=".",
        help="Output directory (default: current directory)."
    )

    args = parser.parse_args()

    chunk_file(
        input_file=args.input,
        sample_name=args.sample,
        chunk_size=args.chunk_size,
        output_dir=args.output,
    )


if __name__ == "__main__":
    main()
