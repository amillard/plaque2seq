#!/usr/bin/env python3

import os
import glob
import subprocess
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

# Default log file
LOG_FILE = "checkv_processing.log"

def log(message: str):
    """Append timestamped message to LOG_FILE and print it."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, "a") as f:
        f.write(f"[{ts}] {message}\n")
    print(message)

def run_cmd(cmd: str):
    """Run a shell command, logging stdout/stderr, and raising on failure."""
    log(f" Running: {cmd}")
    proc = subprocess.run(
        cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, executable="/bin/bash"
    )
    log(f" STDOUT:\n{proc.stdout}")
    log(f" STDERR:\n{proc.stderr}")
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, proc.stdout, proc.stderr)

def run_checkv_on_fasta(fasta: str, threads: int, db_dir: str) -> None:
    """
    Run CheckV end_to_end on a single contigs.fasta.
    """
    out_dir = os.path.join(os.path.dirname(fasta), "checkv")
    if os.path.isdir(out_dir):
        log(f" Skipping {fasta}: {out_dir!r} exists.")
        return

    cmd = f"checkv end_to_end {fasta} {out_dir} -t {threads} -d {db_dir}"
    log(f"  CheckV on {fasta}")
    try:
        run_cmd(cmd)
        tmp = os.path.join(out_dir, "tmp")
        if os.path.isdir(tmp):
            log(f"🧹 Removing tmp dir {tmp}")
            subprocess.run(f"rm -rf {tmp}", shell=True, check=True)
        log(f"✅ CheckV finished for {fasta}")
    except subprocess.CalledProcessError as e:
        log(f"❌ CheckV failed for {fasta} (rc={e.returncode})")

def find_contigs(base_dir: str) -> list:
    """
    Return all contigs.fasta paths within base_dir/*/contigs.fasta
    and base_dir/*/*/contigs.fasta
    standard way these contigs are stored
    """
    p1 = glob.glob(os.path.join(base_dir, "*", "contigs.fasta"))
    p2 = glob.glob(os.path.join(base_dir, "*", "*", "contigs.fasta"))
    return sorted(p1 + p2)

def process_directory(
    base_dir: str,
    threads: int = 16,
    workers: int = None,
    db_dir: str = "/scratch/entro/shared/checkv-db-v1.5",
    log_file: str = None
):
    """
    Find all contigs.fasta files and run CheckV in parallel.
    """
    global LOG_FILE
    if log_file:
        LOG_FILE = log_file

    log(f" Scanning {base_dir} for contigs.fasta")
    fastas = find_contigs(base_dir)
    log(f" Found {len(fastas)} files")

    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_fasta = {
            executor.submit(run_checkv_on_fasta, f, threads, db_dir): f for f in fastas
        }

        for future in as_completed(future_to_fasta):
            fasta = future_to_fasta[future]
            try:
                future.result()
            except Exception as e:
                log(f"❌ Exception for {fasta}: {e}")

# === CLI ===
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run CheckV on all contigs.fasta files in a directory tree.")
    parser.add_argument("base_dir", help="Directory under which to search for contigs.fasta")
    parser.add_argument("-t", "--threads", type=int, default=16,
                        help="Threads per CheckV run [default: 16]")
    parser.add_argument("-w", "--workers", type=int, default=None,
                        help="Max parallel CheckV jobs [default: all cores]")
    parser.add_argument("-d", "--db_dir", default="/scratch/entro/shared/checkv-db-v1.5",
                        help="Path to CheckV database directory")
    parser.add_argument("--log_file", default=None,
                        help="Path to log file (default: checkv_processing.log)")

    args = parser.parse_args()

    process_directory(
        base_dir=args.base_dir,
        threads=args.threads,
        workers=args.workers,
        db_dir=args.db_dir,
        log_file=args.log_file,
    )
