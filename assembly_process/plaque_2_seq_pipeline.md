# Plaque2Seq: A Phage Genome Assembly and Curation Pipeline

This pipeline processes full-length sequencing reads into curated viral contigs through a series of modular steps. Each script handles one part of the process, from read fragmentation to CheckV quality control and Pilon polishing.



can recreate the environment with the yml file 



## Step-by-step Usage

### **Step 1** – Split Full Reads into 300 bp reads

```bash
python plaque_2_seq_1.py -d full_reads/ -m bc2sample.tsv -c 300 -w 8
```

- `-d`: Input directory with full-length reads
- `-m`: Barcode-to-sample mapping file
- `-c`: read size (default: 300)
- `-w`: Number of workers (parallel jobs) speeds up the process 

---

### **Step 2** – Normalize Reads to Target Coverage

```bash
python plaque_2_seq_seq_step2.py -t 150 -o ./ ./
```

- `-t`: Target coverage (default: 150)
- `-o`: Output directory
  
  run in the current directory is the easiest thing to 

---

### **Step 3** – Assemble Normalized Reads with SPAdes

```bash
python plaque_2_seq_step3.py -d ./ -t 4 -w 4
```

- `-d`: Directory with `*_norm.fq.gz` files
- `-t`: Threads per assembly job
- `-w`: Number of concurrent assemblies 

less cores and multiple concurrent runs is fastest 4 * 4 is quicker than 16 threads 



---

### **Step 4** – Run CheckV on All Assemblies

```bash
python plaque_2_seq_step4.py -t 4 -w 4 -d ~/shared/checkv-db-v1.5/ ./
```

- `-t`: Threads per CheckV run
- `-w`: Max parallel CheckV jobs
- `-d`: Path to CheckV database
- `./`: Base directory containing `contigs.fasta`

run checkv 

---

### **Step 5** – Concatenate All CheckV Summaries

```bash
python plaque_2_seq_step5.py
```

No arguments needed; gathers `quality_summary.tsv` files.

---

### **Step 6** – Extract Viral Contigs

```bash
python plaque_2_seq_step6a.py --dry-run combined_quality_summary.tsv
```

- `--outdir`: Output directory (optional)
- `--dry-run`: Print what would be extracted, without writing files

To extract contigs:

```bash
python extract_viral.py combined_quality_summary.tsv
```

---

### **Step 7** – Detect and Trim Circular Contig Overlaps

```bash
python plaque_2_seq_step7.py --indir clean/ --outdir trimmed/ --log trim.log
```

- `--indir`: Directory with input FASTA files
- `--outdir`: Directory for trimmed contigs
- `--log`: Log file path

---

### **Step 8** – Map Reads to Viral Contigs and Compute Coverage

```bash
python plaque_2_seq_step8.py \
  --dict bc2sample.tsv \
  --clean_dir clean/ \
  --reads_dir ./ \
  --bam_dir bam/ \
  -w 4
```

- `--dict`: Sample mapping file
- `--clean_dir`: Directory with viral contigs
- `--reads_dir`: Directory with chunked reads
- `--bam_dir`: Output directory for BAMs
- `-w`: Number of workers

---

### **Step 9** – Run Pilon for Assembly Polishing

```bash
python plaque_2_seq_step9.py \
  --dict bc2sample.tsv \
  --clean_dir clean/ \
  --bam_dir bam/ \
  --outdir pilon/
```

- `--dict`: Sample mapping file
- `--clean_dir`: Polished FASTA files
- `--bam_dir`: BAM alignment files
- `--outdir`: Directory for Pilon output

---

## Notes

- Ensure all tools are available in your environment (SPAdes, CheckV, Pilon, etc.)
- Use `--help` for more information on each step.
- Scripts are modular and can be rerun independently. 
- Scripts generally check if the file exists already 
