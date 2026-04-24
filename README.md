# CRISPR Cas9 Primer Designer

This repository provides tools to design optimized primers for sequencing CRISPR Cas9 cut sites, specifically tailored for Illumina sequencing. It automates the primer design process using Primer3, handling matches between target guides and larger genomic template sequences.

## Features
- **Primer3 Integration**: Uses the powerful Primer3 engine for optimized primer selection.
- **Sequence Matching**: Automatically locates gRNA targets within large template sequences.
- **Multi-Sequence Support**: Scans multiple sequences in a single FASTA file to find all occurrences of a cut site.
- **Direct gRNA Input**: Option to provide a gRNA sequence directly via command line without a TSV.
- **Thermodynamic Optimization**: Optimized for Illumina sequencing (MiSeq/NextSeq) with flexible parameters for Tm, GC content, and product size.
- **Multiple Result Handling**: Returns multiple ranked primer pairs per target.

## NEB Tm Calculation

After designing primers, you can automatically calculate optimal annealing temperatures (Ta) and melting temperatures (Tm) for **Q5 Hot Start 2X Master Mix** using the NEB Tm API.

```bash
python neb_tm_calculator.py path/to/your_primers.tsv
```

See [NEB_TM_GUIDE.md](NEB_TM_GUIDE.md) for more detailed options.

## Prerequisites

### Tools
- **Samtools** (Optional): Useful for extracting genomic sequences if you are preparing your own templates from a reference genome.
  ```bash
  brew install samtools
  ```
- **Primer3**: The core engine for primer design. The script attempts to find `primer3_core` in common locations, including `/usr/local/bin/` and within the repository's `primer3/src/` directory.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/robertmunch/CRISPR_primer_designer.git
   cd CRISPR_primer_designer
   ```

2. Install the necessary Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Designing Primers

The main tool is `crispr_primer_designer.py`. You can provide template sequences in two ways:

#### Option A: Automatic Genomic Extraction (Recommended)
Provide a TSV with genomic coordinates and a reference genome FASTA. The script will automatically extract the surrounding region for each target.

**Don't have hg38?** Use `--download-genome` to fetch the **GRCh38 no-ALT analysis set** automatically from NCBI (~833 MB compressed). This version excludes ALT contigs (`chr*_alt`) so that primer uniqueness checks reflect the primary assembly only — sequences that are unique in the primary assembly will not incorrectly fail because they also appear in an ALT contig. The file is cached at `~/.cache/crispr_primer_designer/hg38_no_alt.fa` after the first download. **It is safe to pass `--download-genome` on every run** — if the file is already cached the download is skipped instantly.

```bash
python crispr_primer_designer.py \
  -c cutsites.tsv \
  --download-genome \
  -f 250 \
  -o my_primers
```

If you already have a local copy of hg38 somewhere else, point to it with `-g` instead:

```bash
python crispr_primer_designer.py \
  -c cutsites.tsv \
  -g /path/to/hg38.fa \
  -f 250 \
  -o my_primers
```

#### Option B: Manual FASTA Templates
Provide a TSV with guide information and a pre-prepared FASTA file containing the template sequences (can contain multiple sequences).

```bash
python crispr_primer_designer.py \
  -c cutsites.tsv \
  -t templates.fasta \
  -o my_primers
```

#### Option C: Direct gRNA Input
Provide a gRNA sequence directly via the command line and match it against a multi-sequence FASTA.

```bash
python crispr_primer_designer.py \
  --grna GCTGACTTACCAGATGGGAC \
  -t templates.fasta \
  -o my_primers
```

### Input File Formats

#### 1. Cutsites TSV (`-c`)
A tab-separated file with at least these columns:
- `target_gene`: A label for the target.
- `grna_name`: Unique name for the guide.
- `sequence`: The gRNA sequence (e.g., 20bp).
- `location`: Genomic coordinates (required for Option A) or a simple position (e.g., `chr1:1000-2000` or `1500`).

#### 2. Genome FASTA (`-g`) / Templates FASTA (`-t`)
- **Genome FASTA**: Used with Option A. Must be a indexed FASTA file (e.g., `hg38.fa` with `hg38.fa.fai`).
- **Templates FASTA**: Used with Option B. Standard FASTA file containing DNA regions.

### Command Line Arguments
- `-c`, `--cutsites`: Path to the TSV file with target information.
- `--grna`: Directly provide a gRNA sequence (skips the need for `-c`).
- `--grna-name`: Name/Label for the provided `--grna` (default: "gRNA").
- `-g`, `--genome`: Path to a local genome FASTA for automatic extraction.
- `--download-genome`: Automatically download hg38 from UCSC (cached after first use). Mutually exclusive with `-g`.
- `-f`, `--flank`: Flanking distance around cut site (default: 250bp).
- `-t`, `--templates`: Path to the FASTA file with template sequences (can contain multiple records).
- `--target-distance`: Minimum distance from cut site to primer (default: 50bp).
- `-o`, `--prefix`: Prefix for the output files (default: `primers`).
- `--output_dir`: Directory to store the output results.
- `--check-uniqueness`: (Optional) Verify that primers are unique in the genome (requires `-g`).

## Output Details

The tool generates several files:
- `[prefix]_primers.tsv`: The main results file containing primer sequences, properties (Tm, GC), and amplicon details.
- `[prefix]_snapgene.tsv`: A simplified tab-delimited list (`Name` [tab] `Sequence`) specifically formatted for easy import into SnapGene.
- `[prefix]_bulk_order.tsv`: A TSV formatted for oligo ordering (includes TruSeq adapters, 25nm scale, and STD purification).
- `[prefix]_input.txt`: The raw input generated for Primer3.
- `[prefix]_output.txt`: The raw output from the Primer3 run.

## Troubleshooting

- **No Primers Found**: Try relaxing melting temperature (Tm) constraints in `crispr_primer_designer.py` or providing larger template sequences.
- **Sequence Not Found**: Ensure the gRNA sequences in your `-c` file exist within the sequences provided in your `-t` file.
- **Primer3 Executable Not Found**: Ensure `primer3_core` is in your system PATH or located in the project's `primer3/src/` folder.
