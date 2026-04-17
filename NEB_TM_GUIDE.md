# NEB Tm API Usage Guide

This guide explains how to use the `neb_tm_calculator.py` script to automate Tm and Ta calculations for primers designed with the CRISPR Primer Designer.

## Prerequisites

Ensure you are in the project folder and have the virtual environment activated:

```bash
source .venv/bin/activate
```

## Running the Calculation

After running the main primer design script, you will have a `primers_primers.tsv` (or `.csv`) file. To calculate NEB-specific temperatures, run:

```bash
python neb_tm_calculator.py primers_primers.tsv
```

### Options

- **Custom Output Path**:
  ```bash
  python neb_tm_calculator.py primers_primers.tsv -o my_results_with_tm.tsv
  ```
- **Custom Email**:
  ```bash
  python neb_tm_calculator.py primers_primers.tsv -e your_email@gene.com
  ```
- **Skip SSL Verification** (if you have network issues):
  ```bash
  python neb_tm_calculator.py primers_primers.tsv --no-verify
  ```
- **Generate Batch Processing File** (to manually upload to NEB website):
  ```bash
  python neb_tm_calculator.py primers_primers.tsv --batch-file neb_batch.csv
  ```
  This creates a CSV that can be uploaded to [NEB Tm Calculator Batch Mode](https://tmcalculator.neb.com/#!/batch).

## Output Format

The script creates a new file (default suffix `_with_neb`) with three additional columns:

1.  **NEB_Tm1**: Melting temperature for the left primer.
2.  **NEB_Tm2**: Melting temperature for the right primer.
3.  **NEB_Ta**: Recommended annealing temperature for **Q5 Hot Start 2X Master Mix**.

---
*Note: All calculations use the 'q5-0' product code on the [NEB Tm API](https://tmapi.neb.com/).*
