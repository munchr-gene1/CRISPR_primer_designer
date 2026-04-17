#!/usr/bin/env python3
import csv
import json
import requests
import argparse
import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple

"""
NEB Tm Calculator for CRISPR Primer Designer
Uses the NEB Tm API (https://tmapi.neb.com/) to calculate Tm and Ta values for primers.
Target Polymerase: Q5 Hot Start 2X Master Mix (code: q5)
"""

def get_neb_data(seq_pairs: List[Tuple[str, str]], email: str = "tmapi@neb.com", prod_code: str = "q5-0", conc: float = 0.5, verify_ssl: bool = True):
    """
    Fetch Tm and Ta values from NEB API for a batch of sequence pairs.
    """
    url = "https://tmapi.neb.com/tm/batch"
    payload = {
        "prodcode": prod_code,
        "conc": conc,
        "email": email,
        "seqpairs": seq_pairs
    }
    
    try:
        response = requests.post(url, json=payload, timeout=30, verify=verify_ssl)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error calling NEB API: {e}")
        return None

def process_primer_file(input_file: str, output_file: str = None, email: str = "tmapi@neb.com", verify_ssl: bool = True):
    """
    Read primer file, fetch NEB data, and save updated file.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        print(f"Error: File {input_file} not found.")
        return

    # Determine delimiter
    delimiter = '\t' if input_path.suffix.lower() == '.tsv' else ','
    
    # Read rows
    rows = []
    with open(input_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)
            
    if not rows:
        print("No data found in file.")
        return

    # Extract sequences
    # Assuming column names 'Left_Primer' and 'Right_Primer' based on crispr_primer_designer.py
    seq_pairs = []
    for row in rows:
        left = row.get('Left_Primer')
        right = row.get('Right_Primer')
        if left and right:
            seq_pairs.append([left, right])
        else:
            seq_pairs.append([None, None])

    # Filter out empty pairs for API call
    valid_pairs = [p for p in seq_pairs if p[0] and p[1]]
    
    if not valid_pairs:
        print("No valid primer pairs found to calculate.")
        return

    print(f"Fetching data for {len(valid_pairs)} primer pairs from NEB...")
    api_data = get_neb_data(valid_pairs, email=email, verify_ssl=verify_ssl)
    
    if not api_data or 'data' not in api_data:
        print("Failed to retrieve data from NEB API.")
        return

    # Map API results back to rows
    # The API returns data in the same order as seqpairs
    results_map = {}
    for i, pair in enumerate(valid_pairs):
        res = api_data['data'][i]
        results_map[tuple(pair)] = res

    # Add new fieldnames
    new_fields = ['NEB_Tm1', 'NEB_Tm2', 'NEB_Ta']
    output_fieldnames = fieldnames + [f for f in new_fields if f not in fieldnames]

    # Prepare output file path
    if not output_file:
        output_file = str(input_path.with_name(f"{input_path.stem}_with_neb{input_path.suffix}"))

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames, delimiter=delimiter)
        writer.writeheader()
        
        for row in rows:
            left = row.get('Left_Primer')
            right = row.get('Right_Primer')
            pair = tuple([left, right])
            
            if pair in results_map:
                res = results_map[pair]
                row['NEB_Tm1'] = f"{res.get('tm1', 0):.1f}"
                row['NEB_Tm2'] = f"{res.get('tm2', 0):.1f}"
                row['NEB_Ta'] = f"{res.get('ta', 0):.1f}"
            else:
                row['NEB_Tm1'] = 'N/A'
                row['NEB_Tm2'] = 'N/A'
                row['NEB_Ta'] = 'N/A'
            
            writer.writerow(row)

    print(f"Success! Updated file saved to: {output_file}")

def generate_neb_batch_file(input_file: str, batch_file: str):
    """
    Generate a CSV file for NEB Tm calculator batch processing website.
    Format: Name1, Seq1, Name2, Seq2
    """
    input_path = Path(input_file)
    if not input_path.exists():
        print(f"Error: File {input_file} not found.")
        return

    # Determine delimiter
    delimiter = '\t' if input_path.suffix.lower() == '.tsv' else ','
    
    # Read rows
    rows = []
    with open(input_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            rows.append(row)
            
    if not rows:
        print("No data found in file.")
        return

    with open(batch_file, 'w', newline='') as f:
        writer = csv.writer(f)
        for i, row in enumerate(rows):
            left = row.get('Left_Primer')
            right = row.get('Right_Primer')
            if not left or not right:
                continue
                
            # Try to get descriptive names from Sequence_ID and Pair_Index
            seq_id = row.get('Sequence_ID', f'P{i+1}')
            pair_idx = row.get('Pair_Index')
            
            if pair_idx is not None and pair_idx != '':
                name1 = f"{seq_id}_{pair_idx}_fwd"
                name2 = f"{seq_id}_{pair_idx}_rev"
            else:
                name1 = f"{seq_id}_fwd"
                name2 = f"{seq_id}_rev"
            
            writer.writerow([name1, left, name2, right])

    print(f"Success! NEB batch processing file saved to: {batch_file}")
    print("You can now upload this file to https://tmcalculator.neb.com/#!/batch")

def main():
    parser = argparse.ArgumentParser(description="Calculate NEB Tm and Ta values for primers.")
    parser.add_argument("input_file", help="Path to the primers_primers.tsv or .csv file.")
    parser.add_argument("-o", "--output", help="Optional output file path.")
    parser.add_argument("-e", "--email", default="tmapi@neb.com", help="Email for NEB API usage.")
    parser.add_argument("--no-verify", action="store_true", help="Disable SSL certificate verification.")
    parser.add_argument("--batch-file", help="Generate a CSV file for NEB batch processing and exit.")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    
    if args.batch_file:
        generate_neb_batch_file(args.input_file, args.batch_file)
    else:
        process_primer_file(args.input_file, args.output, args.email, verify_ssl=not args.no_verify)

if __name__ == "__main__":
    main()
