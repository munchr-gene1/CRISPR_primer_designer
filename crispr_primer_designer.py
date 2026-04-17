#!/usr/bin/env python3
"""
CRISPR Cas9 Primer Designer using Primer3

This script designs primers to sequence CRISPR Cas9 cut sites using Primer3.
It takes a TSV file with gRNA information and large template sequences,
then generates primers optimized for Illumina sequencing.

"""

import argparse
import csv
import subprocess
import os
import re
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple
import json
import multiprocessing
import time
from functools import partial


@dataclass
class CutSite:
    """Represents a CRISPR Cas9 cut site"""
    target_gene: str
    grna_name: str
    sequence: str
    location: str  # Can be genomic coordinates or position
    extra_data: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        # Clean the sequence (remove spaces, uppercase)
        self.sequence = self.sequence.strip().upper().replace(" ", "")


@dataclass
class TemplateSequence:
    """Represents a large template sequence region"""
    name: str
    sequence: str
    chromosome: Optional[str] = None
    start_position: Optional[int] = None
    end_position: Optional[int] = None


@dataclass
class Primer3Config:
    """Configuration for Primer3 parameters"""
    # Product size range - optimized for Illumina sequencing
    product_size_range: str = "200-450"  # Illumina MiSeq/NextSeq optimal range

    # Primer size parameters
    primer_min_size: int = 18
    primer_opt_size: int = 20
    primer_max_size: int = 25

    # Melting temperature parameters
    primer_min_tm: float = 57.0
    primer_opt_tm: float = 60.0
    primer_max_tm: float = 63.0
    primer_pair_max_diff_tm: float = 3.0

    # GC content
    primer_min_gc: float = 40.0
    primer_opt_gc: float = 50.0
    primer_max_gc: float = 60.0

    # Number of primers to return
    num_return: int = 5

    # Target flanking - minimum distance from target to primer
    # This ensures the cut site is well-centered in the amplicon
    target_flank_min: int = 50

    # Self-complementarity limits
    primer_max_self_any: float = 8.0
    primer_max_self_end: float = 3.0
    primer_pair_max_compl_any: float = 8.0
    primer_pair_max_compl_end: float = 3.0

    # Thermodynamic parameters
    primer_salt_monovalent: float = 50.0
    primer_salt_divalent: float = 1.5
    primer_dntp_conc: float = 0.6
    primer_dna_conc: float = 50.0

    # Use thermodynamic alignment
    primer_thermodynamic_oligo_alignment: int = 1

    # Additional parameters
    primer_max_ns_accepted: int = 0
    primer_max_poly_x: int = 3
    primer_gc_clamp: int = 1
    primer_max_end_gc: int = 3

    # Explain flag for debugging
    primer_explain_flag: int = 1

    # Additional custom parameters
    extra_params: Dict[str, str] = field(default_factory=dict)


class Primer3InputGenerator:
    """Generates Primer3 input files from cut site and template data"""

    def __init__(self, config: Primer3Config):
        self.config = config

    def find_sequence_in_template(self, query: str, template: str,
                                   allow_reverse_complement: bool = True) -> List[Tuple[int, int, str]]:
        """
        Find a sequence within a template.
        Returns list of (start_position, length, strand) tuples.
        """
        matches = []
        if not query:
            return []
        
        query = query.upper()
        template = template.upper()

        # Search forward strand
        start = 0
        while True:
            pos = template.find(query, start)
            if pos == -1:
                break
            matches.append((pos, len(query), '+'))
            start = pos + 1

        # Search reverse complement if allowed
        if allow_reverse_complement:
            rev_comp = self._reverse_complement(query)
            start = 0
            while True:
                pos = template.find(rev_comp, start)
                if pos == -1:
                    break
                matches.append((pos, len(rev_comp), '-'))
                start = pos + 1

        return matches

    def _reverse_complement(self, seq: str) -> str:
        """Generate reverse complement of a DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                      'N': 'N', 'R': 'Y', 'Y': 'R', 'M': 'K',
                      'K': 'M', 'S': 'S', 'W': 'W'}
        return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))

    def parse_location(self, location_str: str) -> Tuple[Optional[str], Optional[int], Optional[int]]:
        """
        Parse location string to extract chromosome, start, end.
        Handles formats like:
        - "chr1:1000-2000"
        - "1000-2000"
        - "1000,20" (start,length format)
        - "1000" (just position)
        """
        location_str = location_str.strip()

        # Try chr:start-end format
        match = re.match(r'(chr\w+):(\d+)-(\d+)', location_str, re.IGNORECASE)
        if match:
            return match.group(1), int(match.group(2)), int(match.group(3))

        # Try start-end format
        match = re.match(r'(\d+)-(\d+)', location_str)
        if match:
            return None, int(match.group(1)), int(match.group(2))

        # Try start,length format (Primer3 style)
        match = re.match(r'(\d+),(\d+)', location_str)
        if match:
            start = int(match.group(1))
            length = int(match.group(2))
            return None, start, start + length

        # Try just position
        match = re.match(r'(\d+)', location_str)
        if match:
            pos = int(match.group(1))
            return None, pos, pos + 1

        return None, None, None

    def calculate_target(self, cut_site: CutSite, template: TemplateSequence,
                         target_padding: int = 10) -> Optional[Tuple[int, int]]:
        """
        Calculate SEQUENCE_TARGET position for a cut site within a template.
        """
        target = None

        # 1. Try to find the gRNA sequence in the template
        matches = self.find_sequence_in_template(cut_site.sequence, template.sequence)

        if matches:
            match_start, match_len, strand = matches[0]
            if strand == '+':
                cut_position = match_start + match_len - 3
            else:
                cut_position = match_start + 3

            target_start = max(0, cut_position - target_padding)
            target_length = 2 * target_padding + 6
            target = (target_start, target_length)

        # 2. If not found, try to use location information
        if target is None:
            chrom, start, end = self.parse_location(cut_site.location)
            if start is not None:
                if template.start_position is not None:
                    relative_start = start - template.start_position
                    if 0 <= relative_start < len(template.sequence):
                        target_length = (end - start) if end else 20
                        target = (relative_start, target_length + 2 * target_padding)
                else:
                    if 0 <= start < len(template.sequence):
                        target_length = (end - start) if end else 20
                        target = (start, target_length + 2 * target_padding)

        if target:
            return target

        # 3. Fallback: If it's a small template and we can't find the sequence,
        # assume the target is in the center. This is ideal for cases where the
        # target sequence is expected to be centered in the provided template.
        if len(template.sequence) < 2000:
            center = len(template.sequence) // 2
            return (max(0, center - target_padding), 2 * target_padding + 6)

        return None

    def generate_primer3_input(self, cut_site: CutSite, template: TemplateSequence,
                                sequence_id: Optional[str] = None) -> Optional[str]:
        """
        Generate a Primer3 input record for a single cut site.

        Returns the Boulder-IO formatted input string, or None if target couldn't be located.
        """
        # Calculate target position
        target = self.calculate_target(cut_site, template, target_padding=self.config.target_flank_min)

        if target is None:
            print(f"Warning: Could not locate cut site {cut_site.grna_name} in template {template.name}")
            return None

        target_start, target_length = target

        # Ensure target is within bounds
        if target_start < 0:
            target_start = 0
        if target_start + target_length > len(template.sequence):
            target_length = len(template.sequence) - target_start

        # Generate sequence ID
        if sequence_id is None:
            clean_template_name = template.name
            if clean_template_name.endswith('_template'):
                clean_template_name = clean_template_name[:-9]

            if clean_template_name == cut_site.grna_name:
                sequence_id = clean_template_name
            elif cut_site.grna_name in clean_template_name:
                sequence_id = clean_template_name
            elif clean_template_name in cut_site.grna_name:
                sequence_id = cut_site.grna_name
            else:
                sequence_id = f"{clean_template_name}_{cut_site.grna_name}"

        # Build the Primer3 input
        lines = []

        # Sequence tags
        lines.append(f"SEQUENCE_ID={sequence_id}")
        lines.append(f"SEQUENCE_TEMPLATE={template.sequence}")
        lines.append(f"SEQUENCE_TARGET={target_start},{target_length}")

        # Global primer parameters
        lines.append(f"PRIMER_TASK=generic")
        lines.append(f"PRIMER_PICK_LEFT_PRIMER=1")
        lines.append(f"PRIMER_PICK_RIGHT_PRIMER=1")
        lines.append(f"PRIMER_PICK_INTERNAL_OLIGO=0")

        # Size parameters
        lines.append(f"PRIMER_PRODUCT_SIZE_RANGE={self.config.product_size_range}")
        lines.append(f"PRIMER_MIN_SIZE={self.config.primer_min_size}")
        lines.append(f"PRIMER_OPT_SIZE={self.config.primer_opt_size}")
        lines.append(f"PRIMER_MAX_SIZE={self.config.primer_max_size}")

        # Temperature parameters
        lines.append(f"PRIMER_MIN_TM={self.config.primer_min_tm}")
        lines.append(f"PRIMER_OPT_TM={self.config.primer_opt_tm}")
        lines.append(f"PRIMER_MAX_TM={self.config.primer_max_tm}")
        lines.append(f"PRIMER_PAIR_MAX_DIFF_TM={self.config.primer_pair_max_diff_tm}")

        # GC content
        lines.append(f"PRIMER_MIN_GC={self.config.primer_min_gc}")
        lines.append(f"PRIMER_OPT_GC_PERCENT={self.config.primer_opt_gc}")
        lines.append(f"PRIMER_MAX_GC={self.config.primer_max_gc}")

        # Number of results
        lines.append(f"PRIMER_NUM_RETURN={self.config.num_return}")

        # Self-complementarity
        lines.append(f"PRIMER_MAX_SELF_ANY={self.config.primer_max_self_any}")
        lines.append(f"PRIMER_MAX_SELF_END={self.config.primer_max_self_end}")
        lines.append(f"PRIMER_PAIR_MAX_COMPL_ANY={self.config.primer_pair_max_compl_any}")
        lines.append(f"PRIMER_PAIR_MAX_COMPL_END={self.config.primer_pair_max_compl_end}")

        # Thermodynamic parameters
        lines.append(f"PRIMER_SALT_MONOVALENT={self.config.primer_salt_monovalent}")
        lines.append(f"PRIMER_SALT_DIVALENT={self.config.primer_salt_divalent}")
        lines.append(f"PRIMER_DNTP_CONC={self.config.primer_dntp_conc}")
        lines.append(f"PRIMER_DNA_CONC={self.config.primer_dna_conc}")

        # Other parameters
        lines.append(f"PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT={self.config.primer_thermodynamic_oligo_alignment}")
        lines.append(f"PRIMER_MAX_NS_ACCEPTED={self.config.primer_max_ns_accepted}")
        lines.append(f"PRIMER_MAX_POLY_X={self.config.primer_max_poly_x}")
        lines.append(f"PRIMER_GC_CLAMP={self.config.primer_gc_clamp}")
        lines.append(f"PRIMER_MAX_END_GC={self.config.primer_max_end_gc}")
        lines.append(f"PRIMER_EXPLAIN_FLAG={self.config.primer_explain_flag}")

        # Add any extra custom parameters
        for key, value in self.config.extra_params.items():
            lines.append(f"{key}={value}")

        # Record terminator
        lines.append("=")

        return "\n".join(lines)


class SequenceExtractor:
    """Extract sequences from a genome FASTA using samtools faidx"""

    def __init__(self, genome_fa: str):
        self.genome_fa = genome_fa
        self.index_path = f"{genome_fa}.fai"
        self.index = {}

        # Load index if it exists
        if os.path.exists(self.index_path):
            self._load_index()
        else:
            print(f"Warning: Index file {self.index_path} not found.")
            print(f"Attempting to create index with samtools...")
            try:
                subprocess.run(["samtools", "faidx", genome_fa], check=True)
                self._load_index()
            except (subprocess.CalledProcessError, FileNotFoundError):
                print(f"Error: Failed to index genome. Ensure samtools is installed or provide a .fai file.")

    def _load_index(self):
        """Parse the FASTA index (.fai) file"""
        try:
            with open(self.index_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        # name, length, offset, line_bases, line_width
                        self.index[parts[0]] = {
                            'length': int(parts[1]),
                            'offset': int(parts[2]),
                            'line_bases': int(parts[3]),
                            'line_width': int(parts[4])
                        }
        except Exception as e:
            print(f"Error loading index: {e}")

    def extract_region(self, chromosome: str, start: int, end: int) -> Optional[str]:
        """Extract a specific genomic region"""
        region = f"{chromosome}:{start}-{end}"
        
        # Try a few variations of the chromosome name
        chrom_variants = [chromosome]
        if chromosome.lower().startswith("chr"):
            chrom_variants.append(chromosome[3:])
        else:
            chrom_variants.append(f"chr{chromosome}")
            chrom_variants.append(f"CHR{chromosome}")

        # Try samtools first if available
        try:
            # Check if samtools is even in the path first to avoid subprocess error display
            if subprocess.run(["which", "samtools"], capture_output=True).returncode == 0:
                for chrom in chrom_variants:
                    cmd = ["samtools", "faidx", self.genome_fa, f"{chrom}:{start}-{end}"]
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    if result.returncode == 0:
                        lines = result.stdout.strip().split('\n')
                        if len(lines) > 1:
                            return "".join(lines[1:])
        except Exception:
            pass

        # Fallback to pure Python extraction using index
        for chrom in chrom_variants:
            if chrom in self.index:
                idx = self.index[chrom]
                
                # Check bounds
                start = max(1, start)
                end = min(idx['length'], end)
                if start > end:
                    return None
                
                def get_file_offset(pos):
                    # 1-based position to 0-based file offset
                    return idx['offset'] + (pos - 1) + ((pos - 1) // idx['line_bases']) * (idx['line_width'] - idx['line_bases'])
                
                start_offset = get_file_offset(start)
                end_offset = get_file_offset(end)
                
                try:
                    with open(self.genome_fa, 'rb') as f:
                        f.seek(start_offset)
                        data = f.read(end_offset - start_offset + 1)
                        # Remove newlines and decode
                        return data.decode('ascii', errors='ignore').replace('\n', '').replace('\r', '')
                except Exception as e:
                    print(f"Python extraction error for {region} on {chrom}: {e}")
        
        print(f"Error: Chromosome {chromosome} (or variants) not found in index.")
        return None


def _scan_chunk_worker(chunk_info, sequences, return_locations=False):
    """Worker function to count or find sequence occurrences in a file chunk"""
    file_path, start_offset, size, overlap, chunk_chrom, chunk_chrom_start = chunk_info
    
    # chunk_chrom and chunk_chrom_start are used when scanning by chromosome
    # if they are None, we are doing a generic byte-scan (current behavior)
    
    results = {
        'counts': {seq: 0 for seq in sequences},
        'locations': {seq: [] for seq in sequences} if return_locations else None
    }

    try:
        with open(file_path, 'rb') as f:
            f.seek(start_offset)
            # Read size + overlap to handle span-overs
            data = f.read(size + overlap).upper()
            text = data.decode('ascii', errors='ignore')

            if chunk_chrom:
                # We are scanning a specific chromosome part
                # text should be pure sequence data (no headers)
                # But wait, the file might have newlines
                cleaned_text = re.sub(r'\s+', '', text)
                
                # We need to map positions in cleaned_text back to genomic positions
                # This is tricky because of newlines. 
                # Simpler: find in raw text but skip newlines in match length
                for seq in sequences:
                    # Search for sequence, allowing for internal whitespace/newlines
                    # Regex for "A T G" with any whitespace between
                    pattern = "".join([f"{base}[\\s]*" for base in seq])
                    for match in re.finditer(pattern, text):
                        # Match position in the raw chunk text
                        match_start = match.start()
                        
                        # Calculate genomic position if we have the chromosome info
                        # We need to know how many non-sequence characters are before match_start
                        pre_match = text[:match_start]
                        actual_bases_before = len(re.sub(r'\s+', '', pre_match))
                        
                        genomic_pos = chunk_chrom_start + actual_bases_before
                        
                        results['counts'][seq] += 1
                        if return_locations:
                            results['locations'][seq].append(f"{chunk_chrom}:{genomic_pos}")
            else:
                # Legacy generic scan (used for quick uniqueness check)
                # Split by '>' to isolate headers
                parts = text.split('>')
                cleaned_text_parts = []

                for i, part in enumerate(parts):
                    if i == 0 and start_offset > 0:
                        cleaned_text_parts.append(re.sub(r'\s+', '', part))
                    else:
                        if '\n' in part:
                            seq_data = part.split('\n', 1)[1]
                            cleaned_text_parts.append(re.sub(r'\s+', '', seq_data))

                cleaned_text = "".join(cleaned_text_parts)
                for seq in sequences:
                    results['counts'][seq] = cleaned_text.count(seq)

    except Exception as e:
        print(f"Worker error at offset {start_offset}: {e}")

    return results


class ParallelGenomeScanner:
    """Multi-processed scanner to count sequence occurrences or find locations in a FASTA file"""

    def __init__(self, genome_fa: str, num_processes: Optional[int] = None):
        self.genome_fa = genome_fa
        self.num_processes = num_processes or multiprocessing.cpu_count()
        self.chunk_size = 64 * 1024 * 1024  # 64MB chunks
        self.overlap = 1000  # 1kb overlap to handle primers crossing chunks
        
        # Load index for location mapping
        self.index = {}
        index_path = f"{genome_fa}.fai"
        if os.path.exists(index_path):
            with open(index_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        self.index[parts[0]] = {
                            'length': int(parts[1]),
                            'offset': int(parts[2]),
                            'line_bases': int(parts[3]),
                            'line_width': int(parts[4])
                        }

    def count_occurrences(self, sequences: List[str]) -> Dict[str, int]:
        """Count occurrences of multiple sequences in parallel (legacy behavior)"""
        if not sequences:
            return {}

        file_size = os.path.getsize(self.genome_fa)
        chunks = []
        for offset in range(0, file_size, self.chunk_size):
            # (file_path, start_offset, size, overlap, chunk_chrom, chunk_chrom_start)
            chunks.append((self.genome_fa, offset, self.chunk_size, self.overlap, None, None))

        print(f"Scanning genome with {self.num_processes} processes ({len(chunks)} chunks)...")
        start_time = time.time()
        sequences = [s.upper() for s in sequences]

        with multiprocessing.Pool(processes=self.num_processes) as pool:
            worker_func = partial(_scan_chunk_worker, sequences=sequences, return_locations=False)
            results = pool.map(worker_func, chunks)

        total_counts = {seq: 0 for seq in sequences}
        for res in results:
            for seq, count in res['counts'].items():
                total_counts[seq] += count

        duration = time.time() - start_time
        print(f"Scan complete in {duration:.1f}s")
        return total_counts

    def find_locations(self, sequences: List[str]) -> Dict[str, List[str]]:
        """Find all genomic locations of sequences in parallel"""
        if not sequences or not self.index:
            return {seq: [] for seq in sequences}

        # Create jobs based on chromosome index
        chunks = []
        for chrom, info in self.index.items():
            chrom_len = info['length']
            chrom_offset = info['offset']
            line_width = info['line_width']
            line_bases = info['line_bases']
            
            # Approximate byte size of chromosome in file
            chrom_byte_size = (chrom_len // line_bases) * line_width + (chrom_len % line_bases)
            
            # Split large chromosomes into manageable chunks
            for start_base in range(1, chrom_len + 1, 1000000): # 1MB chunks of bases
                end_base = min(chrom_len, start_base + 1000000 - 1)
                
                # Calculate byte offsets
                def get_offset(pos):
                    return chrom_offset + ((pos - 1) // line_bases) * line_width + ((pos - 1) % line_bases)
                
                job_start_offset = get_offset(start_base)
                job_end_offset = get_offset(end_base)
                job_size = job_end_offset - job_start_offset + 1
                
                # Add overlap (about 100 bases is enough for primers)
                job_overlap = 200
                
                chunks.append((self.genome_fa, job_start_offset, job_size, job_overlap, chrom, start_base))

        print(f"Finding locations in {len(self.index)} chromosomes ({len(chunks)} chunks)...")
        start_time = time.time()
        sequences = [s.upper() for s in sequences]

        with multiprocessing.Pool(processes=self.num_processes) as pool:
            worker_func = partial(_scan_chunk_worker, sequences=sequences, return_locations=True)
            results = pool.map(worker_func, chunks)

        all_locations = {seq: [] for seq in sequences}
        for res in results:
            for seq, locs in res['locations'].items():
                all_locations[seq].extend(locs)

        # Deduplicate (overlaps might cause double counting)
        for seq in all_locations:
            all_locations[seq] = sorted(list(set(all_locations[seq])))

        duration = time.time() - start_time
        print(f"Location search complete in {duration:.1f}s")
        return all_locations


class Primer3OutputParser:
    """Parse Primer3 output into structured data"""

    @dataclass
    class PrimerPair:
        """Represents a primer pair result"""
        pair_index: int
        left_sequence: str
        left_position: Tuple[int, int]  # (start, length)
        left_tm: float
        left_gc: float
        right_sequence: str
        right_position: Tuple[int, int]
        right_tm: float
        right_gc: float
        product_size: int
        pair_penalty: float
        amplicon_sequence: str = ""

    @dataclass
    class Primer3Result:
        """Complete result from a Primer3 run"""
        sequence_id: str
        num_returned: int
        primer_pairs: List['Primer3OutputParser.PrimerPair']
        target_gene: str = ""
        location: str = ""
        extra_data: Dict[str, str] = field(default_factory=dict)
        left_explain: str = ""
        right_explain: str = ""
        pair_explain: str = ""
        errors: List[str] = field(default_factory=list)
        warnings: List[str] = field(default_factory=list)

    def parse_output(self, output_text: str) -> List['Primer3OutputParser.Primer3Result']:
        """Parse Primer3 Boulder-IO output into structured results"""
        results = []

        # Split into records (separated by '=')
        records = output_text.strip().split('\n=\n')

        for record in records:
            if not record.strip():
                continue

            result = self._parse_record(record)
            if result:
                results.append(result)

        return results

    def _parse_record(self, record: str) -> Optional['Primer3OutputParser.Primer3Result']:
        """Parse a single Primer3 output record"""
        lines = record.strip().split('\n')
        data = {}

        for line in lines:
            if '=' in line:
                key, value = line.split('=', 1)
                data[key] = value

        if not data:
            return None

        # Extract template sequence for amplicon generation
        template_seq = data.get('SEQUENCE_TEMPLATE', '')

        # Extract basic info
        sequence_id = data.get('SEQUENCE_ID', 'unknown')
        num_returned = int(data.get('PRIMER_PAIR_NUM_RETURNED', 0))

        # Extract errors and warnings
        errors = []
        warnings = []
        if 'PRIMER_ERROR' in data and data['PRIMER_ERROR']:
            errors.append(data['PRIMER_ERROR'])
        if 'PRIMER_WARNING' in data and data['PRIMER_WARNING']:
            warnings.append(data['PRIMER_WARNING'])

        # Extract explain strings
        left_explain = data.get('PRIMER_LEFT_EXPLAIN', '')
        right_explain = data.get('PRIMER_RIGHT_EXPLAIN', '')
        pair_explain = data.get('PRIMER_PAIR_EXPLAIN', '')

        # Extract primer pairs
        primer_pairs = []
        for i in range(num_returned):
            try:
                # Left primer
                left_seq = data.get(f'PRIMER_LEFT_{i}_SEQUENCE', '')
                left_pos_str = data.get(f'PRIMER_LEFT_{i}', '0,0')
                left_pos = tuple(map(int, left_pos_str.split(',')))
                left_tm = float(data.get(f'PRIMER_LEFT_{i}_TM', 0))
                left_gc = float(data.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 0))

                # Right primer
                right_seq = data.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')
                right_pos_str = data.get(f'PRIMER_RIGHT_{i}', '0,0')
                right_pos = tuple(map(int, right_pos_str.split(',')))
                right_tm = float(data.get(f'PRIMER_RIGHT_{i}_TM', 0))
                right_gc = float(data.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 0))

                # Pair info
                product_size = int(data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0))
                pair_penalty = float(data.get(f'PRIMER_PAIR_{i}_PENALTY', 0))

                # Calculate amplicon sequence if template available
                amplicon_seq = ""
                if template_seq:
                    # In Primer3, PRIMER_LEFT is (5_end, length)
                    # PRIMER_RIGHT is (3_end, length)
                    start = left_pos[0]
                    end = right_pos[0]
                    amplicon_seq = template_seq[start : end + 1]

                pair = self.PrimerPair(
                    pair_index=i,
                    left_sequence=left_seq,
                    left_position=left_pos,
                    left_tm=left_tm,
                    left_gc=left_gc,
                    right_sequence=right_seq,
                    right_position=right_pos,
                    right_tm=right_tm,
                    right_gc=right_gc,
                    product_size=product_size,
                    pair_penalty=pair_penalty,
                    amplicon_sequence=amplicon_seq
                )
                primer_pairs.append(pair)
            except (ValueError, KeyError) as e:
                continue

        return self.Primer3Result(
            sequence_id=sequence_id,
            num_returned=num_returned,
            primer_pairs=primer_pairs,
            left_explain=left_explain,
            right_explain=right_explain,
            pair_explain=pair_explain,
            errors=errors,
            warnings=warnings
        )


class CRISPRPrimerDesigner:
    """Main class to orchestrate CRISPR primer design"""

    def __init__(self, config: Optional[Primer3Config] = None,
                 primer3_path: str = "primer3_core"):
        self.config = config or Primer3Config()

        # Try to find primer3_core if not provided as absolute path
        if not os.path.isabs(primer3_path):
            # Calculate path relative to the script
            script_dir = Path(__file__).parent.absolute()
            local_p3 = script_dir / "primer3" / "src" / "primer3_core"
            
            common_paths = [
                str(local_p3),
                "/Users/munchr/bin/CRISPR_primer_designer/primer3/src/primer3_core",
                "/Users/munchr/Documents/GitHub/CRISPR_primer_designer/primer3/src/primer3_core",
                "/usr/local/bin/primer3_core",
                "/usr/bin/primer3_core",
                "/opt/homebrew/bin/primer3_core",
                "primer3_core"
            ]
            for p in common_paths:
                if os.path.exists(p) or subprocess.run(["which", p], capture_output=True).returncode == 0:
                    primer3_path = p
                    break

        self.primer3_path = primer3_path
        self.input_generator = Primer3InputGenerator(self.config)
        self.output_parser = Primer3OutputParser()
        self._id_metadata_map = {}  # Store metadata for each generated sequence_id

    def load_cutsites_from_tsv(self, tsv_path: str) -> List[CutSite]:
        """Load cut sites from a TSV file"""
        cutsites = []

        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')

            # Normalize column names (handle different cases/formats)
            for row in reader:
                # Try to find the right columns
                target_gene = (row.get('Target gene') or row.get('target_gene') or
                              row.get('Gene') or row.get('gene') or '')
                grna_name = (row.get('gRNA name') or row.get('grna_name') or
                            row.get('Name') or row.get('name') or '')
                sequence = (row.get('Sequence') or row.get('sequence') or
                           row.get('gRNA') or row.get('grna') or '')
                location = (row.get('location') or row.get('Location') or
                           row.get('Position') or row.get('position') or '')

                # At minimum, we need a sequence
                if sequence:
                    # Store everything else in extra_data
                    known_cols = {'Target gene', 'target_gene', 'Gene', 'gene',
                                 'gRNA name', 'grna_name', 'Name', 'name',
                                 'Sequence', 'sequence', 'gRNA', 'grna',
                                 'location', 'Location', 'Position', 'position'}
                    extra_data = {k: v for k, v in row.items() if k not in known_cols}

                    cutsites.append(CutSite(
                        target_gene=target_gene,
                        grna_name=grna_name or f"gRNA_{len(cutsites)+1}",
                        sequence=sequence,
                        location=location,
                        extra_data=extra_data
                    ))

        return cutsites

    def load_templates_from_fasta(self, fasta_path: str) -> List[TemplateSequence]:
        """Load template sequences from a FASTA file"""
        templates = []
        current_name = None
        current_seq = []
        current_chrom = None
        current_start = None
        current_end = None

        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_name and current_seq:
                        templates.append(TemplateSequence(
                            name=current_name,
                            sequence=''.join(current_seq),
                            chromosome=current_chrom,
                            start_position=current_start,
                            end_position=current_end
                        ))

                    # Parse new header
                    header = line[1:].strip()
                    current_name = header.split()[0]
                    current_seq = []

                    # Try to extract coordinates from header
                    # Format: >name chr:start-end or >chr:start-end
                    coord_match = re.search(r'(chr[\w\.]+):(\d+)-(\d+)', header, re.IGNORECASE)
                    if coord_match:
                        current_chrom = coord_match.group(1)
                        current_start = int(coord_match.group(2))
                        current_end = int(coord_match.group(3))
                    elif ':' in header and '-' in header:
                         # Try more generic coordinate match
                         gen_match = re.search(r'([\w\.]+):(\d+)-(\d+)', header)
                         if gen_match:
                             current_chrom = gen_match.group(1)
                             current_start = int(gen_match.group(2))
                             current_end = int(gen_match.group(3))
                    else:
                        current_chrom = None
                        current_start = None
                        current_end = None
                elif line:
                    current_seq.append(line.upper())

            # Don't forget the last sequence
            if current_name and current_seq:
                templates.append(TemplateSequence(
                    name=current_name,
                    sequence=''.join(current_seq),
                    chromosome=current_chrom,
                    start_position=current_start,
                    end_position=current_end
                ))

        return templates

    def create_template_from_string(self, name: str, sequence: str,
                                     chromosome: str = None,
                                     start_position: int = None) -> TemplateSequence:
        """Create a template from a sequence string"""
        return TemplateSequence(
            name=name,
            sequence=sequence.upper().replace(" ", "").replace("\n", ""),
            chromosome=chromosome,
            start_position=start_position,
            end_position=start_position + len(sequence) if start_position else None
        )

    def match_cutsites_to_templates(self, cutsites: List[CutSite],
                                     templates: List[TemplateSequence]) -> Dict[str, List[CutSite]]:
        """
        Match cut sites to templates based on sequence matching or location.
        Returns a dict mapping template names to lists of cut sites found in them.
        This version scans ALL templates to find all matches, rather than stopping at the first.
        """
        matches = {t.name: [] for t in templates}
        # Reset metadata map for this run
        self._id_metadata_map = {}

        for cutsite in cutsites:
            templates_found = set()
            
            # 1. Try direct name matching
            for template in templates:
                if cutsite.grna_name == template.name:
                    matches[template.name].append(cutsite)
                    templates_found.add(template.name)

            # 2. Try sequence/location matching
            for template in templates:
                if template.name in templates_found:
                    continue

                found = self.input_generator.find_sequence_in_template(
                    cutsite.sequence, template.sequence
                )

                if found:
                    matches[template.name].append(cutsite)
                    templates_found.add(template.name)
                    continue

                # If not found by sequence, try by location matching
                chrom, start, end = self.input_generator.parse_location(cutsite.location)
                if template.chromosome and chrom:
                    if (template.chromosome.lower() == chrom.lower() and
                        template.start_position is not None and
                        template.start_position <= start <= template.end_position):
                        matches[template.name].append(cutsite)
                        templates_found.add(template.name)

        return matches

    def generate_templates_from_genome(self, cutsites: List[CutSite],
                                       genome_fa: str,
                                       flank: int = 250) -> List[TemplateSequence]:
        """
        Generate TemplateSequence objects by extracting regions from a genome.
        """
        extractor = SequenceExtractor(genome_fa)
        templates = []

        for cutsite in cutsites:
            chrom, start, end = self.input_generator.parse_location(cutsite.location)
            if chrom and start:
                # Calculate flanking region
                extract_start = max(1, start - flank)
                extract_end = (end or start) + flank

                seq = extractor.extract_region(chrom, extract_start, extract_end)
                if seq:
                    # Create a unique template for this cutsite
                    template_name = cutsite.grna_name
                    templates.append(TemplateSequence(
                        name=template_name,
                        sequence=seq,
                        chromosome=chrom,
                        start_position=extract_start,
                        end_position=extract_end
                    ))

        return templates

    def filter_unique_primers(self, results: List[Primer3OutputParser.Primer3Result],
                              genome_fa: str) -> List[Primer3OutputParser.Primer3Result]:
        """
        Filter out primer pairs where either primer is non-unique in the genome.
        """
        if not results:
            return []

        # 1. Collect all sequences to check (including reverse complements)
        # We check both to ensure uniqueness on both strands
        seq_to_check = set()
        primer_to_rc = {}

        for res in results:
            for pair in res.primer_pairs:
                for seq in [pair.left_sequence, pair.right_sequence]:
                    seq_up = seq.upper()
                    seq_to_check.add(seq_up)
                    rc = self.input_generator._reverse_complement(seq_up)
                    seq_to_check.add(rc)
                    primer_to_rc[seq_up] = rc

        # 2. Run parallel scan
        scanner = ParallelGenomeScanner(genome_fa)
        counts = scanner.count_occurrences(list(seq_to_check))

        # 3. Filter results
        filtered_results = []
        total_removed = 0

        for res in results:
            unique_pairs = []
            for pair in res.primer_pairs:
                left = pair.left_sequence.upper()
                right = pair.right_sequence.upper()

                left_rc = primer_to_rc[left]
                right_rc = primer_to_rc[right]

                # Check for uniqueness: total occurrences of (seq + its rc) should be 1
                # If total is 1, it means it only exists once in the genome (either on + or - strand)
                # If total is > 1, it has another binding site.
                left_total = counts.get(left, 0) + counts.get(left_rc, 0)
                right_total = counts.get(right, 0) + counts.get(right_rc, 0)

                if left_total == 1 and right_total == 1:
                    unique_pairs.append(pair)
                else:
                    total_removed += 1

            if unique_pairs:
                res.primer_pairs = unique_pairs
                res.num_returned = len(unique_pairs)
                filtered_results.append(res)
            else:
                # If no pairs are unique, we still keep the result object but with 0 pairs
                res.primer_pairs = []
                res.num_returned = 0
                res.warnings.append("All designed primers were non-unique in the genome and were filtered out.")
                filtered_results.append(res)

        print(f"Uniqueness filtering: Removed {total_removed} non-unique primer pairs")
        return filtered_results

    def report_primer_locations(self, results: List[Primer3OutputParser.Primer3Result],
                                genome_fa: str, output_file: str):
        """
        Find and report all genomic locations for all primers in the results.
        Saves to a TSV file.
        """
        if not results:
            return

        # 1. Collect unique primer sequences
        all_primers = set()
        primer_to_rc = {}
        for res in results:
            for pair in res.primer_pairs:
                for seq in [pair.left_sequence, pair.right_sequence]:
                    seq_up = seq.upper()
                    all_primers.add(seq_up)
                    rc = self.input_generator._reverse_complement(seq_up)
                    all_primers.add(rc)
                    primer_to_rc[seq_up] = rc

        # 2. Find all locations
        scanner = ParallelGenomeScanner(genome_fa)
        locs_found = scanner.find_locations(list(all_primers))

        # 3. Write report
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['Guide_Name', 'Primer_Type', 'Sequence', 'Location', 'Strand'])

            for res in results:
                for pair in res.primer_pairs:
                    # Report for left primer
                    l_seq = pair.left_sequence.upper()
                    l_rc = primer_to_rc[l_seq]
                    
                    # Forward matches
                    for loc in locs_found.get(l_seq, []):
                        writer.writerow([res.sequence_id, 'Left', l_seq, loc, '+'])
                    # Reverse complement matches
                    for loc in locs_found.get(l_rc, []):
                        writer.writerow([res.sequence_id, 'Left', l_seq, loc, '-'])

                    # Report for right primer
                    r_seq = pair.right_sequence.upper()
                    r_rc = primer_to_rc[r_seq]
                    
                    # Forward matches
                    for loc in locs_found.get(r_seq, []):
                        writer.writerow([res.sequence_id, 'Right', r_seq, loc, '+'])
                    # Reverse complement matches
                    for loc in locs_found.get(r_rc, []):
                        writer.writerow([res.sequence_id, 'Right', r_seq, loc, '-'])

        print(f"Saved all primer locations to: {output_file}")

    def _has_poly_x(self, sequence: str, max_poly_x: int = 3) -> bool:
        """Check if sequence has more than max_poly_x identical bases in a row"""
        sequence = sequence.upper()
        for base in "ACGT":
            if base * (max_poly_x + 1) in sequence:
                return True
        return False

    def filter_poly_x_primers(self, results: List[Primer3OutputParser.Primer3Result],
                               max_poly_x: int = 3) -> List[Primer3OutputParser.Primer3Result]:
        """
        Filter out primer pairs where either primer contains a mononucleotide repeat
        longer than max_poly_x (e.g., AAAA if max_poly_x=3).
        """
        if not results:
            return []

        filtered_results = []
        total_removed = 0

        for res in results:
            valid_pairs = []
            for pair in res.primer_pairs:
                left_bad = self._has_poly_x(pair.left_sequence, max_poly_x)
                right_bad = self._has_poly_x(pair.right_sequence, max_poly_x)

                if not left_bad and not right_bad:
                    valid_pairs.append(pair)
                else:
                    total_removed += 1

            # Update the result object
            res.primer_pairs = valid_pairs
            res.num_returned = len(valid_pairs)
            if not valid_pairs and not res.errors:
                res.warnings.append(f"All designed primers contained mononucleotide repeats > {max_poly_x} and were filtered out.")
            
            filtered_results.append(res)

        if total_removed > 0:
            print(f"Poly-X filtering: Removed {total_removed} primer pairs with mononucleotide repeats > {max_poly_x}")
        
        return filtered_results

    def generate_all_inputs(self, cutsites: List[CutSite],
                            templates: List[TemplateSequence],
                            output_file: str = None) -> str:
        """
        Generate Primer3 input for all cut site/template combinations.
        """
        # Match cutsites to templates
        matches = self.match_cutsites_to_templates(cutsites, templates)

        all_inputs = []

        for template in templates:
            matched_cutsites = matches[template.name]

            for cutsite in matched_cutsites:
                input_str = self.input_generator.generate_primer3_input(
                    cutsite, template
                )
                if input_str:
                    all_inputs.append(input_str)

                    # Capture the generated ID and store metadata
                    # ID line is the first line: SEQUENCE_ID=name
                    id_line = input_str.split('\n')[0]
                    seq_id = id_line.split('=', 1)[1]
                    self._id_metadata_map[seq_id] = {
                        'target_gene': cutsite.target_gene,
                        'location': cutsite.location,
                        'extra_data': cutsite.extra_data
                    }

        combined_input = "\n".join(all_inputs)

        if output_file:
            # Ensure the output directory exists if it's a path
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, 'w') as f:
                f.write(combined_input)

        return combined_input

    def run_primer3(self, input_text: str,
                    thermodynamic_params_path: str = None) -> str:
        """Run Primer3 with the given input"""
        cmd = [self.primer3_path]

        if thermodynamic_params_path:
            cmd.extend([f'--p3_settings_file={thermodynamic_params_path}'])

        try:
            result = subprocess.run(
                cmd,
                input=input_text,
                capture_output=True,
                text=True,
                check=True
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            print(f"Primer3 error: {e.stderr}")
            raise
        except FileNotFoundError:
            print(f"Primer3 executable not found at: {self.primer3_path}")
            print("Please install Primer3 or specify the correct path.")
            raise

    def design_primers(self, cutsites: List[CutSite],
                       templates: List[TemplateSequence],
                       output_prefix: str = "primers",
                       output_dir: Optional[str] = None) -> List[Primer3OutputParser.Primer3Result]:
        """
        Main method to design primers for all cut sites.

        Returns parsed Primer3 results and saves output files.
        """
        # Set up output paths
        prefix_path = Path(output_prefix)
        if output_dir:
            out_dir = Path(output_dir)
            # We already checked existence in main, but let's be safe
            out_dir.mkdir(parents=True, exist_ok=True)
            base_output = out_dir / prefix_path.name
        else:
            base_output = prefix_path

        # Generate input
        input_file = f"{base_output}_input.txt"
        primer3_input = self.generate_all_inputs(cutsites, templates, input_file)

        print(f"Generated Primer3 input: {input_file}")
        print(f"Total records: {primer3_input.count('=')}")

        # Run Primer3
        try:
            output_text = self.run_primer3(primer3_input)

            # Save raw output
            output_file = f"{base_output}_output.txt"
            with open(output_file, 'w') as f:
                f.write(output_text)
            print(f"Saved Primer3 output: {output_file}")

            # Parse results
            results = self.output_parser.parse_output(output_text)

            # Map metadata back to results using the internal map
            for res in results:
                meta = self._id_metadata_map.get(res.sequence_id)
                if meta:
                    res.target_gene = meta['target_gene']
                    res.location = meta['location']
                    res.extra_data = meta['extra_data']

            # Filter for mononucleotide repeats (Poly-X)
            # Default is max 3 identical bases in a row (remove 4 or more)
            results = self.filter_poly_x_primers(results, max_poly_x=self.config.primer_max_poly_x)

            # Save parsed results as TSV
            self._save_results_tsv(results, f"{base_output}_primers.tsv")
            
            # Save simplified list for SnapGene
            self._save_snapgene_tsv(results, f"{base_output}_snapgene.tsv")
            
            # Save bulk order for oligo ordering
            self._save_bulk_order_tsv(results, f"{base_output}_bulk_order.tsv")

            return results

        except Exception as e:
            print(f"Error running Primer3: {e}")
            print("Primer3 input has been saved. You can run it manually.")
            return []

    def _save_results_tsv(self, results: List[Primer3OutputParser.Primer3Result],
                          output_file: str):
        """Save primer results to a TSV file"""
        # Ensure directory exists
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            # Determine all extra columns present in results
            extra_cols = []
            for res in results:
                for k in res.extra_data.keys():
                    if k not in extra_cols:
                        extra_cols.append(k)

            # Header
            header = ['Guide_Name', 'Location', 'Sequence_ID', 'Pair_Index'] + extra_cols + [
                'Left_Primer', 'Left_Position', 'Left_Length', 'Left_Tm', 'Left_GC',
                'Right_Primer', 'Right_Position', 'Right_Length', 'Right_Tm', 'Right_GC',
                'Product_Size', 'Pair_Penalty', 'Amplicon_Sequence'
            ]
            writer.writerow(header)

            for result in results:
                for pair in result.primer_pairs:
                    row = [
                        result.target_gene,
                        result.location,
                        result.sequence_id,
                        pair.pair_index
                    ]
                    # Add extra data columns
                    for col in extra_cols:
                        row.append(result.extra_data.get(col, ''))

                    # Add primer data
                    row.extend([
                        pair.left_sequence,
                        pair.left_position[0],
                        pair.left_position[1],
                        f"{pair.left_tm:.1f}",
                        f"{pair.left_gc:.1f}",
                        pair.right_sequence,
                        pair.right_position[0],
                        pair.right_position[1],
                        f"{pair.right_tm:.1f}",
                        f"{pair.right_gc:.1f}",
                        pair.product_size,
                        f"{pair.pair_penalty:.2f}",
                        pair.amplicon_sequence
                    ])
                    writer.writerow(row)

        print(f"Saved primer results: {output_file}")

    def _save_snapgene_tsv(self, results: List[Primer3OutputParser.Primer3Result],
                            output_file: str):
        """
        Save simplified primer list for SnapGene import (Name[tab]Sequence).
        Duplicates by sequence are removed, keeping the first name encountered.
        """
        # Ensure directory exists
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        seen_sequences = set()
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            for result in results:
                for pair in result.primer_pairs:
                    # Left primer
                    if pair.left_sequence not in seen_sequences:
                        left_name = f"{result.sequence_id}_{pair.pair_index}_L"
                        writer.writerow([left_name, pair.left_sequence])
                        seen_sequences.add(pair.left_sequence)
                    
                    # Right primer
                    if pair.right_sequence not in seen_sequences:
                        right_name = f"{result.sequence_id}_{pair.pair_index}_R"
                        writer.writerow([right_name, pair.right_sequence])
                        seen_sequences.add(pair.right_sequence)
        
        print(f"Saved SnapGene primer list: {output_file} ({len(seen_sequences)} unique primers)")

    def _save_bulk_order_tsv(self, results: List[Primer3OutputParser.Primer3Result],
                             output_file: str):
        """Save bulk order TSV for oligo ordering with TruSeq adapters"""
        # TruSeq Adapters
        LEFT_ADAPTER = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        RIGHT_ADAPTER = "GACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
        
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        seen_sequences = set()
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            for result in results:
                for pair in result.primer_pairs:
                    # Left primer
                    full_left = LEFT_ADAPTER + pair.left_sequence
                    if full_left not in seen_sequences:
                        left_name = f"{result.sequence_id}_{pair.pair_index}_L"
                        writer.writerow([left_name, full_left, "25nm", "STD"])
                        seen_sequences.add(full_left)
                    
                    # Right primer
                    full_right = RIGHT_ADAPTER + pair.right_sequence
                    if full_right not in seen_sequences:
                        right_name = f"{result.sequence_id}_{pair.pair_index}_R"
                        writer.writerow([right_name, full_right, "25nm", "STD"])
                        seen_sequences.add(full_right)
        
        print(f"Saved Bulk Order list: {output_file} ({len(seen_sequences)} oligos)")


def create_example_files():
    """Create example input files for testing"""

    # Example TSV file
    tsv_content = """Target gene\tgRNA name\tSequence\tlocation
BRCA1\tBRCA1_gRNA1\tGCTGACTTACCAGATGGGAC\tchr17:43044295-43044314
BRCA1\tBRCA1_gRNA2\tACTGGATCCAGATGACGTAC\tchr17:43045100-43045119
TP53\tTP53_gRNA1\tGTCCAGATGACCAGGTGCAT\tchr17:7577000-7577019
"""

    with open('example_cutsites.tsv', 'w') as f:
        f.write(tsv_content)

    # Example FASTA file with template sequences
    fasta_content = """>BRCA1_region chr17:43044000-43051000
ATGGATTTCTGCTGCTCGCGCTACTCTCTCTCTGTCTGGCCTGGAGGCTATCCAGCGTGAGTCTCTCCTACCCTCCCGCT
GGGCCTGTAGCGGGGCTGTGGTCGCAGGGCCAACATAGGAGAAGAATCTCCCGGTTTGTCTGTCCACGCGCTCTGGCTGC
TGAACGCCCTCTGCTCCCAGCTGGCAGCCAGGGAATGGCAGAAGGCAAGAACCCTCCCCTCCTCACCCTCATCACCGAGC
CCGGCAGAAGCTCCCAGGATGCTGACTTTGTCACCAGTCCGGGAAGCTCTTGGCAGAAGACGCGGCAAGCAGCAGAGCGA
GCTGACTTACCAGATGGGACAGGCTGCCCGCCGCCTCAGGAAGTAAGGAACTGCACGCTGGCGTGGTGGCTCACGCCTGT
AATCCTAGCACTTTGGGAGGCCGAGGCGGGTGGATCATGAGGTCAGGAGATCGAGACCATCCTGGCTAACATGGTGAAAC
>TP53_region chr17:7576500-7578000
GCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGG
GCCCCAGGCCTCTGATTCCTCACTGATTGCTCTTAGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTG
CGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGA
CTGTCCAGATGACCAGGTGCATGGTGCCAGCCCTGGCTCCCCAGAATGCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCA
GCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCA
"""

    with open('example_templates.fasta', 'w') as f:
        f.write(fasta_content)

    print("Created example files:")
    print("  - example_cutsites.tsv")
    print("  - example_templates.fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Design primers for CRISPR Cas9 cut site sequencing using Primer3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with TSV and FASTA files
  python crispr_primer_designer.py -c cutsites.tsv -t templates.fasta

  # Specify custom product size range for different Illumina platforms
  python crispr_primer_designer.py -c cutsites.tsv -t templates.fasta --product-size "150-300"

  # Create example files for testing
  python crispr_primer_designer.py --create-examples

  # Generate input only (don't run Primer3)
  python crispr_primer_designer.py -c cutsites.tsv -t templates.fasta --input-only
        """
    )

    parser.add_argument('-c', '--cutsites',
                        help='TSV file with cut site information')
    parser.add_argument('--grna',
                        help='Directly provide gRNA sequence (instead of -c TSV)')
    parser.add_argument('--location',
                        help='Directly provide genomic location (e.g. chr15:1000-2000)')
    parser.add_argument('--grna-name', default='gRNA',
                        help='Name for the gRNA (default: gRNA)')
    parser.add_argument('-t', '--templates',
                        help='FASTA file with template sequences (can contain multiple sequences)')
    parser.add_argument('-g', '--genome',
                        help='Path to genome FASTA file (e.g. hg38.fa) for automatic template extraction')
    parser.add_argument('-f', '--flank', type=int, default=250,
                        help='Flanking distance around cut site when extracting from genome (default: 250bp)')
    parser.add_argument('-p', '--prefix', default='primers',
                        help='Output file prefix (default: primers)')
    parser.add_argument('-o', '--output_dir',
                        help='Directory to put output files')
    parser.add_argument('--primer3-path', default='primer3_core',
                        help='Path to primer3_core executable')

    # Primer3 parameters
    parser.add_argument('--product-size', default='200-450',
                        help='Product size range (default: 200-450 for Illumina)')
    parser.add_argument('--primer-min-size', type=int, default=18,
                        help='Minimum primer size (default: 18)')
    parser.add_argument('--primer-opt-size', type=int, default=20,
                        help='Optimal primer size (default: 20)')
    parser.add_argument('--primer-max-size', type=int, default=25,
                        help='Maximum primer size (default: 25)')
    parser.add_argument('--min-tm', type=float, default=57.0,
                        help='Minimum melting temperature (default: 57.0)')
    parser.add_argument('--opt-tm', type=float, default=60.0,
                        help='Optimal melting temperature (default: 60.0)')
    parser.add_argument('--max-tm', type=float, default=63.0,
                        help='Maximum melting temperature (default: 63.0)')
    parser.add_argument('--num-return', type=int, default=5,
                        help='Number of primer pairs to return (default: 5)')
    parser.add_argument('--target-distance', type=int, default=50,
                        help='Minimum distance from cut site to primer (default: 50bp)')
    parser.add_argument('--max-poly-x', type=int, default=3,
                        help='Maximum identical nucleotides in a row (default: 3, removes 4 or more)')
    parser.add_argument('--check-uniqueness', action='store_true',
                        help='Check that primers are unique in the genome (requires --genome)')
    parser.add_argument('--list-locations', action='store_true',
                        help='List all genomic locations for each primer (requires --genome)')

    # Other options
    parser.add_argument('--create-examples', action='store_true',
                        help='Create example input files')
    parser.add_argument('--input-only', action='store_true',
                        help='Only generate Primer3 input, do not run Primer3')
    parser.add_argument('--config-json',
                        help='JSON file with Primer3 configuration')

    args = parser.parse_args()

    # Create examples if requested
    if args.create_examples:
        create_example_files()
        return

    # Check required arguments
    if not args.cutsites and not args.grna and not args.location:
        parser.error("Either --cutsites, --grna, or --location must be provided")

    if not args.templates and not args.genome:
        parser.error("Either --templates or --genome must be provided")

    if (args.check_uniqueness or args.list_locations) and not args.genome:
        parser.error("--check-uniqueness and --list-locations require --genome to be provided")

    # Create configuration
    if args.config_json:
        with open(args.config_json) as f:
            config_dict = json.load(f)
        config = Primer3Config(**config_dict)
    else:
        config = Primer3Config(
            product_size_range=args.product_size,
            primer_min_size=args.primer_min_size,
            primer_opt_size=args.primer_opt_size,
            primer_max_size=args.primer_max_size,
            primer_min_tm=args.min_tm,
            primer_opt_tm=args.opt_tm,
            primer_max_tm=args.max_tm,
            num_return=args.num_return,
            target_flank_min=args.target_distance,
            primer_max_poly_x=args.max_poly_x
        )

    # Initialize designer
    designer = CRISPRPrimerDesigner(config=config, primer3_path=args.primer3_path)

    # Check output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = None

    # Load cut sites
    if args.grna:
        print(f"Using provided gRNA sequence: {args.grna}")
        cutsites = [CutSite(
            target_gene=args.grna_name,
            grna_name=args.grna_name,
            sequence=args.grna,
            location=args.location or ""
        )]
    elif args.location:
        print(f"Using provided genomic location: {args.location}")
        cutsites = [CutSite(
            target_gene=args.grna_name,
            grna_name=args.grna_name,
            sequence="",
            location=args.location
        )]
    else:
        print(f"Loading cut sites from: {args.cutsites}")
        cutsites = designer.load_cutsites_from_tsv(args.cutsites)
    print(f"  Loaded {len(cutsites)} cut site(s)")

    # Load or generate templates
    templates = []
    if args.templates:
        print(f"Loading templates from: {args.templates}")
        templates = designer.load_templates_from_fasta(args.templates)
        print(f"  Loaded {len(templates)} templates")

    if not templates and args.genome:
        print(f"Extracting genomic templates from: {args.genome} (flank={args.flank}bp)")
        templates = designer.generate_templates_from_genome(cutsites, args.genome, args.flank)
        print(f"  Generated {len(templates)} templates")

        # Save generated templates for reference
        if templates:
            templates_filename = f"{args.prefix}_templates.fasta"
            templates_path = output_dir / templates_filename if output_dir else Path(templates_filename)
            with open(templates_path, 'w') as f:
                for t in templates:
                    f.write(f">{t.name} {t.chromosome}:{t.start_position}-{t.end_position}\n{t.sequence}\n")
            print(f"  Saved templates to: {templates_path}")

    if args.input_only:
        # Just generate the input file
        filename = f"{args.prefix}_input.txt"
        if output_dir:
            input_file = output_dir / filename
        else:
            input_file = Path(filename)

        designer.generate_all_inputs(cutsites, templates, str(input_file))
        print(f"\nGenerated Primer3 input file: {input_file}")
        print("Run Primer3 manually with:")
        print(f"  primer3_core < {input_file} > {input_file.parent / (args.prefix + '_output.txt')}")
    else:
        # Full pipeline
        results = designer.design_primers(cutsites, templates, args.prefix, args.output_dir)

        # Optional uniqueness check
        if args.check_uniqueness and args.genome:
            print("\nVerifying primer uniqueness across the genome...")
            results = designer.filter_unique_primers(results, args.genome)

            # Save the filtered results again to ensure the TSV is updated
            out_file = f"{args.prefix}_primers.tsv"
            snap_file = f"{args.prefix}_snapgene.tsv"
            bulk_file = f"{args.prefix}_bulk_order.tsv"
            if args.output_dir:
                out_file = str(Path(args.output_dir) / out_file)
                snap_file = str(Path(args.output_dir) / snap_file)
                bulk_file = str(Path(args.output_dir) / bulk_file)
            
            designer._save_results_tsv(results, out_file)
            designer._save_snapgene_tsv(results, snap_file)
            designer._save_bulk_order_tsv(results, bulk_file)

        # Optional location reporting
        if args.list_locations and args.genome:
            print("\nFinding all genomic locations for designed primers...")
            loc_file = f"{args.prefix}_locations.tsv"
            if args.output_dir:
                loc_file = str(Path(args.output_dir) / loc_file)
            designer.report_primer_locations(results, args.genome, loc_file)

        # Print summary
        print("\n" + "="*60)
        print("RESULTS SUMMARY")
        print("="*60)

        total_pairs = 0
        for result in results:
            if result.errors:
                print(f"\n{result.sequence_id}: ERRORS - {'; '.join(result.errors)}")
            elif result.num_returned == 0:
                print(f"\n{result.sequence_id}: No primers found")
                if result.pair_explain:
                    print(f"  Explanation: {result.pair_explain}")
            else:
                print(f"\n{result.sequence_id}: {result.num_returned} primer pair(s)")
                total_pairs += result.num_returned
                for pair in result.primer_pairs[:1]:  # Show first pair
                    print(f"  Best pair (penalty: {pair.pair_penalty:.2f}):")
                    print(f"    Left:  {pair.left_sequence} (Tm={pair.left_tm:.1f}°C)")
                    print(f"    Right: {pair.right_sequence} (Tm={pair.right_tm:.1f}°C)")
                    print(f"    Product size: {pair.product_size}bp")

        print(f"\nTotal primer pairs designed: {total_pairs}")


if __name__ == "__main__":
    main()
