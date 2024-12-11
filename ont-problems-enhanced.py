#!/usr/bin/env python3

import argparse
from datetime import datetime
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from pathlib import Path
import logging
import re
from collections import Counter
from typing import List, Dict, Set, Tuple

class KmerAnalyzer:
    def __init__(self, k_lengths: List[int], end_window: int = 1000):
        self.k_lengths = k_lengths
        self.end_window = end_window
        self.kmer_counts = {k: Counter() for k in k_lengths}
        self.end_kmers = {}
        
    def find_kmers(self, sequence: str, contig_id: str) -> List[Dict]:
        """Find all k-mers and their positions in the sequence."""
        kmer_positions = []
        seq_len = len(sequence)
        
        for k in self.k_lengths:
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                self.kmer_counts[k][kmer] += 1
                
                # Check if k-mer is near sequence ends
                distance_to_start = i
                distance_to_end = seq_len - (i + k)
                is_near_end = (distance_to_start <= self.end_window or 
                             distance_to_end <= self.end_window)
                
                if is_near_end:
                    if contig_id not in self.end_kmers:
                        self.end_kmers[contig_id] = Counter()
                    self.end_kmers[contig_id][kmer] += 1
                
                # Get reverse complement
                rc_kmer = str(Seq(kmer).reverse_complement())
                
                kmer_positions.append({
                    'contig': contig_id,
                    'start': i,
                    'end': i + k,
                    'kmer': kmer,
                    'length': k,
                    'is_reverse': False,
                    'near_end': is_near_end,
                    'distance_to_nearest_end': min(distance_to_start, distance_to_end)
                })
                
        return kmer_positions

class ProblemRegionDetector:
    def __init__(self, min_homopolymer_length=4, min_repeat_length=3):
        self.min_homopolymer_length = min_homopolymer_length
        self.min_repeat_length = min_repeat_length
        
    def find_homopolymers(self, sequence: str, contig_id: str) -> List[Dict]:
        """Find homopolymer runs in sequence."""
        homopolymers = []
        pattern = f'([ATCG]){{{self.min_homopolymer_length},}}'
        
        for match in re.finditer(pattern, str(sequence)):
            homopolymers.append({
                'contig': contig_id,
                'start': match.start(),
                'end': match.end(),
                'length': match.end() - match.start(),
                'sequence': match.group(),
                'type': 'homopolymer'
            })
        return homopolymers
    
    def find_repeats(self, sequence: str, contig_id: str) -> List[Dict]:
        """Find di- and tri-nucleotide repeats."""
        repeats = []
        patterns = [
            (f'([ATCG]{{2}}){{{self.min_repeat_length},}}', 'dimer'),
            (f'([ATCG]{{3}}){{{self.min_repeat_length},}}', 'trimer')
        ]
        
        for pattern, repeat_type in patterns:
            for match in re.finditer(pattern, str(sequence)):
                if len(set(match.group())) < len(match.group()):
                    repeats.append({
                        'contig': contig_id,
                        'start': match.start(),
                        'end': match.end(),
                        'length': match.end() - match.start(),
                        'sequence': match.group(),
                        'type': f'{repeat_type}_repeat'
                    })
        return repeats

def write_bed_file(regions: List[Dict], output_path: str):
    """Write regions to BED format (0-based)."""
    with open(output_path, 'w') as f:
        for region in regions:
            f.write(f"{region['contig']}\t{region['start']}\t{region['end']}\t"
                   f"{region.get('type', 'kmer')}_{region.get('sequence', region.get('kmer', ''))}\t"
                   f"{region.get('length', 0)}\t{'+' if not region.get('is_reverse', False) else '-'}\n")

def write_gff_file(regions: List[Dict], output_path: str):
    """Write regions to GFF format (1-based)."""
    with open(output_path, 'w') as f:
        f.write("##gff-version 3\n")
        for region in regions:
            attributes = f"ID={region.get('type', 'kmer')}_{region['start']};"\
                        f"Sequence={region.get('sequence', region.get('kmer', ''))}"
            f.write(f"{region['contig']}\tONT_problems\t{region.get('type', 'kmer')}\t"
                   f"{region['start'] + 1}\t{region['end']}\t.\t"
                   f"{'+' if not region.get('is_reverse', False) else '-'}\t.\t{attributes}\n")

def write_gtf_file(regions: List[Dict], output_path: str):
    """Write regions to GTF format (1-based)."""
    with open(output_path, 'w') as f:
        for region in regions:
            attributes = f"sequence \"{region.get('sequence', region.get('kmer', ''))}\";"
            f.write(f"{region['contig']}\tONT_problems\t{region.get('type', 'kmer')}\t"
                   f"{region['start'] + 1}\t{region['end']}\t.\t"
                   f"{'+' if not region.get('is_reverse', False) else '-'}\t.\t{attributes}\n")

def process_sequences(file_path: str, k_lengths: List[int], output_format: str = 'bed'):
    """Process sequences and identify problematic regions and k-mers."""
    detector = ProblemRegionDetector()
    kmer_analyzer = KmerAnalyzer(k_lengths)
    
    all_regions = []
    open_func = gzip.open if str(file_path).endswith('.gz') else open
    
    with open_func(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            
            # Find problematic regions
            homopolymers = detector.find_homopolymers(sequence, record.id)
            repeats = detector.find_repeats(sequence, record.id)
            kmers = kmer_analyzer.find_kmers(sequence, record.id)
            
            all_regions.extend(homopolymers + repeats + kmers)
    
    return all_regions, kmer_analyzer

def main():
    parser = argparse.ArgumentParser(description='Detect problematic regions for ONT sequencing')
    parser.add_argument('--input', required=True, help='Input FASTA file (can be gzipped)')
    parser.add_argument('--output-prefix', default='ont_problems', help='Output file prefix')
    parser.add_argument('--kmer-lengths', type=int, nargs='+', default=[6, 7],
                       help='K-mer lengths to analyze')
    parser.add_argument('--output-format', choices=['bed', 'gff', 'gtf'], default='bed',
                       help='Output format for regions')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    # Process sequences
    logger.info(f"Processing sequences from {args.input}")
    all_regions, kmer_analyzer = process_sequences(args.input, args.kmer_lengths)
    
    # Create output directory
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Write regions file in requested format
    output_path = f"{args.output_prefix}_regions.{args.output_format}"
    if args.output_format == 'bed':
        write_bed_file(all_regions, output_path)
    elif args.output_format == 'gff':
        write_gff_file(all_regions, output_path)
    else:  # gtf
        write_gtf_file(all_regions, output_path)
    
    # Save problems to CSV
    problems_df = pd.DataFrame(all_regions)
    problems_df.to_csv(f"{args.output_prefix}_problems.csv", index=False)
    
    # Generate metadata summary
    problem_counts = problems_df['type'].value_counts().to_dict()
    metadata = {
        'date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'input_file': args.input,
        'total_regions': len(all_regions),
        'problem_counts': problem_counts,
        'kmer_stats': {
            k: dict(counts.most_common(5)) 
            for k, counts in kmer_analyzer.kmer_counts.items()
        },
        'end_kmers': {
            contig: dict(counts.most_common(5))
            for contig, counts in kmer_analyzer.end_kmers.items()
        }
    }
    
    with open(f"{args.output_prefix}_metadata.txt", 'w') as meta_out:
        for key, value in metadata.items():
            meta_out.write(f"{key}:\n{value}\n\n")
    
    logger.info("Processing complete")
    
    # Display summary
    print("\nTop k-mers at contig ends:")
    for contig, counts in kmer_analyzer.end_kmers.items():
        print(f"\n{contig}:")
        for kmer, count in counts.most_common(5):
            print(f"  {kmer}: {count}")
    
    print("\nProblem region counts by type:")
    for problem_type, count in problem_counts.items():
        print(f"  {problem_type}: {count}")

if __name__ == "__main__":
    main()
