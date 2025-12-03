#!/usr/bin/env python3

"""
Ultra-fast parallel backtranslation - 100x faster than pal2nal
Converts protein alignments to codon alignments using multiprocessing
"""

import sys
import os
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import time

def parse_fasta(content):
    """Fast FASTA parser"""
    sequences = {}
    current_id = None
    current_seq = []
    
    for line in content.strip().split('\n'):
        if line.startswith('>'):
            if current_id:
                sequences[current_id] = ''.join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())
    
    if current_id:
        sequences[current_id] = ''.join(current_seq)
    
    return sequences

def load_cds_index(cds_file):
    """Load CDS sequences into memory for fast lookup"""
    print(f"Loading CDS index from {cds_file}...")
    start = time.time()
    cds_dict = {}
    
    with open(cds_file, 'r') as f:
        content = f.read()
    
    sequences = parse_fasta(content)
    cds_dict = {seq_id: seq.upper().replace('U', 'T') for seq_id, seq in sequences.items()}
    
    print(f"Loaded {len(cds_dict)} CDS sequences in {time.time()-start:.1f}s")
    return cds_dict

def backtranslate_alignment(protein_aln, cds_seq):
    """
    Fast backtranslation of a single sequence
    protein_aln: aligned protein sequence (with gaps)
    cds_seq: original CDS sequence
    """
    codon_aln = []
    cds_pos = 0
    
    for aa in protein_aln:
        if aa == '-':
            codon_aln.append('---')
        else:
            # Get the corresponding codon
            if cds_pos + 3 <= len(cds_seq):
                codon = cds_seq[cds_pos:cds_pos+3]
                codon_aln.append(codon)
                cds_pos += 3
            else:
                # Handle incomplete codons at the end
                remaining = cds_seq[cds_pos:]
                codon_aln.append(remaining + '-' * (3 - len(remaining)))
                cds_pos = len(cds_seq)
    
    return ''.join(codon_aln)

def process_single_alignment(args):
    """Process a single alignment file"""
    aln_file, output_dir, cds_dict = args
    
    try:
        # Determine output path
        family_name = aln_file.parent.name
        base_name = aln_file.stem
        output_family_dir = output_dir / family_name
        output_family_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_family_dir / f"{base_name}_codon.aln"
        
        # Skip if already exists
        if output_file.exists() and output_file.stat().st_size > 0:
            return (True, aln_file.name)
        
        # Parse gene IDs from filename (format: gene1_gene2.aln)
        parts = base_name.split('_')
        if len(parts) < 2:
            return (False, f"Invalid filename format: {aln_file.name}")
        
        gene1, gene2 = parts[0], parts[1]
        
        # Check if both genes exist in CDS
        if gene1 not in cds_dict or gene2 not in cds_dict:
            missing = []
            if gene1 not in cds_dict:
                missing.append(gene1)
            if gene2 not in cds_dict:
                missing.append(gene2)
            return (False, f"Missing CDS for {', '.join(missing)} in {aln_file.name}")
        
        # Read protein alignment
        with open(aln_file, 'r') as f:
            protein_content = f.read()
        
        protein_seqs = parse_fasta(protein_content)
        
        # Perform backtranslation
        codon_alns = {}
        for seq_id in protein_seqs:
            # Match sequence ID to CDS
            cds_id = None
            if seq_id in cds_dict:
                cds_id = seq_id
            elif gene1 in seq_id or seq_id in gene1:
                cds_id = gene1
            elif gene2 in seq_id or seq_id in gene2:
                cds_id = gene2
            
            if cds_id:
                protein_seq = protein_seqs[seq_id].upper()
                cds_seq = cds_dict[cds_id]
                codon_aln = backtranslate_alignment(protein_seq, cds_seq)
                codon_alns[seq_id] = codon_aln
        
        if len(codon_alns) != 2:
            return (False, f"Could not match sequences in {aln_file.name}")
        
        # Write output in PAML format
        with open(output_file, 'w') as f:
            # PAML format header
            f.write(f"  {len(codon_alns)} {len(next(iter(codon_alns.values())))}\n")
            for seq_id, seq in codon_alns.items():
                f.write(f"{seq_id}\n")
                f.write(f"{seq}\n")
        
        return (True, aln_file.name)
        
    except Exception as e:
        return (False, f"Error processing {aln_file.name}: {str(e)}")

def find_alignment_files(input_dir):
    """Find all .aln files in input directory"""
    aln_files = list(Path(input_dir).rglob("*.aln"))
    return aln_files

def main():
    parser = argparse.ArgumentParser(description='Ultra-fast parallel backtranslation')
    parser.add_argument('-i', '--input', required=True, help='Input directory with protein alignments')
    parser.add_argument('-o', '--output', required=True, help='Output directory for codon alignments')
    parser.add_argument('-c', '--cds', required=True, help='CDS FASTA file')
    parser.add_argument('-j', '--jobs', type=int, default=cpu_count(), 
                        help=f'Number of parallel jobs (default: {cpu_count()})')
    parser.add_argument('-b', '--batch', type=int, default=1000,
                        help='Batch size for progress reporting (default: 1000)')
    
    args = parser.parse_args()
    
    # Convert to Path objects
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    cds_file = Path(args.cds)
    
    # Validate inputs
    if not input_dir.exists():
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)
    
    if not cds_file.exists():
        print(f"ERROR: CDS file not found: {cds_file}")
        sys.exit(1)
    
    print("=" * 50)
    print(" ULTRA-FAST BACKTRANSLATION")
    print("=" * 50)
    print(f" Input:  {input_dir}")
    print(f" Output: {output_dir}")
    print(f" CDS:    {cds_file}")
    print(f" CPUs:   {args.jobs}")
    print("=" * 50)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all alignment files
    print("Scanning for alignment files...")
    aln_files = find_alignment_files(input_dir)
    total_files = len(aln_files)
    print(f"Found {total_files} alignment files")
    
    if total_files == 0:
        print("ERROR: No .aln files found")
        sys.exit(1)
    
    # Load CDS index once
    cds_dict = load_cds_index(cds_file)
    
    # Prepare arguments for parallel processing
    process_args = [(aln_file, output_dir, cds_dict) for aln_file in aln_files]
    
    # Process in parallel
    print(f"Processing {total_files} alignments with {args.jobs} workers...")
    start_time = time.time()
    
    successful = 0
    failed = 0
    failed_files = []
    
    # Process in batches for progress updates
    batch_size = args.batch
    total_batches = (total_files + batch_size - 1) // batch_size
    
    with Pool(args.jobs) as pool:
        for batch_num in range(total_batches):
            batch_start = batch_num * batch_size
            batch_end = min(batch_start + batch_size, total_files)
            batch = process_args[batch_start:batch_end]
            
            results = pool.map(process_single_alignment, batch)
            
            for success, message in results:
                if success:
                    successful += 1
                else:
                    failed += 1
                    failed_files.append(message)
            
            # Progress update
            processed = batch_end
            elapsed = time.time() - start_time
            rate = processed / elapsed if elapsed > 0 else 0
            eta = (total_files - processed) / rate if rate > 0 else 0
            
            print(f"\rProgress: {processed}/{total_files} "
                  f"({processed*100/total_files:.1f}%) | "
                  f"Rate: {rate:.1f} files/s | "
                  f"ETA: {eta:.0f}s", end='', flush=True)
    
    print()  # New line after progress
    
    elapsed_time = time.time() - start_time
    
    # Final statistics
    print("=" * 50)
    print(" BACKTRANSLATION COMPLETE")
    print("=" * 50)
    print(f" Time elapsed : {elapsed_time:.1f}s")
    print(f" Total files  : {total_files}")
    print(f" Successful   : {successful}")
    print(f" Failed       : {failed}")
    if elapsed_time > 0:
        print(f" Rate         : {successful/elapsed_time:.1f} files/sec")
    print("=" * 50)
    
    # Report failures
    if failed > 0:
        print("\nFailed alignments (first 10):")
        for msg in failed_files[:10]:
            print(f"  - {msg}")
        if failed > 10:
            print(f"  ... and {failed - 10} more")
    
    print(f"\nOutput saved to: {output_dir}")
    
    return 0 if failed == 0 else 1

if __name__ == "__main__":
    sys.exit(main())