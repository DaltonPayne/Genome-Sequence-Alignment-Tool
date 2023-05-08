#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner

def smith_waterman_alignment(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    
    alignments = aligner.align(seq1, seq2)
    
    return alignments

def read_fasta_file(filename):
    with open(filename, 'r') as f:
        sequences = [str(record.seq) for record in SeqIO.parse(f, 'fasta')]
    
    return sequences

def main():
    if len(sys.argv) != 3:
        print("Usage: python genome_alignment_tool.py <fasta_file1> <fasta_file2>")
        sys.exit(1)

    fasta_file1 = sys.argv[1]
    fasta_file2 = sys.argv[2]

    seq1 = read_fasta_file(fasta_file1)
    seq2 = read_fasta_file(fasta_file2)

    alignments = smith_waterman_alignment(seq1, seq2)

    for alignment in alignments:
        print(alignment)

if __name__ == "__main__":
    main()
