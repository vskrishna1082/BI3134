#! /usr/bin/env python3
import scoreseq

seq_a = "ARNDC"
seq_b = "ACKMP"
substitution_matrix = "Protein_BLOSUM62.csv"

seq = scoreseq.SequencePair(seq_a, seq_b, substitution_matrix)
matrix = seq.smith_watermann()
print(matrix)
