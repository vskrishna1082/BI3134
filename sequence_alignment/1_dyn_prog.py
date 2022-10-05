#! /usr/bin/env python3
# Align nucleotide sequences with global dynamic programming
# (!)returns only one of the best scoring alignments

import timeit
import numpy as np
from collections import deque

# Provide a sequence string
seq_a= "AATATTATATCCT"
seq_b= "ATgtacatcct"

# Define a Substitution Matrix and Gap Penality
subs_matrix = np.array([[2,0,1,0],[0,2,0,1],[1,0,2,0],[0,1,0,2]]) 
gp = 4 # matrix generation fails when gp=0

# Convert string of input to numbered list of nucleotides
nuLib=["A","T","G","C"] # Library of symbols
def numSeq(x):
    def nuNum(x):
        for idx, i in enumerate(nuLib):
            if x==i:
                return idx
                break # stop at first match
    return list(map(nuNum, [*x.upper()]))

# Define a function that maximizes score and returns (score,traceback)
def maxScoreBack(score_NW, score_N, score_W, subs_score):
    tBack=[0,0,0] # define empty traceback
    nscore_NW = score_NW + subs_score # get new scores
    nscore_N = score_N - gp
    nscore_W = score_W - gp
    nscores=[nscore_NW,nscore_N,nscore_W]
    maxNscore=max(nscores)# get max nscore
    for idx, i in enumerate(nscores):
        if i == maxNscore:
            tBack[idx]=1
    tBack=makeOctal(tBack)
    return [maxNscore,tBack]

# Convert a list of 3 bools into (fake) octal (fn:traceback)
def makeOctal(x):
    bi_str=''.join(str(i) for i in x)
    return int(bi_str, 2) # not an octal, but works

# Create the scoring and traceback matrices
def dyn_prog_matrix(seq1, seq2):
    ls1 = len(seq1)
    ls2 = len(seq2)
    # add top-row and first-column to array
    s_m_tr = np.arange(-gp, -gp*(ls1+1), -gp) # scoringMatrix topRow
    t_m_tr = np.ones(shape=(ls1)) # tracebackMatrix topRow
    st_m_tr = np.r_[[[s_m_tr]],[[t_m_tr]]] # concatenate both

    s_m_fc = np.arange(0, -gp*(ls2+1), -gp).reshape(-1,1) # sM firstClmn
    t_m_fc = 2*np.ones(shape=(ls2+1)).reshape(-1,1) # tM firstClmn
    t_m_fc[0]=0
    st_m_fc = np.r_[[s_m_fc],[t_m_fc]] # concatenate both

    # create random np array of size 2xNxM
    st_m= np.zeros(shape=(2, ls2,ls1)) # Score(NxM)xTrace(NxM)
    st_m = np.c_[st_m_fc, np.concatenate((st_m_tr, st_m),axis=1)] # final 2xN+1xM+1
    st_m = st_m.astype(int)

    # dynamic programming recursive relation
    for k in list(np.arange(len(seq1))+1):
        for j in list(np.arange(len(seq2))+1):
            result=maxScoreBack(st_m[0][j-1][k-1],st_m[0][j-1][k],st_m[0][j][k-1],subs_matrix[seq2[j-1],seq1[k-1]])
            st_m[0][j][k]=result[0]
            st_m[1][j][k]=result[1]
    return(st_m)

# Trace best alignment from scoring and traceack matrices
def trace_alignment(st_matrix, seq1, seq2):
    a_seq1 = deque() # aligned sequences as lists
    a_seq2 = deque()
    pos_x = len(seq1)
    pos_y = len(seq2)
    curr_tback = st_matrix[1][pos_y][pos_x]
    while curr_tback != 0:
        #print(f"{pos_x},{pos_y} tb= {curr_tback}") # debugging
        if curr_tback == 4:
            a_seq1.appendleft(nuLib[seq1[pos_x-1]])
            a_seq2.appendleft(nuLib[seq2[pos_y-1]])
            pos_y+=-1
            pos_x+=-1
            curr_tback = st_matrix[1][pos_y][pos_x]
        elif curr_tback == 1 or curr_tback == 5: # gap in second sequence
            pos_x+=-1
            a_seq1.insert(0, nuLib[seq1[pos_x]])
            a_seq2.insert(0, '-')
            curr_tback = st_matrix[1][pos_y][pos_x]
        elif curr_tback == 2 or curr_tback == 3 or curr_tback == 6 or curr_tback == 7: # gap in first sequence
            pos_y+=-1
            a_seq1.insert(0, '-')
            a_seq2.insert(0, nuLib[seq2[pos_y]])
            curr_tback = st_matrix[1][pos_y][pos_x]
    return ''.join(a_seq1)+"\n"+''.join(a_seq2)

def mytestfunc():
    nseq_a = numSeq(seq_a)
    nseq_b = numSeq(seq_b)
    stmatrix=dyn_prog_matrix(nseq_a,nseq_b)
    aligned=trace_alignment(stmatrix,nseq_a,nseq_b)
    print(aligned)
    print("Score: "+str(stmatrix[0][len(nseq_b)][len(nseq_a)]))
print(f"Time taken: {timeit.timeit(mytestfunc, number=1)}")
