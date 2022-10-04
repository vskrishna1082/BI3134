#! /usr/bin/env python3
# Align nucleotide sequences with dynamic programming
import timeit
import numpy as np

# Define a Substitution Matrix and Gap Penality
subs_matrix = np.array([[2,0,1,0],[0,2,0,1],[1,0,2,0],[0,1,0,2]]) 
gp = 4

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
    t_m_tr = np.zeros(shape=(ls1)) # tracebackMatrix topRow
    st_m_tr = np.r_[[[s_m_tr]],[[t_m_tr]]] # concatenate both

    s_m_fc = np.arange(0, -gp*(ls2+1), -gp).reshape(-1,1) # sM firstClmn
    t_m_fc = np.zeros(shape=(ls2+1)).reshape(-1,1) # tM firstClmn
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

seq_a= "ATATAAGAGCTCAGA"
seq_b= "TATAATCAGCAGTATATAT"
nseq_a = numSeq(seq_a)
nseq_b = numSeq(seq_b)
maxScoreBack(-6, 0, 0, 2)
def mytestfunc():
    stmatrix=dyn_prog_matrix(nseq_a,nseq_b)
mytestfunc()
print(timeit.timeit(mytestfunc, number=1))