#! /usr/bin/env python3
# Align nucleotide sequences with global dynamic programming

import timeit
import numpy as np
from collections import deque

# Provide a sequence string
seq_a= "ATGTATGCTGCTGCATTGTCTGATCG"
seq_b= "AATGTAATGCGCTGTCGATCGTAGCT"
alignments=[]

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
seq1 = numSeq(seq_a)
seq2 = numSeq(seq_b)
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
def trace_alignment(st_matrix):
    a_seq1 = deque() # aligned sequences as lists
    a_seq2 = deque()
    x = np.size(st_matrix,axis=2)-1
    y =  np.size(st_matrix,axis=1)-1
    curr_tback = st_matrix[1][y][x]

    def move_diagonal(al_seq1,al_seq2, pos_x, pos_y):
        al_seq1.appendleft(nuLib[seq1[pos_x-1]])
        al_seq2.appendleft(nuLib[seq2[pos_y-1]])
        pos_y+=-1
        pos_x+=-1
        tback = st_matrix[1][pos_y][pos_x]
        return (tback,pos_x,pos_y)

    def move_left(al_seq1,al_seq2, pos_x, pos_y):
        pos_x+=-1
        al_seq1.appendleft(nuLib[seq1[pos_x]])
        al_seq2.appendleft('-')
        tback = st_matrix[1][pos_y][pos_x]
        return (tback,pos_x,pos_y)

    def move_up(al_seq1,al_seq2, pos_x, pos_y):
        pos_y+=-1
        al_seq1.appendleft('-')
        al_seq2.appendleft(nuLib[seq2[pos_y]])
        tback = st_matrix[1][pos_y][pos_x]
        return (tback,pos_x,pos_y)

    def branch_diagonal():
        branch = [a_seq1.copy(),a_seq2.copy()]
        branch[0].appendleft(nuLib[seq1[x-1]])
        branch[1].appendleft(nuLib[seq2[y-1]])
        sub_align = trace_alignment(st_matrix[:,:y,:x])
        branch[0]=sub_align[0]+branch[0] 
        branch[1]=sub_align[1]+branch[1]
        alignments.append(tuple(branch))

    def branch_left():
        branch = [a_seq1.copy(),a_seq2.copy()]
        branch[0].appendleft(nuLib[seq1[x-1]])
        branch[1].appendleft('-')
        sub_align = trace_alignment(st_matrix[:,:y+1,:x])
        branch[0]=sub_align[0]+branch[0] 
        branch[1]=sub_align[1]+branch[1]
        alignments.append(tuple(branch))

    while curr_tback != 0:
#        print(f"{x},{y} tb= {curr_tback}") # debugging
        if curr_tback == 4:
            (curr_tback,x,y) = move_diagonal(a_seq1,a_seq2,x,y)

        elif curr_tback == 1: 
            (curr_tback,x,y) = move_left(a_seq1,a_seq2,x,y)

        elif curr_tback == 2:
            (curr_tback,x,y) = move_up(a_seq1,a_seq2,x,y)

        elif curr_tback == 3:
            branch_left()
            (curr_tback,x,y) = move_up(a_seq1,a_seq2,x,y)

        elif curr_tback == 5:
            branch_diagonal()
            (curr_tback,x,y) = move_left(a_seq1,a_seq2,x,y)

        elif curr_tback == 6:
            branch_diagonal()
            (curr_tback,x,y) = move_up(a_seq1,a_seq2,x,y)

        elif curr_tback == 7:
            branch_diagonal()
            branch_left()
            (curr_tback,x,y) = move_up(a_seq1,a_seq2,x,y)
    return (a_seq1,a_seq2)

def mytestfunc():
    stmatrix=dyn_prog_matrix(seq1,seq2)
    aligned=trace_alignment(stmatrix)
    alignments.append(aligned)
    for idx, i in enumerate(alignments):
        print('\nSequence '+str(idx+1)+':\n'+''.join(i[0])+'\n'+''.join(i[1]))
    print("Score: "+str(stmatrix[0][len(seq2)][len(seq1)]))
print(f"Time taken: {timeit.timeit(mytestfunc, number=1)}")
