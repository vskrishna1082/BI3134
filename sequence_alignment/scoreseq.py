#! /usr/bin/env python3
# Generate a score and traceback matrix from two sequences

import numpy as np
import symloader

class SequencePair: 
    """
    Functions and attributes for a pair of sequences to be aligned
    """
    def __init__(self, seq_a, seq_b, substitution_matrix):
        data = symloader.LoadSymbols(substitution_matrix)

        # numerize sequence string
        def numSeq(x):
            def nuNum(x):
                for idx, i in enumerate(data.symLibrary):
                    if x==i:
                        return idx
                        break
            return list(map(nuNum, [*x.upper()]))

        self.num_seq1 = numSeq(seq_a)
        self.num_seq2 = numSeq(seq_b)
        self.subs_matrix = data.subScores
        self.gp = data.gapPenality

    # construct a matrix of required shape (0s for to-be-computed positions)
    def make_matrix(self, method='global'):
        seq1 = self.num_seq1
        seq2 = self.num_seq2

        ls1 = len(seq1)
        ls2 = len(seq2)

        # create array of size 2xNxM
        dyn_matrix = np.zeros(shape=(2, ls1, ls2))
        # generate the first columns
        if self.gp != 0: score_firstcol = np.arange(-self.gp, -self.gp*(ls1+1), -self.gp).reshape(-1,1)
        else: score_firstcol = np.zeros(shape = ls1).reshape(-1,1)
        trace_firstcol = 2*np.ones(shape=(ls1)).reshape(-1,1)
        firstcol = np.r_[[score_firstcol],[trace_firstcol]]

        # generate the top rows
        if self.gp != 0: score_toprow = np.arange(0, -self.gp*(ls2+1), -self.gp)
        else: score_toprow = np.zeros(shape = ls2+1)
        trace_toprow = np.ones(shape=(ls2+1))
        trace_toprow[0]=0
        toprow = np.r_[[[score_toprow],[trace_toprow]]]

        # concatenate top row, first column and core matrix
        firstcol_core = np.concatenate((firstcol, dyn_matrix), axis=2)
        dyn_matrix = np.concatenate((toprow, firstcol_core), axis = 1)
        dyn_matrix = dyn_matrix.astype(int)
        if method == 'local':
            dyn_matrix[:,0] = 0
            dyn_matrix[:,:,0] = 0
        return dyn_matrix

    def align(self, method='global'):
    # quick algorithm to get maximum score etc.
        def maxScoreBack(score_NW, score_N, score_W, subs_score):
            tBack=[0,0,0]
            nscore_NW = score_NW + subs_score
            nscore_N = score_N -self.gp
            nscore_W = score_W -self.gp

            def makeOctal(x):
                bi_str=''.join(str(i) for i in x)
                return int(bi_str, 2)

            if method == 'local':
                nscores=[nscore_NW,nscore_N,nscore_W, 0]
            else:
                nscores=[nscore_NW,nscore_N,nscore_W]
            maxNscore=max(nscores)# get max nscore
            for idx, i in enumerate(nscores[:2]):
                if i == maxNscore and maxNscore != 0:
                    print(str(i)+'='+str(maxNscore))
                    tBack[idx]=1
                print(tBack)
            tBack=makeOctal(tBack)
            return [maxNscore,tBack]

        st_m = self.make_matrix(method=method)
        for k in list(np.arange(len(self.num_seq1))+1):
            for j in list(np.arange(len(self.num_seq2))+1):
                result=maxScoreBack(st_m[0][j-1][k-1],st_m[0][j-1][k],st_m[0][j][k-1],self.subs_matrix[self.num_seq1[j-1],self.num_seq2[k-1]])
                st_m[0][j][k] = result[0]
                st_m[1][j][k] = result[1]
        return(st_m)
