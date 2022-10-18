#! /usr/bin/env python3
# Generate a score and traceback matrix from two sequences

import numpy as np
import symloader

class SequencePair: 
    """
    Takes two string of sequences and a substitution_data object as input 
    and returns a dynamic programming matrix computed using NeedleMann Wunsch
    (global alignment) algorithm.
    
    Provides 'scoreMatrix', 'traceMatrix', 'maxScore' and 'sequencePair'.
    """
    def __init__(self, seq_a, seq_b, substitution_matrix):
        data = symloader.LoadSymbols(substitution_matrix)
        library = data.symLibrary

        # numerize sequence string
        def numSeq(x):
            def nuNum(x):
                for idx, i in enumerate(library):
                    if x==i:
                        return idx
                        break
            return list(map(nuNum, [*x.upper()]))

        self.seq_a = seq_a
        self.seq_b = seq_b
        self.num_seq1 = numSeq(seq_a)
        self.num_seq2 = numSeq(seq_b)
        self.subs_matrix = data.subScores
        self.gap_pen = data.gapPenality

    # construct a matrix of required shape (0s for to-be-computed positions)
    def make_matrix(self):
        seq1 = self.num_seq1
        seq2 = self.num_seq2
        gp = self.gap_pen

        ls1 = len(seq1)
        ls2 = len(seq2)

        # create array of size 2xNxM
        dyn_matrix = np.zeros(shape=(2, ls1, ls2))
        # generate the first columns
        if gp != 0: score_firstcol = np.arange(-gp, -gp*(ls1+1), -gp).reshape(-1,1)
        else: score_firstcol = np.zeros(shape = ls1).reshape(-1,1)
        trace_firstcol = 2*np.ones(shape=(ls1)).reshape(-1,1)
        firstcol = np.r_[[score_firstcol],[trace_firstcol]]

        # generate the top rows
        if gp != 0: score_toprow = np.arange(0, -gp*(ls2+1), -gp)
        else: score_toprow = np.zeros(shape = ls2+1)
        trace_toprow = np.ones(shape=(ls2+1))
        trace_toprow[0]=0
        toprow = np.r_[[[score_toprow],[trace_toprow]]]

        # concatenate top row, first column and core matrix
        firstcol_core = np.concatenate((firstcol, dyn_matrix), axis=2)
        dyn_matrix = np.concatenate((toprow, firstcol_core), axis = 1)
        dyn_matrix = dyn_matrix.astype(int)
        return dyn_matrix

    def needleman_wunsch(self):
    # quick algorithm to get maximum score etc.
        def maxScoreBack(score_NW, score_N, score_W, subs_score):
            gp = self.gap_pen
            tBack=[0,0,0]
            nscore_NW = score_NW + subs_score
            nscore_N = score_N -gp
            nscore_W = score_W -gp
            nscores=[nscore_NW,nscore_N,nscore_W]
            maxNscore=max(nscores)# get max nscore
            def makeOctal(x):
                bi_str=''.join(str(i) for i in x)
                return int(bi_str, 2)
            for idx, i in enumerate(nscores):
                if i == maxNscore:
                    tBack[idx]=1
            tBack=makeOctal(tBack)
            return [maxNscore,tBack]
        st_m = self.make_matrix()
        for k in list(np.arange(len(self.num_seq1))+1):
            for j in list(np.arange(len(self.num_seq2))+1):
                result=maxScoreBack(st_m[0][j-1][k-1],st_m[0][j-1][k],st_m[0][j][k-1],self.subs_matrix[self.num_seq1[j-1],self.num_seq2[k-1]])
                st_m[0][j][k] = result[0]
                st_m[1][j][k] = result[1]
        return(st_m[0])

    def smith_watermann(self):
    # quick algorithm to get maximum score etc.
        def maxScoreBack(score_NW, score_N, score_W, subs_score):
            gp = self.gap_pen
            tBack=[0,0,0]
            nscore_NW = score_NW + subs_score
            nscore_N = score_N -gp
            nscore_W = score_W -gp
            nscores=[nscore_NW,nscore_N,nscore_W, 0]
            maxNscore=max(nscores)# get max nscore
            def makeOctal(x):
                bi_str=''.join(str(i) for i in x)
                return int(bi_str, 2)
            if maxNscore != 0:
                for idx, i in enumerate(nscores[:-1]):
                    if i == maxNscore:
                        tBack[idx]=1
                tBack=makeOctal(tBack)
            else: tBack = 0
            return [maxNscore,tBack]
        st_m = self.make_matrix()
        st_m[:,0] = 0
        st_m[:,:,0] = 0
        for k in list(np.arange(len(self.num_seq1))+1):
            for j in list(np.arange(len(self.num_seq2))+1):
                result=maxScoreBack(st_m[0][j-1][k-1],st_m[0][j-1][k],st_m[0][j][k-1],self.subs_matrix[self.num_seq1[j-1],self.num_seq2[k-1]])
                st_m[0][j][k] = result[0]
                st_m[1][j][k] = result[1]
        return(st_m)
