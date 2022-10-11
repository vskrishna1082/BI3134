#! /usr/bin/env python3
# Symbol Loader for sequence alignment (dynamic programming)
import pandas as pd
import os

class LoadSymbols:
    """Get substitution scores and gap penality from file"""
    def __init__(self, matrix_fileName):
        matrix_file = os.path.join('substitution_matrices', matrix_fileName)
        f = open(matrix_file, 'r')
        f.seek(14)
        df = pd.read_csv(matrix_file,
                         sep=',',
                         comment='#',
                         index_col=0,
                         skipinitialspace=True
                         )
        self.gapPenality = int(f.read(1))
        self.symLibrary = list(df.columns)
        self.subScores = df.to_numpy()
