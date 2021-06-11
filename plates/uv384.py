import sys 
import os
import re 
from string import ascii_uppercase
import numpy as np
import pandas as pd 
from scipy.ndimage import gaussian_filter1d
from tqdm import  tqdm

import plates # ???

class Block:
    def __init__(self, data, control = None, concs = None, smiles = None):
        self.data = data
        self.control = control
        self.concs = concs
        self.smiles = smiles
    @property
    def df(self):
        return self.data
    @property
    def ctrl(self):
        return self.control
    @property
    def norm(self):
        df = self.df
        if self.control is not None:
            ctrl = self.control 
            df = df.sub(ctrl) # issue - crtl.shape != df.shape sometimes
            df = df.dropna(axis=1)
        df = df.sub(df.iloc[:,-1], axis=0) # anchor at 800 nm
        return df
    @property
    def diff(self):
        norm = self.norm
        return norm.sub(norm.iloc[0,:])
    @property
    def response(self):
        diff = self.diff
        response = diff.loc[:,390].abs() + diff.loc[:,420].abs()
        if self.concs is not None:
            assert len(self.concs) == len(response)
            response.index = self.concs
        return response
    @property
    def mm(self):
        assert self.concs is not None, 'Need some concentrations bud!'
        response = self.response
        return plates.analysis.MichaelisMenten(response, pd.Series(response.index))

    def report(self):
        pass



class UV384(plates.Plate):
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None):
        super().__init__(path, name, parser)
        self.concs = concs 
        if smiles is None:
            self.smiles = [None] * 24
        else:
            self.smiles = smiles
        if names is None:
            self.names = range(24)
        else:
            self.names = names
        self.control = None
    @property
    def blocknums(self):
        pass
    @property
    def blocks(self):
        if self.blocknums is not None:
            return [Block(self.locblock(i), 
                        concs = self.concs,
                        smiles = j) for i, j in zip(self.blocknums, self.smiles)]
    def locblock(self, num):
        pass

    def block(self, n, smiles = None):
        if smiles is None:
            if self.smiles is not None:
                smiles = self.smiles[n]
            else:
                smiles = None
        if self.control is not None:
            return Block(data = self.locblock(n),
                        control = self.control.locblock(n), 
                        concs = self.concs,
                        smiles = smiles)
        else:
            return Block(self.locblock(n), 
                    concs = self.concs,
                        smiles = smiles)

class UV384m1(UV384):
    # hand assay v1 (controls in each plate) plate obj
    # each col is 1 compound and contains a control in every other well
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None):
        super().__init__(path, name, concs, smiles, names, parser )
        # this is getting confusing
    @property
    def blocknums(self):
        return range(1, len(self.df) // 16)

    def locblock(self, num):
        col = self.col(num)
        return col[::2], col[1::2] # samples, controls
    @property
    def blocks(self):
        if self.blocknums is not None:
            return [self.block(i, 
                        smiles = j) for i, j in zip(self.blocknums, self.smiles)]
    def block(self, n, smiles = None):
        if smiles is None:
            if self.smiles is not None:
                smiles = self.smiles[n]
            else:
                smiles = None
        data, control = self.locblock(n)
        print(data)
        return Block(data = data,
                    control = control, 
                    concs = self.concs,
                    smiles = smiles)

class UV384m2(UV384):
    # hand assay v2 (controls in seperate plate) plate obj
    # compounds in all wells, each col: odd = block.1 even = block.2
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None):
        super().__init__(path, name, concs, smiles, names, parser )
    @property
    def blocknums(self):
        return range(1, (len(self.df) // 8) + 1)

    def locblock(self, num):
        col = (num + 1) // 2 
        assert col < 25
        if num % 2 == 1:
            rows = ascii_uppercase[:16:2]
        else:
            rows = ascii_uppercase[1:16:2]
        wells = [f'{i}{col}' for i in rows]
        return self.df.loc[wells,:]
    

class UV384m3(UV384m2):
    # echo asssay
    # no controls, they're in a seperate plate
    # each column has two blocks 
    # one in the first 8 rows and the other in the next 8
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None):
        super().__init__(path, name, concs, smiles, names, parser)
    @property
    def blocknums(self):
        return list(range(1, len(self.df) // 8))
    def locblock(self, num):
        # first 8 or last 8
        col = (num + 1) // 2 
        if num % 2 == 1:
            rows = ascii_uppercase[:8]
        else:
            rows = ascii_uppercase[8:16]
        wells = [f'{i}{col}' for i in rows]
        return self.df.loc[wells,:]

class UV384m4(UV384m3):
    # sideways plocks!
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None, control = None):
        super().__init__(path, name, concs, smiles, names, parser)
        self.control = control
    def locblock(self, num):
        # first 8 or last 8
        row = ascii_uppercase[(num - 1)  // 3 ]
        if num % 3 == 1:
            cols = range(1,9)
        elif num % 3 == 2:
            cols = range(9,17)
        elif num % 3 == 0:
            cols = range(17,25)
        wells = [f'{row}{i}' for i in cols]
        return self.df.loc[wells,:]

def smooth(data):
    normData = gaussian_filter1d(data, 1, axis = 1) # gaussian smoothing
    return pd.DataFrame(normData, index = data.index, columns = data.columns)
