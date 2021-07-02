import sys 
import os
import re 
from string import ascii_uppercase
import warnings
import numpy as np
import pandas as pd 
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from tqdm import  tqdm

import plates # ???

class Block:
    def __init__(self, data, control = None, test_concs = None, ctrl_concs = None, smiles = None, smooth = False, norm410 = False):

        self.data = data
        self.control = control
        self.test_concs = test_concs
        self.ctrl_concs = ctrl_concs
        self.smiles = smiles
        self.smooth = smooth
        self.norm410 = norm410
        self.wells = self.data.index.tolist()
    @property
    def df(self):
        if self.control is None:
            return self.data
        else:
            return self.data, self.control
    @property
    def norm(self):
        df = self.data
        df = df.sub(df.iloc[:,-1], axis=0) # anchor at 800 nm
        if self.smooth:
            df = smooth(df)

        if self.control is not None:
            # build mask
            ctrl = self.control
            ctrl = ctrl.sub(ctrl.iloc[:,-1], axis=0)
            if self.smooth:
                ctrl = smooth(ctrl)

            #missing = [j != i for i,(j,k) in enumerate(zip(self.ctrl_concs, self.test_concs))]
            #if sum(missing) > 0:
            #    interp = interp1d(self.ctrl_concs, ctrl.T, kind = 'linear', fill_value='extrapolate')
            #    if not interp.bounds_error:
            #        ctrl_interp = pd.concat([pd.Series(interp(i), name = k, index = ctrl.columns)\
            #                if j else ctrl.loc[k,:] for i,j,k in zip(self.test_concs, missing, ctrl.index)],
            #                axis=1).T
            #        if sum(ctrl_interp.isna()) == 0:
            #            ctrl = ctrl_interp

            df = df.sub(ctrl)
            df = df.dropna(axis=1)
        if self.norm410:
            # note - in this case, vmax is only comparable within runs that also norm410
            df = df.div(df.loc[:,410], axis=0) * 0.3 # so it fits in the axes
        return df
    @property
    def diff(self):
        norm = self.norm
        return norm.sub(norm.iloc[0,:])
    @property
    def response(self):
        diff = self.diff
        response = diff.loc[:,390].abs() + diff.loc[:,420].abs()
        if self.test_concs is not None:
            assert len(self.test_concs) == len(response)
            #if len(set(self.test_concs)) > 1: # in case they're all the same, then they can still be plotted
            response.index = self.test_concs
        return response
    @property
    def mm(self):
        assert self.test_concs is not None, 'Need some concentrations bud!'
        response = self.response
        # index makes problems in mm
        return plates.analysis.MichaelisMenten(self.test_concs.reset_index(drop=True), 
                    response.reset_index(drop=True))
    def __repr__(self):
        return f'\n{"-"*10}\nblock num: {self.block_num}\nplate: {self.plate_name}\npath:{self.plate_path}\nwells: {self.wells}\nsmiles: {self.smiles}\n{"-"*10}'

    def report(self):
        pass

class BlocksConatiner:
    def __init__(self, blocks):
        self.blocks = blocks
    def __len__(self):
        return len(self.blocks)
    def __getitem__(self, i):
        return self.blocks[i]
    def __repr__(self):
        return f'number of blocks: {len(self)}\nplate: {self.plate}'
    @property
    def plate(self):
        return set([i.plate_name for i in self.blocks])


class UV384(plates.Plate):
    def __init__(self, path, name = None, parser = None):
        super().__init__(path, name, parser)
    @property
    def blocknums(self):
        pass
    def block(self, n):
        # just return df slice
        pass
    @property
    def blocks(self):
        return [self.block(i) for i in self.blocknums]

class UV384m1(UV384):
    # hand assay v1 (controls in each plate) plate obj
    # each col is 1 compound and contains a control in every other well
    def __init__(self, path, name = None, parser = None):
        super().__init__(path, name, parser )
    @property
    def blocknums(self):
        return range(1, len(self.df) // 16)
    def block(self, num):
        col = self.col(num)
        return col[::2], col[1::2] # samples, controls
               
               
               

class UV384m2(UV384):
    # hand assay v2 (controls in seperate plate) plate obj
    # compounds in all wells, each col: odd = block.1 even = block.2
    def __init__(self, path, name = None, parser = None):
        super().__init__(path, name, parser )
    @property
    def blocknums(self):
        return range(1, (len(self.df) // 8) + 1)

    def block(self, num):
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
    def __init__(self, path, name = None, parser = None):
        super().__init__(path, name, parser)
    @property
    def blocknums(self):
        return list(range(1, len(self.df) // 8))
    def block(self, num):
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
    def __init__(self, path, name = None, parser = None):
        super().__init__(path, name, parser)
    def block(self, num):
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
    normData = gaussian_filter1d(data, 3, axis = 1) # gaussian smoothing
    return pd.DataFrame(normData, index = data.index, columns = data.columns)
