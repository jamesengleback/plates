import sys 
import os
import re 
from string import ascii_uppercase
import numpy as np
import pandas as pd 
from scipy.ndimage import gaussian_filter1d
from tqdm import  tqdm

import plates



class UV384(plates.Plate):
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None):
        super().__init__(path, name, parser)
        if concs is None:
            self.concs = np.array(range(8))
        else:
            self.concs = concs 
        if smiles is None:
            self.smiles = [None] * 24
        else:
            self.smiles = smiles
        if names is None:
            self.names = range(24)
        else:
            self.names = names

class UV384m1(UV384):
    # hand assay v1 (controls in each plate) plate obj
    # each col is 1 compound and contains a control in every other well
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None):
        super().__init__(path, name, concs, smiles, names, parser )
    @property
    def blocks(self):
        return range(1, len(self.df) // 16)

    def block(self, num):
        col = self.col(num)
        return col[::2], col[1::2] # samples, controls

    def report(self, save_dir = None, concs = None, smiles = None, names = None):
        if concs is None:
            concs = np.array(range(1,9))
        if smiles is None:
            smiles = [None for i in range(len(self.blocks))]
        if names is None:
            names = [None for i in range(len(self.blocks))]
        if save_dir is None:
            save_dir = f"{os.path.basename(self.path)}-report"
        os.makedirs(save_dir, exist_ok=True)
        dfs = []
        for i, j in tqdm(enumerate(self.blocks), total = len(self.blocks)):
            samples, controls = self.block(j)
            x = plates.reports.report(samples = samples, 
                    controls = controls, 
                    concs = concs, 
                    save_path = os.path.join(save_dir, f'{i}.png'), 
                    smiles = smiles[i],
                    names = names)
            dfs.append(x)
        report = pd.concat(dfs, axis = 1, join='inner')
        report.to_csv(os.path.join(save_dir, 'metrics.csv'))
        return report

class UV384m2(UV384):
    # hand assay v2 (controls in seperate plate) plate obj
    # compounds in all wells, each col: odd = block.1 even = block.2
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None):
        super().__init__(path, name, concs, smiles, names, parser )
    @property
    def blocks(self):
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

    def report(self, controlPlate = None, save_dir = None, concs = None, smiles = None, names = None):
        if concs is None:
            concs = np.array(range(1,9))
        if smiles is None:
            smiles = [None for i in range(len(self.blocks))]
        if names is None:
            names = [None for i in range(len(self.blocks))]
        if save_dir is None:
            save_dir = f"{os.path.basename(self.path).split('.')[0]}-report"

        os.makedirs(save_dir, exist_ok=True)
        dfs = []
        for i, (j, k) in tqdm(enumerate(zip(self.blocks, smiles)), total = len(self.blocks)):
            samples = self.block(j)
            if controlPlate is None:
                controls = None
            else:
                controls = controlPlate.block(j)
            x = plates.reports.report(samples = samples, 
                    controls = controls,
                    concs = concs, 
                    save_path = os.path.join(save_dir, f'{j}.png'), 
                    smiles = k)
            dfs.append(x)
        report = pd.concat(dfs, axis = 1, join='inner')
        report.to_csv(os.path.join(save_dir, 'metrics.csv'))
        return report
    

class UV384m3(UV384m2):
    # echo asssay
    # no controls, they're in a seperate plate
    # each column has two blocks 
    # one in the first 8 rows and the other in the next 8
    def __init__(self, path, name = None, concs = None, smiles = None, names = None, parser = None):
        super().__init__(path, name, concs, smiles, names, parser)
    @property
    def blocks(self):
        return list(range(1, len(self.df) // 8))
    def block(self, num):
        # first 8 or last 8
        col = (num + 1) // 2 
        if num % 2 == 1:
            rows = ascii_uppercase[8:16]
        else:
            rows = ascii_uppercase[:8]
        wells = [f'{i}{col}' for i in rows]
        return self.df.loc[wells,:]

class transforms:
    def smooth(data):
        normData = gaussian_filter1d(data, 1, axis = 1) # gaussian smoothing
        return pd.DataFrame(normData, index = data.index, columns = data.columns)

    def normalize(samples, controls = None, smooth = False):
        # samples and controls
        if smooth:
            samples = transforms.smooth(samples.astype(float))
            if controls is not None:
                controls = transforms.smooth(controls.astype(float))
        
        # normalize at 800 nm
        samples = samples.subtract(samples[800], axis = 0).reset_index(drop=True)

        # subtract controls
        if controls is not None:
            controls = controls.subtract(controls[800], axis = 0).reset_index(drop=True)
            samples = samples.subtract(controls, axis = 0)

        # scale at 405 nm inflection point
        inflection = samples.loc[:,405]
        scaling = 1 / inflection
        return samples.multiply(scaling, axis = 0)

    def difference(data):
        # normalized block
        return data.subtract(data.loc[data.index[0],:], axis = 1)

    def response(data):
        # difference block
        return abs(data.loc[:,420]) + abs(data.loc[:, 390])

