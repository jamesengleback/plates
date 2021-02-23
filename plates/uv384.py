import sys 
import os
import re 
import numpy as np
import pandas as pd 
from scipy.ndimage import gaussian_filter1d

import plates



class UV384m1(plates.Plate):
    # hand assay plate obj
    def __init__(self, path, name = None, concs = None, smiles = None, names = None):
        super().__init__(path, name)
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
    
    @property
    def blocks(self):
        return list(set([int(re.findall('[A-Z]([0-9]+)', i)[0])  for i in self.df.index]))

    def block(self, num):
        col = self.col(num)
        return col[::2], col[1::2] # samples, controls

class transforms:
    def smooth(data):
        normData = gaussian_filter1d(data, 1, axis = 1) # gaussian smoothing
        return pd.DataFrame(normData, index = data.index, columns = data.columns)

    def normalize(data1, data2, smooth = False):
        # samples and controls
        if smooth:
            samples = transforms.smooth(data1)
            controls = transforms.smooth(data2)
        else:
            samples = data1
            controls = data2
        
        # normalize at 800 nm
        samples = samples.subtract(samples[800], axis = 0).reset_index(drop=True)
        controls = controls.subtract(controls[800], axis = 0).reset_index(drop=True)

        # subtract controls
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

