import sys 
import os
import re 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from rdkit import Chem
from rdkit.Chem import Draw
from tqdm import tqdm

class Plate:
    # generic
    def __init__(self, csv, name = None):
        self.path = csv 
        if name is None:
            self.name = csv.split('.')[0]
        else:
            self.name = name

    @property
    def meta(self):
        # dict
        with open(self.path,'r') as f:
             m = ''.join([f.readline() for i in range(6)])
        metadata = {'date': re.findall('Date: (\d+/\d+/\d+)', m)[0],
                   'time': re.findall('Time: (\d+:\d+:\d+)', m)[0],
                   'protocol':re.findall('Test name: (.+?),',m)[0],
                   'plateID':re.findall('ID1: (.+)', m)[0]}
        return metadata

    @property
    def df(self):
        # clean
        df = pd.read_csv(self.path, skiprows=7)
        df.index = df.iloc[:,0] # wells
        df.drop(['Unnamed: 0','Wavelength [nm]'],axis=1, inplace=True)
        df.dropna(axis=1, inplace=True)
        df.columns = df.columns.astype(int) # wavelengths
        return df.astype(float)

    def row(self, letter):
        idx = pd.Series(self.df.index).str.contains(letter.upper())
        return self.df.loc[idx.to_list(), :]

    def col(self, num):
        idx = [int(re.findall('[A-Z]([0-9]+)', i)[0]) == int(num)  for i in self.df.index]
        return self.df.loc[idx, :]


class AssayPlate(Plate):
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

    def normalizeBlock(self, num, smooth = False):
        # block number
        samples, controls = [i.subtract(i[800], axis = 0).reset_index(drop=True) for i in self.block(num)] # normalize at 800 nm
        samples = samples.subtract(controls, axis = 0)
        # scale at 405 nm inflection point
        inflection = samples.loc[:,405]
        scaling = 1 / inflection
        normalizedSamples = samples.multiply(scaling, axis = 0)
        if smooth:
            normalizedSamples = gaussian_filter1d(normalizedSamples, 1, axis = 1) # gaussian smoothing
            normalizedSamples = pd.DataFrame(normalizedSamples, index = samples.index, columns = samples.columns)
        return normalizedSamples

    def differenceBlock(self, num):
        # block num
        data = self.normalizeBlock(num)
        return data.subtract(data.loc[0,:], axis = 1)

    def responseBlock(self, num):
        # block num
        data = self.differenceBlock(num)
        return data.loc[:,420] - data.loc[:, 390]

    def report(self, save_path = None):
        if save_path is None:
            if self.name is None:
                save_path = self.csv
            else:
                save_path = self.name
        output = pd.DataFrame([], columns = ['name', 'smiles', 'kd', 'vmax'])
        os.makedirs(save_path, exist_ok=True)
        for i, j, k  in tqdm(zip(self.blocks, self.smiles, self.names),
                total = len(self.blocks)):
            report(plate = self, 
                num = i,
                concs = self.concs,
                smiles = j,
                name=k, 
                save_path= os.path.join(save_path, f'{k}.png'))
            kd, vmax = MichaelisMenten(self.concs, self.responseBlock(i))

            output = output.append(pd.DataFrame([[k, j, kd, vmax]], 
                columns = ['name', 'smiles', 'kd', 'vmax']))
        output.reset_index(inplace=True, drop=True)
        output.to_csv(os.path.join(save_path, f'{self.name}-report.csv'))


def subplotTraces(data, ax, concs = None):
    # dataframe format: Plate.df
    # concs - array-like
    if concs is None:
        colors = plt.cm.inferno(np.linspace(0,1,len(data)))
    else:
        assert len(concs) == len(data)
        minmaxscale = lambda x : (x - min(x)) / max(x)
        colors = plt.cm.inferno(minmaxscale(concs))
    for i, j in zip(data.index, colors):
        ax.plot(data.loc[i, 300:], c = j, lw = 1)
    if concs is None:
        ax.legend(range(len(data)), title = 'Trace Number', loc = 'right')
    else:
        ax.legend(labels = list(concs), title = '[Ligand] µM', loc = 'right')
    ax.set_xlabel('Wavelength nm')
    ax.set_ylabel('Absorbance')
    ax.set_ylim(-0.5,2)
    ax.set_xlim(300, 800)
    ax.set_xticks(np.linspace(300, 800, 11))
    ax.vlines(390,-0.5,2, linestyle = '--', lw = 0.5, color = '0.5')
    ax.vlines(420,-0.5,2, linestyle = '--', lw = 0.5, color = '0.5')

def subplotMichaelisMenten(x, y, ax):
    y = y.replace(np.inf, 0) # error handling
    km, vmax = MichaelisMenten(x,y)
    ax.scatter(x,y)
    mm = lambda x, km, vmax : (x * vmax) / (km + x)
    xx = np.linspace(min(x), max(x), 100)
    ax.plot(xx, mm(xx, km, vmax))
    ax.set_ylim(min(y) - 0.1, max(y) * 1.1)
    ax.set_xlabel('[Ligand] µM')
    ax.set_ylabel('response')

def subplotText(ax, dictionary):
    s = ''.join([f'{i} = {round(j, 3)}\n' if type(j) == float else f'{i} = {j}\n' for i, j in zip(dictionary.keys(), dictionary.values())])
    ax.text(0.5,0.5, s, ha='center')
    ax.axis('off')


def MichaelisMenten(x,y):
    y = y.replace(np.inf, 0) # error handling
    mm = lambda x, km, vmax : (x * vmax) / (km + x)
    (km, vmax), covariance = curve_fit(mm, x, y, bounds=((0,0),(np.inf, np.inf)))
    return km, vmax

def report(plate, num, concs = None, name = None, smiles = None, save_path = None):
    if concs is None:
        concs= np.array(range(8))
    fig = plt.figure(figsize = (10,10))
    if smiles is None:
        grid = plt.GridSpec(2,2)
    else:
        grid = plt.GridSpec(3,2)

    ax1 = fig.add_subplot(grid[0,0])
    subplotTraces(plate.normalizeBlock(num, smooth = True), ax1, concs)

    ax2 = fig.add_subplot(grid[1,0])
    subplotTraces(plate.differenceBlock(num), ax2, concs)

    ax3 = fig.add_subplot(grid[0,1])
    subplotMichaelisMenten(concs, plate.responseBlock(num), ax3)

    ax4 = fig.add_subplot(grid[1,1])
    km, vmax = MichaelisMenten(concs, plate.responseBlock(num))
    labels = {'km':round(km, 3) ,'vmax': round(vmax, 3)}
    if name is not None:
        labels['Name'] = name
    subplotText(ax4, labels)
    
    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        ax5 = fig.add_subplot(grid[2,0])
        im = Draw.MolToImage(mol)
        ax5.imshow(im)
        ax5.axis('off')
        
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()
