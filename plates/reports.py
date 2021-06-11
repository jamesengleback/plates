import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem 
from rdkit.Chem import Draw
from tqdm import  tqdm

def report(block, name = None, save_path = None):
    # plot
    if block.smiles is None:
        fig, ax = plt.subplots(2,2, figsize=(10,10))
    else:
        fig, ax = plt.subplots(3,2, figsize=(15,10))

    norm = block.norm
    diff = block.diff
    response = block.response
    mm = block.mm
    km, vmax, c, rsq, = mm['km'], mm['vmax'], mm['c'], mm['rsq']

    concs = block.concs
    if concs is None:
        concs = range(len(samples))
    subplotTraces(norm, ax[0][0], concs = concs)
    subplotTraces(diff, ax[0][1], concs = concs)
    subplotMichaelisMenten(concs, response, ax[1][0], km, vmax, c, rsq)
    title = {'km':round(km, 3), 'vmax':round(vmax, 3), 'rsq':rsq}
    if name is not None:
        title['name'] = name
    subplotText(ax[1][1], title)
    if block.smiles is not None:
        try:
            subplotSMILES(ax[2][0], block.smiles)
        except Exception as e:
            print(e)
        ax[2][1].axis('off')
    plt.tight_layout()
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
    plt.close()



def subplotTraces(data, ax, concs = None, anomalies = None):
    # dataframe format: Plate.df
    # concs:  array-like
    if len(data) > 0:
        if concs is None:
            colors = plt.cm.inferno(np.linspace(0,1,len(data)))
        else:
            #concs = concs[data.index]
            assert len(concs) == len(data)
            minmaxscale = lambda x : (x - min(x)) / max(x)
            colors = plt.cm.inferno(minmaxscale(concs))
        for i, j in zip(data.index, colors):
            ax.plot(data.loc[i, 300:], c = j, lw = 1)
    if anomalies is not None:
        for i in anomalies.index:
            ax.plot(anomalies.loc[i, 300:], c = '0.5', lw = 0.5)
    if concs is None:
        ax.legend(range(len(data)), title = 'Trace Number', loc = 'right')
    else:
        ax.legend([round(i,2) for i in concs], title = '[Ligand] µM', loc = 'right')
    ax.set_xlabel('Wavelength nm')
    ax.set_ylabel('Absorbance')
    ax.set_ylim(-0.5,2)
    ax.set_xlim(300, 800)
    ax.set_xticks(np.linspace(300, 800, 11))
    ax.vlines(390,-0.5,2, linestyle = '--', lw = 0.5, color = '0.5')
    ax.vlines(420,-0.5,2, linestyle = '--', lw = 0.5, color = '0.5')

def subplotMichaelisMenten(x, y, ax, km, vmax, c, rsq):
    y = y.replace(np.inf, 0) # error handling
    ax.scatter(x,y)
    mm = lambda x, km, vmax : ((x * vmax) / (km + x)) + c
    xx = np.linspace(min(x), max(x), 100)
    ax.plot(xx, mm(xx, km, vmax))
    ax.set_ylim(min(y) - 0.1, max(y) * 1.1)
    ax.set_xlabel('[Ligand] µM')
    ax.set_ylabel('response')

def subplotText(ax, dictionary):
    s = ''.join([f'{i} = {round(j, 3)}\n' if type(j) == float else f'{i} = {j}\n' for i, j in zip(dictionary.keys(), dictionary.values())])
    ax.text(0.5,0.5, s, ha='center')
    ax.axis('off')

def subplotSMILES(ax, smiles, name = None):
        mol = Chem.MolFromSmiles(smiles)
        im = Draw.MolToImage(mol)
        ax.imshow(im)
        ax.axis('off')

