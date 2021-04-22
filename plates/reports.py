import os
import pandas as pd
import matplotlib.pyplot as plt
import  multiprocessing
from tqdm import  tqdm
import plates

def report(samples, controls = None, concs = None, name = None,  smiles = None, save_path = None):
    # plot
    if smiles is None:
        fig, ax = plt.subplots(2,2, figsize=(10,10))
    else:
        fig, ax = plt.subplots(3,2, figsize=(15,10))

    norm = plates.transforms.normalize(samples, controls, smooth = True)
    diff = plates.transforms.difference(norm)
    response = plates.transforms.response(diff)
    km, vmax, c, rsq = plates.anal.MichaelisMenten(concs, response)

    if concs is None:
        concs = range(len(samples))
    plates.subplotTraces(norm, ax[0][0], concs = concs)
    plates.subplotTraces(diff, ax[0][1], concs = concs)
    plates.subplotMichaelisMenten(concs, response, ax[1][0], km, vmax, c, rsq)
    title = {'km':round(km, 3), 'vmax':round(vmax, 3), 'rsq':rsq}
    if name is not None:
        title['name'] = name
    plates.subplotText(ax[1][1], title)
    if smiles is not None:
        try:
            plates.subplotSMILES(ax[2][0], smiles)
        except Exception as e:
            print(e)
        ax[2][1].axis('off')
    plt.tight_layout()
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
    plt.close()

    df = pd.DataFrame([[name, smiles, km, vmax, c, rsq]], columns = ['name', 'smiles','km','vmax','c','rsq'])
    return  df 
