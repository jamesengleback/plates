import os
import pandas as pd
import matplotlib.pyplot as plt
import  multiprocessing
from tqdm import  tqdm
import plates

def report(plate, blockNum, concs, name = None,  smiles = None, save_path = None):
    # plot
    if smiles is None:
        fig, ax = plt.subplots(2,2)
    else:
        fig, ax = plt.subplots(3,2)

    norm = plates.transforms.normalize(*plate.block(blockNum), smooth = True)
    diff = plates.transforms.difference(norm)
    response = plates.transforms.response(diff)
    km, vmax, c, rsq = plates.anal.MichaelisMenten(concs, response)

    plates.subplotTraces(norm, ax[0][0], concs = concs)
    plates.subplotTraces(diff, ax[0][1], concs = concs)
    plates.subplotMichaelisMenten(concs, response, ax[1][0], km, vmax, c, rsq)
    title = {'km':round(km, 3), 'vmax':round(vmax, 3), 'rsq':rsq}
    if name is not None:
        title['name'] = name
    plates.subplotText(ax[1][1], title)
    if smiles is not None:
        plates.subplotSMILES(ax[2][0], 'c1ccccc1')
        ax[2][1].axis('off')
    plt.tight_layout()
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
    plt.close()

    df = pd.DataFrame([[name, smiles, km, vmax, c, rsq]], columns = ['name', 'smiles','km','vmax','c','rsq'])
    return  df 

def reportPlate(plate, save_dir, smiles = None, names = None, concs = None):
    if smiles is None:
        smiles = plate.smiles # can be None
    if names is None:
        names = plate.names 
    if concs is None:
        concs = plate.concs

    assert len(smiles) == len(plate.blocks)
    assert len(names) == len(plate.blocks)
    assert len(concs) == 8 # might change

    os.makedirs(save_dir, exist_ok=True)

    dfs = []
    for i in tqdm(plate.blocks):
        dfs.append(report(plate = plate, 
            blockNum = i, 
            concs = concs, 
            smiles = 'c1ccccc1', 
            save_path = os.path.join(save_dir, f'{i}.png')))
    '''
    ## todo: parralel
    # issue - cannot pickle function
    helper = lambda num : report(plate = plate,
                                blockNum = num,
                                concs = concs,
                                smiles = 'CCCCCC=O',
                                save_path = os.path.join(save_dir, f'{num}.png'))

    with multiprocessing.Pool() as pool:
        results = pool.map(helper, plate.blocks)
        pool.join()
    print(list(results))
    '''
    df = pd.concat(dfs).reset_index(drop=True)
    # save
    df.to_csv(os.path.join(save_dir), 'report.csv')
    return df
