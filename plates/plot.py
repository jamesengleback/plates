import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.optimize import curve_fit
from rdkit import Chem
from rdkit.Chem import Draw
from tqdm import tqdm
import plates

def subplotTraces(data, ax, concs = None, anomalies = None):
    # dataframe format: Plate.df
    # concs:  array-like
    if len(data) > 0:
        if concs is None:
            colors = plt.cm.inferno(np.linspace(0,1,len(data)))
        else:
            concs = concs[data.index]
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

def report(plate, num, concs = None, name = None, smiles = None, save_path = None):
    if concs is None:
        concs= np.array(range(8))

    samples, controls = plate.block(num)
    anomalies_idx = anomaly_detection(samples) # filter
    anomalies = samples.loc[anomalies_idx,:]
    samples = samples.drop(anomalies.index)
    controls = controls.loc[samples.index, :]

    if len(samples) > 0:
        print(samples)
        traces= plate.normalizeBlock(samples, controls, smooth = True)
        print(traces)
        difference = plate.differenceBlock(traces)
        response = plate.responseBlock(difference)
        km, vmax, c, rsq = plates.anal.MichaelisMenten(concs, response)

        fig = plt.figure(figsize = (10,10))
        if smiles is None:
            grid = plt.GridSpec(2,2)
        else:
            grid = plt.GridSpec(3,2)
        ax1 = fig.add_subplot(grid[0,0])
        subplotTraces(traces, ax1, concs, anomalies = anomalies)
        ax2 = fig.add_subplot(grid[1,0])
        subplotTraces(difference, ax2, concs)
        ax3 = fig.add_subplot(grid[0,1])
        subplotMichaelisMenten(concs, response, ax3)
        ax4 = fig.add_subplot(grid[1,1])
        labels = {'km':round(km, 3) ,'vmax': round(vmax, 3), 'r squared':rsq}
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
    else:
        pass

class plot:
    # throw away??
    def report(save_path = None):

        output = pd.DataFrame([], columns = ['name', 'smiles', 'kd', 'vmax', 'r squared'])
        os.makedirs(save_path, exist_ok=True)
        # todo: parralel
        def helper(block):
            pass 

        for i, j, k  in tqdm(zip(self.blocks, self.smiles, self.names),
                total = len(self.blocks)):
            report(plate = self, 
                num = i,
                concs = self.concs,
                smiles = j,
                name=k, 
                save_path= os.path.join(save_path, f'{k}.png'))
            kd, vmax, rsq = plates.anal.MichaelisMenten(self.concs, self.response(i))

            output = output.append(pd.DataFrame([[k, j, kd, vmax, rsq]], 
                columns = ['name', 'smiles', 'kd', 'vmax', 'r squared']))
        output.reset_index(inplace=True, drop=True)
        output.to_csv(os.path.join(save_path, f'{self.name}-report.csv'))


