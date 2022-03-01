from string import ascii_uppercase
import re 
import numpy as np
from einops import rearrange
import pandas as pd 

class Plate:
    # generic
    def __init__(self, 
                 path, 
                 fmt=384,
                 name=None,
                 parser='bmg-standard'):
        assert parser in {'bmg-standard', 'bmg-ascii'}, \
                "parser must be one of {'bmg-standard', 'bmg-ascii'}"
        assert fmt in {96, 384}
        self.parser = parser
        self.path = path 
        self.fmt = fmt 
        self.dataframe = None
        if name is None:
            self.name = path.split('.')[0]
        else:
            self.name = name
    def __len__(self):
        return len(self.df)
    def __iter__(self):
        for i in self.wells:
            yield self.df.loc[i,:]
    def __getitem__(self, idx):
        pass
    def __repr__(self):
        return ''
    @property
    def meta(self): # returns dict
        # todo - ascii parser
        with open(self.path,'r') as f:
             m = ''.join([f.readline() for i in range(6)])
        if self.parser == 'bmg-standard':
            return {'date': re.findall('Date: (\d+/\d+/\d+)', m)[0],
                    'time': re.findall('Time: (\d+:\d+:\d+)', m)[0],
                    'protocol':re.findall('Test name: (.+?),',m)[0],
                    'plateID':re.findall('ID1: (.+)', m)[0]}
        elif self.parser == 'bmg-ascii':
            # this format has no metadata
            return {}
        else:
            raise Exception('No parser: options are None (defau')
    @property
    def df(self):
        if self.dataframe is None: # caching
            if self.parser == 'bmg-standard':
                df = self.read_standard_BMG_csv(self.path)
                df.index = self.formatWells(df.index)
                get_letter = lambda x : re.findall('([A-Z]+)', x)[0]
                get_num = lambda x : int(re.findall('([0-9]+)', x)[0]) 
                key_fn = lambda x: get_num(x) + (ascii_uppercase.index(get_letter(x))/16)
                idx = df.index.to_list()
                idx.sort(key=key_fn)
                df = df.loc[idx,:]
                self.dataframe = df 
            elif self.parser == 'bmg-ascii':
                with open(self.path) as f:
                    __data = [i for i in f.read().split('\n\n') if  i != ''  ]
                _data = np.loadtxt(self.path, delimiter=',')
                # data.shape eg : (551, 24, 16) (wavelengths, row, col)
                data = rearrange(_data, '(w r) c -> w c r', w=len(__data))
                df = pd.DataFrame(rearrange(data, 'w c r -> (c r) w'))
                if self.fmt == 96:
                    df.index = [i for i, j in zip(\
            [f'{k}{l}' for k in ascii_uppercase[:9] for l in range(1,13)],
                            df.index)]
                elif self.fmt == 384:
                    df.index = [i for i, j in zip(\
            [f'{k}{l}' for k in ascii_uppercase[:17] for l in range(1,25)],
                            df.index)]
                return df

            else:
                raise Exception('Parser not recognized in plates.Plate.df')
        return self.dataframe
    @property
    def wells(self):
        return self.formatWells(self.df.index)
    def formatWells(self, wells): # A01 -> A1
        fn = lambda  x : f"{re.findall('[A-Z]+', x)[0]}{int(re.findall('[0-9]+',x)[0])}"
        return [fn(i) for i in wells]
    def read_standard_BMG_csv(self, path):
        df = pd.read_csv(path, skiprows=7)
        df.index = df.iloc[:,0] # wells
        df.drop(['Unnamed: 0','Wavelength [nm]'],axis=1, inplace=True)
        df.dropna(axis=1, inplace=True)
        df=df.replace('overflow', 3.5)
        df.columns = df.columns.astype(int) # wavelengths
        return df.astype(float)
    def read_BMG_ascii(self, path):
        # todo - robustness with plate sizes
        with open(path, 'r') as f:
            data = f.read()
        data = data.split('\n\n') # each wavelength in layout 16 * 24 
        data = [i.split('\n') for i in data] # rows in each wavelength block
        data = [np.array([i.split(',') for i in j]) for j in data]
        data = np.stack([i for i in data if i.shape != (1,1)], axis=0)
        data = np.swapaxes(data, 0,2) # shape (24, 16, 551)
        wells = [f'{i}{j}' for j in range(1,25) for i in ascii_uppercase[:16]]

        return pd.DataFrame([j for i in data for j in i],
                index = wells,
                columns = range(250,801)) # todo - check wavelengths
    def row(self, letter):
        idx = pd.Series(self.df.index).str.contains(letter.upper())
        return self.df.loc[idx.to_list(), :]

    def col(self, num):
        idx = [int(re.findall('[A-Z]([0-9]+)', i)[0]) == int(num)  for i in self.df.index]
        return self.df.loc[idx, :]

def rolling_var(x, ksize = 50):
    return np.array([x[i:i+ksize].std() for i in range(0, len(x) - ksize, round(len(x) / ksize))])

def anomaly_detection(traces):
    # pd.DataFrame, return idx of anomalies
    variances = [rolling_var(traces.loc[i,:], ksize=20) for i in traces.index]
    filter_idx = np.array([i.mean() > 0.01 for i in variances])
    return filter_idx

