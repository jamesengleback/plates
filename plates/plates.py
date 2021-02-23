import re 
import numpy as np
import pandas as pd 

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

def rolling_var(x, ksize = 50):
    return np.array([x[i:i+ksize].std() for i in range(0, len(x) - ksize, round(len(x) / ksize))])

def anomaly_detection(traces):
    # pd.DataFrame, return idx of anomalies
    variances = [rolling_var(traces.loc[i,:], ksize=20) for i in traces.index]
    filter_idx = np.array([i.mean() > 0.01 for i in variances])
    return filter_idx

