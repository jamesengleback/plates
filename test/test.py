from unittest import TestCase
import plates

class testPlate(TestCase):
    def testInitBMGStandardCsv(self):
        p = plates.Plate('data/dm-d1,2.CSV')
        meta = p.meta
        df = p.df
        fmtdWells = p.formatWells(df.index)
        rowA = p.row('A')
        col1 = p.col(1)
    def testInitBMGascii(self):
        p = plates.Plate('data/3march-echo-wt-assay-ascii.csv',
                         parser='bmg-ascii')
        parser = p.parser
        path = p.path
        dataframe = p.dataframe
        parser = p.parser
        meta = p.meta
        df = p.df
        fmtdWells = p.formatWells(df.index)
        rowA = p.row('A')
        col1 = p.col(1)

# class Plate:
# class Block:
# class BlocksConatiner:
# class UV384(plates.Plate):
# class UV384m1(UV384):
# class UV384m2(UV384):
# class UV384m3(UV384m2):
# class UV384m4(UV384m3):
# def r_squared(yi,yj):
# def MichaelisMenten(x,y):
# def __init__(self, csv, name = None, parser = None):
# def meta(self):
# def df(self):
# def formatWells(self, wells):
# def read_standard_BMG_csv(self, path):
# def read_BMG_ascii(self, path):
# def row(self, letter):
# def col(self, num):
# def rolling_var(x, ksize = 50):
# def anomaly_detection(traces):
# def report(block, name = None, save_path = None):
# def subplotTraces(data, ax, concs = None, anomalies = None):
# def subplotMichaelisMenten(x, y, ax, km, vmax, rsq):
# def subplotText(ax, dictionary):
# def subplotSMILES(ax, smiles, name = None):
# def __init__(self, data, control = None, test_concs = None, ctrl_concs = None, smiles = None, smooth = False, norm410 = False):
# def df(self):
# def norm(self):
# def diff(self):
# def response(self):
# def mm(self):
# def __repr__(self):
# def report(self):
# def __init__(self, blocks):
# def __len__(self):
# def __getitem__(self, i):
# def __repr__(self):
# def plate(self):
# def __init__(self, path, name = None, parser = None):
# def blocknums(self):
# def block(self, n):
# def blocks(self):
# def __init__(self, path, name = None, parser = None):
# def blocknums(self):
# def block(self, num):
# def __init__(self, path, name = None, parser = None):
# def blocknums(self):
# def block(self, num):
# def __init__(self, path, name = None, parser = None):
# def blocknums(self):
# def block(self, num):
# def __init__(self, path, name = None, parser = None):
# def block(self, num):
# def smooth(data):
