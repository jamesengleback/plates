import numpy as np
import pandas as pd 
from scipy.optimize import curve_fit

def r_squared(yi,yj):
    residuals = yi - yj
    sum_sq_residual = sum(residuals ** 2)
    sum_sq_total = sum((yi - yi.mean()) ** 2) # check this!!!
    return 1 - (sum_sq_residual / sum_sq_total)

def MichaelisMenten(x,y):
    y = y.replace(np.inf, 0) # error handling
    mm = lambda x, km, vmax, c : ((x * vmax) / (km + x)) + c
    (km, vmax, c), covariance = curve_fit(mm, x, y, bounds=((0, 0,0),(np.inf, np.inf, np.inf)))
    
    yh = mm(x, km, vmax, c)
    rsq = r_squared(y, yh)
    return km, vmax, c, rsq

