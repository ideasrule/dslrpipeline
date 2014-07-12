import pyfits
import numpy as np
import sys
from scipy import polyfit, polyval
from scipy.stats import linregress

def reject_outliers_iter(data, m=5.):
    data = data.flatten()
    return data[abs(data - np.mean(data)) < m*np.std(data)]

def reject_outliers(data):
    prev_len = 0
    data = data.flatten()
    for i in range(5):
        data = reject_outliers_iter(data)
        if len(data) == prev_len: return data
        prev_len = len(data)


filelist = sys.argv[1:]

for ch in range(4):
    sum_data = []
    xs = []
    ys = []
    for i in range(len(filelist)):
        data = pyfits.open(filelist[i])[0].data[ch][0:1728,73:]
#        data = reject_outliers(data)
        if i==0: sum_data = data
        else: sum_data += data
        mean_data = sum_data/(i+1)
        xs.append(1./(i+1))
        ys.append(np.var(mean_data))

    slope, intercept, r, pvalue, see = linregress(xs,ys)
    print slope, intercept, r
