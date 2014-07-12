import pyfits
import sys
import numpy as np

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

hdu=pyfits.open(sys.argv[1])[0]
vars=[]

for ch in range(4):
    data = hdu.data[ch]
    vars.append(np.var(reject_outliers(data[0:1728,1600:])))
print vars[0],vars[1],vars[2],vars[3]
#print np.var(hdu.data[
