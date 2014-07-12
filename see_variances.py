import pyfits
import sys
import numpy as np

def reject_outliers(data, m=5.):
    data = data.flatten()
    return data[abs(data - np.mean(data)) < m*np.std(data)]

hdulist=pyfits.open(sys.argv[1], mode="update")
hdu = hdulist[0]

vars = []
for ch in range(4):
    data = reject_outliers(hdu.data[ch][0:1728,73:])
    vars.append(np.var(data))

print vars
