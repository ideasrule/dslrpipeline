import pyfits
import sys
import numpy as np

def reject_outliers(data, m=5.):
    data = data.flatten()
    return data[abs(data - np.mean(data)) < m*np.std(data)]

def get_var(data):
    prev_len = 0
    data = data.flatten()
    for i in range(5):
        data = reject_outliers(data)
        variance = np.var(data)
        if len(data) == prev_len: return variance
        prev_len = len(data)

hdulist=pyfits.open(sys.argv[1], mode="update")
hdu = hdulist[0]

vars = []
print sys.argv[1]

for ch in range(4):
    vars.append(get_var(hdu.data[ch][0:1728,32:68]))
#print vars
#exit(-1)
hdu.data = np.array(hdu.data,dtype='uint16')
hdu.header["CCDTEMP"] = np.median(vars)
hdu.header["BITPIX"]=16
hdulist.flush()
hdulist.close()
