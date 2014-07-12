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

bias=pyfits.open(sys.argv[1])[0].data
dark=pyfits.open(sys.argv[2])[0].data
target_hdu=pyfits.open(sys.argv[3])[0]
print target_hdu.header["ccdtemp"]

Tabs = target_hdu.header["CCDTEMP"] + 273.15
factor = np.exp(-12800*(1./Tabs - 1./293.15))
print factor

target = target_hdu.data
target = target - bias - dark
target = reject_outliers(target[0][0:1728,73:])

print np.var(target)

#target_hdu.data=target
#target_hdu.writeto("dark_subtracted.fits")


