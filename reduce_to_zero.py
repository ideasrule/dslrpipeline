import numpy as np
import pyfits
import sys

def subtract_top(twoDdata):
    means = []
    reduced_means = []
    twoDdata = np.transpose(twoDdata)
    #1744 to 1756 inclusive
    for row in twoDdata:
        offset = np.mean(row[1731:1756])
        means.append(np.mean(row[0:1730]))
#        offset = np.median(row[0:1744])
        row -= offset
        reduced_means.append(np.mean(row[0:1730]))
    twoDdata =  np.transpose(twoDdata)
    print np.std(means), np.std(reduced_means)

hdu=pyfits.open(sys.argv[1])[0]

subtract_top(hdu.data[0])
subtract_top(hdu.data[1])
subtract_top(hdu.data[2])
subtract_top(hdu.data[3])
hdu.writeto("master_bias2.fits")
