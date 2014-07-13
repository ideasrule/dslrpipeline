#!/usr/bin/python
import pyfits
import numpy as np
from optparse import OptionParser

def extract(filename, output, plane):
    hdulist = pyfits.open(filename)
    if len(hdulist) == 4:
        data = np.array(hdulist[plane].data, dtype='uint16')
        newhdu = pyfits.PrimaryHDU(data, hdulist[plane].header)
        newhdu.writeto(output, clobber=True)
        return
    hdu=hdulist[0]
    hdu.data = np.array(hdu.data[plane], dtype='uint16')
    hdu.writeto(output, clobber=True)

def coadd(filename, output):
    hdulist = pyfits.open(filename)
    if len(hdulist)==1:
        data = np.array(hdulist[0].data, dtype='uint16')
    else:
        data = np.array((hdulist[0].data, hdulist[1].data, hdulist[2].data,
            hdulist[3].data))

    result = data[0] + data[1] + data[2] + data[3]
    result = np.array(result, dtype='uint16')
    newhdu = pyfits.PrimaryHDU(result, hdulist[0].header)
    newhdu.header["COLORCH"] = "added"
    newhdu.header["COLOR"] = "added"
    newhdu.writeto(output, clobber=True)

#not used for now
'''
def expand(filename, output):
    #not used
    hdulist = pyfits.open(filename)
    hdu = hdulist[0]
    data = np.array(hdu.data)
    big_arr = np.zeros((data.shape[1]*2, data.shape[2]*2), dtype='uint16')
    for i in range(data.shape[1]):
        for j in range(data.shape[2]):
            big_arr[2*i][2*j] = data[2][i][j]
            big_arr[2*i][2*j+1] = data[0][i][j]
            big_arr[2*i+1][2*j] = data[3][i][j]
            big_arr[2*i+1][2*j+1] = data[1][i][j]
    hdu.data = big_arr
    hdu.writeto(output, clobber=True)
'''

parser = OptionParser()
parser.add_option('-o', dest='output', help='output filename',
                  default='output.fits')
parser.add_option('-c', dest='channel',  type='int',
                  help='Channel to extract (0-3); omit to coadd')
(options, args) = parser.parse_args()

if len(args) != 1:
    print "Must give input filename!"
    exit(-1)
if options.channel == None:
    coadd(args[0], options.output)
else:
    extract(args[0], options.output, options.channel)
