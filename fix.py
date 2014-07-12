#!/usr/bin/python
import os
import sys
import pyfits
import numpy as np

cr2_file = sys.argv[1]
fits_file = cr2_file[0:-4] + ".fits"

if not os.path.isfile(cr2_file) or not os.path.isfile(fits_file):
    exit(-1)

os.system("rawtran -c all -o !temp.fits " + cr2_file)
data = pyfits.open("temp.fits")[0].data

hdulist = pyfits.open(fits_file, mode="update")
hdulist[0].data = np.array(data, dtype='uint16')
hdulist.flush()

