import numpy as np
import pyfits
import sys

hdu=pyfits.open(sys.argv[1])[0]
temp = hdu.header["ccdtemp"] + 273.15
val = hdu.data[0][161][2134] - 2048

print 1/temp, np.log(val)
