import sys
import pyfits
import numpy as np

hdu = pyfits.open(sys.argv[1])[0]
for i in range(4):
    hdu.data[i] -= np.median(hdu.data[i][0:1730,0:29])
hdu.writeto(sys.argv[2])
