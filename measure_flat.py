import pyfits
import sys
import numpy as np

hdu=pyfits.open(sys.argv[1])[0]

print np.mean(hdu.data[0][0:1730,0:29]), np.mean(hdu.data[0][0:1730,73:-1])
