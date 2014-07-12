#!/usr/bin/python
import numpy as np
import sys
import pyfits

hdu = pyfits.open(sys.argv[1])[0]
data = hdu.data
#print np.mean(data[0][0:1740,0:29]), np.mean(data[1][0:1740,0:29]), \
 #   np.mean(data[2][0:1740,0:29]), np.mean(data[3][0:1740,0:29])
print np.mean(data[0][0:1740,31:71]), np.mean(data[1][0:1740,31:71]), \
    np.mean(data[2][0:1740,31:71]), np.mean(data[3][0:1740,31:71])
