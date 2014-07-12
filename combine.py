import numpy as np
import pyfits
import os

i=1
while 1==1:
    if not os.path.isfile("masterdark_3d0_" + str(i) + ".fits"): break
    data=[]
    for ch in range(4):
        filename = "masterdark_3d" + str(ch) + "_" + str(i) + ".fits"
        hdu=pyfits.open(filename)[0]
        data.append(hdu.data)
    hdu.data = np.array(data)
    hdu.writeto("masterdark_3d_" + str(i) + ".fits", clobber=True)
    i += 1
