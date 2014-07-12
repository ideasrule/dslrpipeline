import pyfits
import sys
import numpy as np
import os

ch = int(sys.argv[1])
filelist = sys.argv[2:]
subpixfile = "subpix" + str(ch) + ".fits"
offsets = 0.25*np.array([[1,-1], [1,1], [-1,-1], [-1,1]])

for i in range(len(filelist)):
    photfile = filelist[i]
    fitsfile = filelist[i][:-5] + ".fits"
    tmpfile = filelist[i][:-5] + ".tmp"
    psffile = filelist[i][:-5] + ".psf"
    outfile = filelist[i][:-5] + ".subpixphot"

    call_str = "awk '{print $1,$2+" + str(offsets[ch][0]) + ",$3+" + str(offsets[ch][1]) + "}' " + photfile + " > " + tmpfile
    print call_str
    os.system(call_str)

    call_str = "cat " + tmpfile + " | ./SubPixPhot " + fitsfile + " --sdk-poly=" + psffile + " --input-columns id,x,y --output-columns id,x,y,S,D,K,bg,bg_err,flux,flux_err,flag --subpix=" + subpixfile + " --gain 2.1 -a 1.6 -a 1.8 -a 2.0 -a 2.2 -a 2.4 -a 2.6 -a 2.8 -a 3.0 -a 3.2 -a 3.4"

    print call_str
    os.system(call_str)

