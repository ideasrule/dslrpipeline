import pyfits
import sys
import numpy as np
import os

ch = int(sys.argv[1])
filelist = sys.argv[2:]
subpixfile = "subpix" + str(ch) + ".fits"
offsets = 0.25*np.array([[1,-1], [1,1], [-1,-1], [-1,1]])

prevfile = None

for i in range(len(filelist)):
    photfile = filelist[i]
    fitsfile = filelist[i][:-5] + ".fits"
    tmpfile = filelist[i][:-5] + ".tmp"
    psffile = filelist[i][:-5] + ".psf"

    call_str = "awk '{print $1,$2+" + str(offsets[ch][0]) + ",$3+" + str(offsets[ch][1]) + "}' " + photfile + " > " + tmpfile
    print call_str
    os.system(call_str)

    call_str = "cat " + tmpfile + " | ~hatuser/SVN/HATpipe/source/subpixel_sensitivity/FitPSF " + fitsfile + " --input-columns id,x,y --output-columns id,x,y,S,D,K,A,flux,flux_err,bg,bg_err,nbgpix,flag,sn,chi2,npix --subpix=" + subpixfile + " --gain 2.1 --psf-model=sdk --fit-order 5 --maxS 10"
    if prevfile != None: call_str += " --initial-guess=" + prevfile + " > " + psffile
    else: call_str += " > " + psffile
    print call_str
    os.system(call_str)

    prevfile = psffile
