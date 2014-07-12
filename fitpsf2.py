import pyfits
import sys
import numpy as np
import os

filelist = sys.argv[1:]
offsets = 0.25*np.array([[1,-1], [1,1], [-1,-1], [-1,1]])

for filename in filelist:
    try:
        ch = int(filename[-6])
    except ValueError:
        continue
    subpixfile = "subpix" + str(ch) + ".fits"
    photfile = filename
    fitsfile = filename[:-5] + ".fits"
    tmpfile = filename[:-5] + ".tmp"
    psffile = filename[:-5] + ".psf"
    if not os.path.isfile(psffile): continue

    call_str = "awk '{print $1,$2+" + str(offsets[ch][0]) + ",$3+" + str(offsets[ch][1]) + "}' " + photfile + " > " + tmpfile
    print call_str
    os.system(call_str)

    call_str = "cat " + tmpfile + " | ~hatuser/SVN/HATpipe/source/subpixel_sensitivity/FitPSF " + fitsfile + " --input-columns id,x,y --output-columns id,x,y,S,D,K,A,flux,flux_err,bg,bg_err,nbgpix,flag,sn,chi2 --subpix=" + subpixfile + " --gain 2.1 --fit-order 5 --maxS 10 --psf-model=sdk --initial-guess=" + psffile + " > " + filename[:-5] + ".psf2"
    print call_str
    os.system(call_str)

