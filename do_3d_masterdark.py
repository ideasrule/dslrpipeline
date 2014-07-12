import sys
import random
import os
import numpy as np
import pyfits

biasdata = pyfits.open(sys.argv[1])[0].data
filelist = sys.argv[2:]

tmpfolder = "tmp/dark" + str(random.randint(0,10000000)) + "/"
print tmpfolder
os.mkdir(tmpfolder)

for filename in filelist:
    #if it doesn't meet requirements, ignore
    hdu=pyfits.open(filename)[0]
    if hdu.header["imagetyp"] != "dark" or hdu.header["aborted"]==1 or \
            hdu.header["exptime"] < 179:
        continue
    print filename
    #decompose into 4 files, pretend they're separate
    for ch in range(4):
        basename = os.path.basename(filename)[:-6]
        newname = tmpfolder + basename + str(ch) + ".fits"
        call_str = "mcolor.py -c " + str(ch) + " -o " + newname + \
            " " + filename
        os.system(call_str)

#for each channel, invoke do_masterdark.py

for ch in range(4):
    biasname = tmpfolder + "masterbias" + str(ch) + ".fits"
    newhdu = pyfits.PrimaryHDU(biasdata[ch])
    newhdu.writeto(biasname, clobber=True)
    darkname = "masterdark_3d" + str(ch) + ".fits"

    call_str = "python ~hatuser/HATpipebin/include/HATpipepy/Actions/do_masterdark.py --multi-night --verbosity DEBUG --master-bias=" + biasname + " HATnet " + tmpfolder + "*_" + str(ch) + ".fits" + " -o " + darkname
    print call_str
    os.system(call_str)

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

