import sys
import pyfits
import random
import numpy as np
import os

def choose_best_dark(filelist):
    darks=[]
    i=1
    while 1==1:
        filename="masterdark_3d_" + str(i) + ".fits"
        if os.path.isfile(filename): darks.append(filename)
        else: break
        i += 1

    imagetemps = []
    for filename in filelist:
        imagetemps.append(pyfits.open(filename)[0].header["ccdtemp"])
    median_temp = np.median(imagetemps)

    min_temp_diff = 1000
    min_dark = None
    for filename in darks:
        hdu=pyfits.open(filename)[0]
        darktemp = hdu.header["tempmed"]
        diff = abs(darktemp - median_temp)
        if diff < min_temp_diff:
            min_temp_diff = diff
            min_dark = filename
    return min_dark

masterbias = sys.argv[1]
output = sys.argv[2]
filelist = sys.argv[3:]
if len(filelist)==0: exit(-1)

masterdark = choose_best_dark(filelist)
biasdata = pyfits.open(masterbias)[0].data
darkdata = pyfits.open(masterdark)[0].data

tmpfolder = "tmp/flat" + str(random.randint(0,10000000)) + "/"
os.mkdir(tmpfolder)

for filename in filelist:
    #decompose into 4 files, pretend they're separate
    for ch in range(4):
        basename = os.path.basename(filename)[:-6]
        newname = tmpfolder + basename + str(ch) + ".fits"
        call_str = "mcolor.py -c " + str(ch) + " -o " + newname + \
            " " + filename
        os.system(call_str)

#for each channel, invoke do_masterflat.py
flatdata = []
for ch in range(4):
    biasname = tmpfolder + "master_bias" + str(ch) + ".fits"
    darkname = tmpfolder + "master_dark" + str(ch) + ".fits"
    flatname = tmpfolder + "master_flat" + str(ch) + ".fits"

    newhdu = pyfits.PrimaryHDU(biasdata[ch])
    newhdu.writeto(biasname, clobber=True)
    newhdu.data = darkdata[ch]
    newhdu.writeto(darkname, clobber=True)

    call_str = "python ~hatuser/HATpipebin/include/HATpipepy/Actions/do_masterflat.py --master-bias=" + biasname + " --master-dark=" + darkname + " " + \
        "HATnet " + tmpfolder + "*_" + str(ch) + ".fits" + " -o " + flatname

    print call_str
    os.system(call_str)
    chdata = pyfits.open(flatname)[0].data
    flatdata.append(chdata)

newhdu = pyfits.PrimaryHDU(np.array(flatdata))
newhdu.writeto(output, clobber=True)
