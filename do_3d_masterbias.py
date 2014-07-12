import sys
import random
import os
import numpy as np
import pyfits

filelist = sys.argv[1:]

tmpfolder = "tmp/bias" + str(random.randint(0,10000000))
os.mkdir(tmpfolder)

for filename in filelist:
    #decompose into 4 files, pretend they're separate
    for ch in range(4):
        basename = os.path.basename(filename)[:-6]
        newname = tmpfolder + "/" + basename + str(ch) + ".fits"
        call_str = "mcolor.py -c " + str(ch) + " -o " + newname + \
            " " + filename
        os.system(call_str)

#for each channel, invoke do_masterbias.py
biasdata = []
for ch in range(4):
    biasname = tmpfolder + "/master_bias" + str(ch) + ".fits"
    call_str = "python ~hatuser/HATpipebin/include/HATpipepy/Actions/do_masterbias.py HATnet " + tmpfolder + "/*_" + str(ch) + ".fits" + " -o " + biasname
    print call_str
    os.system(call_str)
    chdata = pyfits.open(biasname)[0].data
    biasdata.append(chdata)

newhdu = pyfits.PrimaryHDU(np.array(biasdata))
newhdu.writeto("masterbias_3d.fits")
