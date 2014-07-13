import sys
import random
import os
import numpy as np
import pyfits
import utils

filelist = sys.argv[1:]

tmpfolder = "tmp/dark" + str(random.randint(0,10000000)) + "/"
os.mkdir(tmpfolder)

newfilelist = []

for filename in filelist:
    #if it doesn't meet requirements, ignore
    hdu=pyfits.open(filename)[0]
    if hdu.header["imagetyp"] != "dark" or hdu.header["aborted"]==1 or \
            hdu.header["exptime"] < 179:
        continue
    newfilelist.append(filename)

utils.split_by_channel(newfilelist, tmpfolder)

#for each channel, invoke do_masterdark.py

command_queue = []
for ch in range(4):
    biasname = "masterbias" + str(ch) + ".fits"
    darkname = "masterdark" + str(ch) + ".fits"

    call_str = "python ~hatuser/HATpipebin/include/HATpipepy/Actions/do_masterdark.py --multi-night --verbosity DEBUG --master-bias=" + biasname + " HATnet " + tmpfolder + "*_" + str(ch) + ".fits" + " -o " + darkname
    command_queue.append(call_str)

utils.run_all(command_queue)

