import sys
import pyfits
import random
import numpy as np
import os
import utils

def choose_best_dark(filelist):
    darks=[]
    i=1

    while 1==1:
        filename="masterdark0_" + str(i) + ".fits"
        if not os.path.isfile(filename): break
        i += 1

    imagetemps = []
    for filename in filelist:
        imagetemps.append(pyfits.open(filename)[0].header["ccdtemp"])
    median_temp = np.median(imagetemps)

    min_temp_diff = 1000
    min_dark = None

    for index in range(1, i):
        filename = "masterdark0_" + str(index) + ".fits"
        hdu=pyfits.open(filename)[0]
        darktemp = hdu.header["tempmed"]
        diff = abs(darktemp - median_temp)
        if diff < min_temp_diff:
            min_temp_diff = diff
            min_dark = index
    return min_dark

outputdir = sys.argv[1]
filelist = sys.argv[2:]
if len(filelist)==0: exit(-1)

masterdark_index = choose_best_dark(filelist)

tmpfolder = "tmp/flat" + str(random.randint(0,10000000)) + "/"
os.mkdir(tmpfolder)

utils.split_by_channel(filelist, tmpfolder, False)

#for each channel, invoke do_masterflat.py
#command_queue = []
for ch in range(4):
    biasname = "masterbias" + str(ch) + ".fits"
    darkname = "masterdark" + str(ch) + "_" + str(masterdark_index) + ".fits"
    flatname = outputdir + "masterflat" + str(ch) + ".fits"

    call_str = "python ~hatuser/HATpipebin/include/HATpipepy/Actions/do_masterflat.py --master-bias=" + biasname + " --master-dark=" + darkname + " " + \
        "HATnet " + tmpfolder + "*_" + str(ch) + ".fits" + " -o " + flatname
    os.system(call_str)
#    command_queue.append(call_str)

#utils.run_all(command_queue)
