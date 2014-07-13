import sys
import pyfits
import random
import numpy as np
import os
import utils

def run_command(command):
    os.system(command)

def choose_best_dark(filename, numdarks):
    #returns index of best dark to use
    imagetemp = int(pyfits.open(filename)[0].header["ccdtemp"])

    min_temp_diff = 1000
    min_dark = None
    for i in range(numdarks):
        filename = "masterdark0_" + str(i+1) + ".fits"
        hdu=pyfits.open(filename)[0]
        darktemp = int(hdu.header["tempmed"])
        diff = abs(darktemp - imagetemp)
        if diff < min_temp_diff:
            min_temp_diff = diff
            min_dark = i
    return min_dark

def run_ficalib(filelist, darkindex, masterflatdir):
    tmpfolder = "tmp/calib" + str(random.randint(0,10000000)) + "/"
    print tmpfolder
    os.mkdir(tmpfolder)

    ficalib_inputs = []
    ficalib_outputs = []
    command_queue = []
    for ch in range(4):
        input_str = ""
        output_str = ""
        for filename in filelist:
            basename = os.path.basename(filename)[:-6]
            newname = tmpfolder + basename + str(ch) + ".fits"
            call_str = "mcolor.py -c " + str(ch) + " -o " + newname + \
                " " + filename
            command_queue.append(call_str)
            input_str += " " + newname + " "
            output_str += " RED/" + filename[:-6] + str(ch) + ".fits"
        ficalib_inputs.append(input_str)
        ficalib_outputs.append(output_str)

    utils.run_all(command_queue)
    command_queue = []

    print "Splitting complete!"
    for ch in range(4):
        biasname = "masterbias" + str(ch) + ".fits"
        darkname = "masterdark" + str(ch) + "_" + str(darkindex+1) + ".fits"
        flatname = masterflatdir + "masterflat" + str(ch) + ".fits"

        call_str = "ficalib -iB " + biasname + " -iD " + darkname + " -iF " +\
            flatname + " " + ficalib_inputs[ch] + " -o " + ficalib_outputs[ch]\
            + " " + "--saturation 14500 --image 71:0:2671:1730 --trim -g 2.1"
        command_queue.append(call_str)
    utils.run_all(command_queue)

masterflatdir = sys.argv[1]
filelist = sys.argv[2:]

images_by_dark=[]
darks=[]

i=1
while 1==1:
    filename="masterdark0_" + str(i) + ".fits"
    if not os.path.isfile(filename): break
    i += 1

numdarks = i-1

for i in range(numdarks):
    images_by_dark.append([])

#now group by dark
for filename in filelist:
    dark_index = choose_best_dark(filename, numdarks)
    images_by_dark[dark_index].append(filename)

for i in range(len(images_by_dark)):
    run_ficalib(images_by_dark[i], i, masterflatdir)
