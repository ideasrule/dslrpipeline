import sys
import pyfits
import random
import numpy as np
import os

def choose_best_dark(filename, darks):
    #returns index of best dark to use
    imagetemp = int(pyfits.open(filename)[0].header["ccdtemp"])

    min_temp_diff = 1000
    min_dark = None
    for i in range(len(darks)):
        filename = darks[i]
        hdu=pyfits.open(filename)[0]
        darktemp = int(hdu.header["tempmed"])
        diff = abs(darktemp - imagetemp)
        if diff < min_temp_diff:
            min_temp_diff = diff
            min_dark = i
    return min_dark

def run_ficalib(filelist, dark, bias, flat):
    darkdata = pyfits.open(dark)[0].data
    biasdata = pyfits.open(bias)[0].data
    flatdata = pyfits.open(flat)[0].data

    tmpfolder = "tmp/calib" + str(random.randint(0,10000000)) + "/"
    print tmpfolder
    os.mkdir(tmpfolder)

    ficalib_inputs = []
    ficalib_outputs = []

    for ch in range(4):
        input_str = ""
        output_str = ""
        for filename in filelist:
            basename = os.path.basename(filename)[:-6]
            newname = tmpfolder + basename + str(ch) + ".fits"
            call_str = "mcolor.py -c " + str(ch) + " -o " + newname + \
                " " + filename
            os.system(call_str)
            input_str += " " + newname + " "
            output_str += " RED/" + filename[:-6] + str(ch) + ".fits"
        ficalib_inputs.append(input_str)
        ficalib_outputs.append(output_str)
    #splitting complete!
    print "Splitting complete!"
    print ficalib_inputs
    print ficalib_outputs

    for ch in range(4):
        biasname = tmpfolder + "masterbias" + str(ch) + ".fits"
        newhdu=pyfits.PrimaryHDU(biasdata[ch])
        newhdu.writeto(biasname, clobber=True)

        darkname = tmpfolder + "masterdark" + str(ch) + ".fits"
        newhdu=pyfits.PrimaryHDU(darkdata[ch])
        newhdu.writeto(darkname)

        flatname = tmpfolder + "masterflat" + str(ch) + ".fits"
        newhdu = pyfits.PrimaryHDU(flatdata[ch])
        newhdu.writeto(flatname)

        call_str = "ficalib -iB " + biasname + " -iD " + darkname + " -iF " +\
            flatname + " " + ficalib_inputs[ch] + " -o " + ficalib_outputs[ch]\
            + " " + "--saturation 14500 --image 71:0:2671:1730 --trim -g 2.1"
        print call_str
        os.system(call_str)

masterbias = sys.argv[1]
masterflat = sys.argv[2]
filelist = sys.argv[3:]

images_by_dark=[]
darks=[]

i=1
while 1==1:
    filename="masterdark_3d_" + str(i) + ".fits"
    if os.path.isfile(filename): darks.append(filename)
    else: break
    i += 1

for i in range(len(darks)):
    images_by_dark.append([])

#now group by dark
for filename in filelist:
    dark_index = choose_best_dark(filename, darks)
    images_by_dark[dark_index].append(filename)

for i in range(len(images_by_dark)):
    run_ficalib(images_by_dark[i], darks[i], masterbias, masterflat)
