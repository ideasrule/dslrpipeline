import sys
import random
import os
import numpy as np
import pyfits
import utils

from utils import split_by_channel
from utils import run_command

filelist = sys.argv[1:]

tmpfolder = "tmp/bias" + str(random.randint(0,10000000))
os.mkdir(tmpfolder)

utils.split_by_channel(filelist, tmpfolder)

#for each channel, invoke do_masterbias.py
command_queue=[]

for ch in range(4):
    biasname = "masterbias" + str(ch) + ".fits"
    call_str = "python ~hatuser/HATpipebin/include/HATpipepy/Actions/do_masterbias.py HATnet " + tmpfolder + "/*_" + str(ch) + ".fits" + " -o " + biasname
    command_queue.append(call_str)

utils.run_all(command_queue)
