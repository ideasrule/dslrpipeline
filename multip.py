#!/usr/bin/python
from multiprocessing import Pool
import sys
import numpy as np
import os
from optparse import OptionParser

class CallObj:
    def __init__(self,command, arg):
        self.command = command
        self.arg = arg

def run_program(call_obj):
    call_str = call_obj.command + " " + call_obj.arg
    print "Starting: " + call_str
    os.system(call_str)

if __name__=='__main__':
    parser = OptionParser()
    parser.add_option("-p", dest="processes", default=None, type=int)

    (options, args) = parser.parse_args()
    
    if options.processes == None:
        workers = Pool()
    else:
        workers = Pool(options.processes)

    command = args[0]
    call_objs = [CallObj(command,arg) for arg in args[1:]]

    workers.map(run_program, call_objs)
    workers.close()
    workers.join()
