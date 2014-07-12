from multiprocessing import Pool
import sys
import numpy as np
import os
from optparse import OptionParser

def run_subpix(filename):
    #filename must be a .phot file
    ch = int(filename[-6])
    print filename, ch
    subpixfile = "subpix" + str(ch) + ".fits"
    offsets = 0.25*np.array([[1,-1], [1,1], [-1,-1], [-1,1]])

    photfile = filename
    fitsfile = filename[:-5] + ".fits"
    psffile = filename[:-5] + ".psf"
    outfile = filename[:-5] + ".subpixphot"

    call_str = "awk '{print $1,$2,$3,$4,$5,$6,$7,$10}' " + psffile + " | ./SubPixPhot " + fitsfile + " --input-columns id,x,y,S,D,K,A,bg --output-columns id,x,y,S,D,K,bg,bg_err,flux,flux_err,flag --subpix=" + subpixfile + " --gain=2.1 -a 1.6 -a 1.8 -a 2.0 -a 2.2 -a 2.4 -a 2.6 -a 2.8 -a 3.0 -a 3.2 -a 3.4 -o " + outfile

    print call_str
    os.system(call_str)

if __name__=='__main__':
    parser = OptionParser()
    parser.add_option("-p", dest="processes", default=None, type=int)

    (options, args) = parser.parse_args()
    
    start = int(args[0])
    end = int(args[1])
    filelist=args[2:]
    if end > len(filelist): end = len(filelist)
    filelist = filelist[start:end]

    if options.processes == None:
        workers = Pool()
    else:
        workers = Pool(options.processes)
    workers.map(run_subpix, filelist)
    workers.close()
    workers.join()
