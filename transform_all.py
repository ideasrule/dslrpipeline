from multiprocessing import Pool
import sys
import pyfits

def transform_once(filename):
    header = pyfits.open(filename)[0].header
    tmpfile = filename[:-5] + ".tmp"
    print tmpfile
    if header["object"] == "G09363730_180":
        print "180"
    else:
        print "113"

if __name__=='__main__':
    filelist = sys.argv[1:]
    workers = Pool(40)
    workers.map(transform_once, filelist)
    workers.close()
    workers.join()
