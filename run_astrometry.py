from multiprocessing import Pool
import pyfits
import sys
import os
from optparse import OptionParser

def run_once(filename):
    #first, need to run fistar
    call_str = "fistar --format x,y,flux,s,d,k -s flux -o " + filename[:-5] +\
        ".fistar " + filename
    os.system(call_str)
    xs=[]
    ys=[]
    fluxes=[]
    count = 0
    for line in open(filename[:-5] + ".fistar"):
        if count > 6000: break
        elements = line.split()
        xs.append(elements[0])
        ys.append(elements[1])
        fluxes.append(elements[2])
        count += 1
    col1 = pyfits.Column(name='x', format='E', array=xs)
    col2 = pyfits.Column(name='y', format='E', array=ys)
    col3 = pyfits.Column(name='fluxes', format='E', array=fluxes)
    cols=pyfits.ColDefs([col1,col2,col3])
    tbhdu=pyfits.new_table(cols)
    tbhdu.writeto(filename[:-5], clobber=True)

    header=pyfits.open(filename)[0].header
    #if this isn't a valid frame, do nothing
    if header["imagetyp"] != "object" or header["aborted"] == 1 or \
            header["exptime"] < 179:
        print filename + " is not a valid frame!"
        return
    ra=header["ra"]*15
    dec=header["dec"]
    radius = 7 #look less than 7 degrees from this RA,Dec
    fov=float(header["fov"])
    width=str(header["naxis1"])
    height=str(header["naxis2"])

    call_str = "solve-field --overwrite -p --scale-low " + str(fov-1) + " --scale-high " + str(fov+1) + " --scale-units degwidth --parity pos -N none --ra " + str(ra) + " --dec " + str(dec) + " --radius " + str(radius) + " -t 3 -w " + \
        width + " -e " + height + " " + filename[:-5]
    os.system(call_str)

    #now "verify WCS" to fit to all stars in all index files
    #q is min. quad size--0.01 to find stars in all index files
    if not os.path.isfile(filename[:-5] + ".solved"): return
    call_str = "solve-field -p -N none --continue -t 3 -q 0.01 -V " + \
        filename[:-5] + ".wcs " + " -w " + width + " -e " + height + " " + \
        filename[:-5]
    print call_str
    os.system(call_str)

if __name__=='__main__':
    filelist = sys.argv[1:]
    workers = Pool(40)
    workers.map(run_once, filelist)
    workers.close()
    workers.join()
