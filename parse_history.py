from os.path import basename
from glob import glob
import pyfits

f = open("output")
lines = []
for line in f:
    line = line.strip()
    if len(lines)>0 and lines[-1][-1]=="\\":
        lines[-1] = lines[-1][:-1] + line[2:]
    else:
        lines.append(line)

elements =  lines[3].split()
elements = elements[6:-2]
for e in elements:    
    name = basename(e)
    name = name[:-6] + "2" + ".fits"
    filename = glob("10-2014*/" + name)[0]
    print filename
#    print pyfits.open(filename)[0].header["ccdtemp"]
#    print name
