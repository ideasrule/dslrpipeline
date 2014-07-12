import pyfits
import sys

hdulist=pyfits.open(sys.argv[1])
tbdata = hdulist[1].data

for row in tbdata:
    print row[0], ' '.join(map(str, row[1:]))
