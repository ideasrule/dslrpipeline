import pyfits
import sys

hdu1 = pyfits.open(sys.argv[1])[0]
hdu2 = pyfits.open(sys.argv[2])[0]

hdu1.data = hdu2.data/hdu1.data
hdu1.writeto("diff.fits", clobber=True)
