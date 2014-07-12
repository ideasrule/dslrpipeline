import sys
import pyfits
import numpy as np

masterbias = sys.argv[1]

biasdata = pyfits.open(masterbias)[0].data[0]
filelist = sys.argv[2:]
temps = []

newfilelist = []
for filename in filelist:
#    hdu=pyfits.open(filename)[0]
 #   if int(hdu.header["aborted"])==0 and int(hdu.header["exptime"])==180:
  #      temps.append(hdu.header["CCDTEMP"])
    newfilelist.append(filename)

print temps
#reftemp = np.median(temps) + 273.15 #temperature in K
#print reftemp
sum_darks = []

for filename in newfilelist:
    print filename
    hdu=pyfits.open(filename)[0]
    Tabs = hdu.header["CCDTEMP"] + 273.15
#    factor = np.exp(12800*(1./Tabs - 1./reftemp))
#    corrected = (hdu.data - biasdata)*factor
#    print hdu.data.shape
 #   print biasdata.shape
    corrected = hdu.data[0:1734,71:2673] - biasdata
    if len(sum_darks)==0: sum_darks = corrected
    else: sum_darks += corrected
    
outhdu = pyfits.PrimaryHDU(sum_darks/len(newfilelist))
outhdu.writeto("master_dark_uncor.fits")



#convert all frames to
