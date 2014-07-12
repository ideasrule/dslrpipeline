import pyfits
import sys
import os

f=open(sys.argv[1])

ix=[]
iy=[]
mags=[]

for line in f:
    elements = line.split()
    ix.append(float(elements[0]))
    iy.append(float(elements[1]))
    mags.append(float(elements[2]))
    print elements

#print ix,iy,mags
col1 = pyfits.Column(name='x', format='E', array=ix)
col2 = pyfits.Column(name='y', format='E', array=iy)
col3 = pyfits.Column(name='mag', format='E', array=mags)
cols=pyfits.ColDefs([col1,col2,col3])
tbhdu=pyfits.new_table(cols)
tbhdu.writeto(sys.argv[2], clobber=True)
