import pyfits
import sys
import os

f=open(sys.argv[1])

names=[]
ras=[]
decs=[]

for line in f:
    elements = line.split()
 #   print elements
    if elements[0]=='#': continue
    names.append(elements[0])
    ras.append(elements[1])
    decs.append(elements[2])

col1 = pyfits.Column(name='name', format='20A', array=names)
col2 = pyfits.Column(name='RA', format='E', array=ras)
col3 = pyfits.Column(name='Dec', format='E', array=decs)
cols=pyfits.ColDefs([col1,col2,col3])
tbhdu=pyfits.new_table(cols)
tbhdu.writeto(sys.argv[2], clobber=True)
