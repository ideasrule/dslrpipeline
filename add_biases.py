import pyfits

def add(filelist):
    master = None
    for filename in filelist:
        hdu = pyfits.open(filename)[0]
        if master==None: master = hdu
        else: master.data += hdu.data
    master.data /= len(filelist)
    master.writeto("master_dark.fits")


f = open("file_objs")
filelist = []
for line in f:
    elements = line.split()
    if len(elements) != 2: continue
    if elements[1]=="nightdark" or elements[1]=="mondark" or elements[1]=="dark": filelist.append(elements[0])

add(filelist)
print filelist
