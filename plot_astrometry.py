#Take _0.fits
#Find _0.catalog, _1.catalog, etc
#grep for the star, awk

import matplotlib.pyplot as plt
import sys
import os

def find_coords(filename):

    call_str = "grep HAT-181-0000237 " + filename + " | awk '{print $3,$4;}'"\
        + ">tmpfile"
#    call_str = "grep HAT-180-0000507 " + filename + " | awk '{print $3,$4;}'"\
#        + ">tmpfile"
    os.system(call_str)
    f = open("tmpfile")
    line = f.readline()
    elements=line.split()
    return [float(elements[0]), float(elements[1])]


xs=[]
ys=[]
colors=['g','b','r','black']
for ch in range(4):
    basename = sys.argv[1][:-6]
    newname = basename + str(ch) + ".catalog"
    coords = find_coords(newname)
    plt.scatter(coords[1],coords[0],c=colors[ch])
#    xs.append(coords[1])
 #   ys.append(coords[0])


#plt.scatter(xs,ys)
plt.show()

#print find_coords(sys.argv[1])
