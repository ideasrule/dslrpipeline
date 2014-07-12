import sys
import pyfits
import matplotlib.pyplot as plt
import numpy as np

hdulist = pyfits.open(sys.argv[1])
tbdata = hdulist[1].data

xs=[]
ys=[]
index_xs=[]
index_ys=[]

for row in tbdata:
    xs.append(row[0])
    ys.append(row[1])
    index_xs.append(row[4])
    index_ys.append(row[5])

xs=np.array(xs)
ys=np.array(ys)
index_xs=np.array(index_xs)
index_ys=np.array(index_ys)

diffx = index_xs - xs
diffy = index_ys - ys
#print len(diffx), len(diffy)
print np.std(diffx), np.std(diffy)
#diffx = diffx[abs(diffx) < 1]
#diffy = diffy[abs(diffy) < 1]

print np.std(diffx), np.std(diffy)
#print len(diffx), len(diffy)

plt.scatter(xs, diffx)
plt.xlabel("X (pixels)")
plt.ylabel("Diff between catalogue and detected star x")
plt.figure()
plt.scatter(xs, diffy)
plt.xlabel("X (pixels)")
plt.ylabel("Diff between catalogue and detected star y")
plt.show()
