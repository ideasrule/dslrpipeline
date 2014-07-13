import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.linalg import lstsq
from mpl_toolkits.mplot3d.axes3d import Axes3D
import random
import scipy.ndimage.filters as f

starname=sys.argv[1]

all_dict=dict()
counts=dict()
medians=[]
mediandevs=[]
weights=[]

def compute_dev(ch, filename, m=5., order=4):
    JDs=[]
    vals = []

    for line in open(filename):
        elements=line.split()
        val=float(elements[2])
        JDs.append(float(elements[1]))
        vals.append(val)
#        if len(vals) > 1100: break

    JDs=np.array(JDs)
    vals=np.array(vals)
    med = np.median(vals)
    deviation = np.median(np.abs(vals - med))
    medians.append(med)
    mediandevs.append(deviation)
    weights.append(1./deviation**2)
    for i in range(len(JDs)):
        if JDs[i] not in all_dict:
            all_dict[JDs[i]] = vals[i]
            counts[JDs[i]] = 1
        else:
            all_dict[JDs[i]] += vals[i]
            counts[JDs[i]] += 1

starname = sys.argv[1]
mode = sys.argv[2]

for i in range(4):
    filename = mode + str(i) + "_" + starname + ".epdlc"
    compute_dev(i, filename)

#print counts.values()
JDs = np.array([key for key in all_dict if counts[key]==4])
vals = np.array([all_dict[key] for key in all_dict if counts[key]==4])
vals /= 4
#vals /= np.sum(weights)

med = np.median(vals)
deviation = np.median(np.abs(vals - med))
medians.append(med)
mediandevs.append(deviation)

for i in range(len(medians)):
    sys.stdout.write(str(medians[i]) + " " + str(mediandevs[i]) + " ")
sys.stdout.write("\n")



