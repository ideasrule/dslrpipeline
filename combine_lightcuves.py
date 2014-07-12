import matplotlib.pyplot as plt
from scipy.linalg import lstsq
import sys
import numpy as np

f=open(sys.argv[1])
matrix=[]
b=[]

vals0=[]
vals1=[]
vals2=[]
vals3=[]

count=0
final=[]
for line in f:
    count += 1
#    if count > 40: break
    elements=line.split()
    ch0=10**(-float(elements[0])/2.5)
    ch1=10**(-float(elements[1])/2.5)
    ch2=10**(-float(elements[2])/2.5)
    ch3=10**(-float(elements[3])/2.5)
    if np.isnan(ch0) or np.isnan(ch1) or np.isnan(ch2) or np.isnan(ch3): continue
    vals0.append(ch0)
    vals1.append(ch1)
    vals2.append(ch2)
    vals3.append(ch3)
#    print ch0,ch1,ch2,ch3
    matrix.append([ch1, ch2, -1])
    b.append(-ch0 - ch3)
    combined=ch0 + ch3 #+ 0.305*ch1 + 0.153*ch2
    final.append(combined)
#    print ch0 + ch3
print lstsq(matrix, b)
plt.plot(vals0, label="0")
plt.plot(vals1,label="1")
plt.plot(vals2,label="2")
plt.plot(vals3,label="3")
plt.plot(final, label="final")
plt.legend()
plt.show()

