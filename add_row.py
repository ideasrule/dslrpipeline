import sys
import numpy as np

f=open(sys.argv[1])
for line in f:
    elements=line.split()
    if len(elements)==0: break
    temp = float(elements[1])
    factor = np.exp(-12800/(273.15 + temp))*1e19
    print factor,elements[2],elements[3],elements[4],elements[5]
    
