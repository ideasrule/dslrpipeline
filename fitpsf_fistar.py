import numpy as np
import sys
import math

def fit(filename, bcol=3, width=2601., height=1734., order=5):
    matrix = []
    vector = []

    for line in open(filename):
        elements=line.split()
        xscale=(float(elements[0])-width/2.)/(width/2.)
        yscale=(float(elements[1])-height/2.)/(height/2.)
        vector.append(float(elements[bcol]))

        matrix_row=[]
        for totpow in range(order+1):
            for ypow in range(totpow+1):
                xpow = totpow - ypow
                matrix_row.append(xscale**xpow*yscale**ypow/math.factorial(xpow)/math.factorial(ypow))
        matrix.append(matrix_row)
    sol,res,rank,s = np.linalg.lstsq(matrix, vector)
    return sol

def writepsf(filename):
    #filename should be .phot
    infilename = filename[:-5] + ".fistar"
    outfile = open(filename[:-5] + ".psf", "w")
    outfile.write("#PSF coefficients: ")
    arr = []
    for bcol in range(3,6):
        result = fit(infilename,bcol)
        for r in result:
            arr.append(r)
    outfile.write(' '.join(map(str, arr)))
    outfile.close()

for photfile in sys.argv[1:]:
    writepsf(photfile)

