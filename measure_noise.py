#!/usr/bin/python
import pyfits
import numpy as np
import sys
import os
import re
from scipy import polyfit, polyval
from scipy.stats import linregress
import matplotlib.pyplot as plt

def reject_outliers(data, m=5.):
    data = data.flatten()
    print "avg,std", np.mean(data), np.std(data), np.max(data)
    return data[abs(data - np.mean(data)) < m*np.std(data)]

def analyze_twoim(filename1, filename2):
    hdu1 = pyfits.open(filename1)[0]
    hdu2 = pyfits.open(filename2)[0]
  #  data1 = hdu1.data[764:964,1273:1473] - np.median(hdu1.data[0:1730,0:29])
   # data2 = hdu2.data[764:964,1273:1473] - np.median(hdu2.data[0:1730,0:29])
    data1 = hdu1.data[0:1728,73:-1]
    data2 = hdu2.data[0:1728,73:-1]

    avg = (data1 + data2)/2 - 2048.
    diff = data1 - data2
              
    return (np.mean(avg), np.var(diff)/2.)

maxexp = int(sys.argv[1])

for ch in ['R', 'G1', 'G2', 'B']:
    print ch
    means = []
    vars = []
    for exp in range(1, maxexp, 2):
        filelist = []
        mean_collection = []
        var_collection = []
        for i in range(9):
            #consider images i and i+1
            filename1 = "gain" + str(i) + "_" + str(exp) + ch + ".fits"
            filename2 = "gain" + str(i+1) + "_" + str(exp) + ch + ".fits"
            (mean,var) = analyze_twoim(filename1, filename2)
            mean_collection.append(mean)
            var_collection.append(var)

        (mean,var) = np.median(mean_collection), np.median(var_collection)
#        print mean
#        print mean,var
        means.append(mean)
        vars.append(var)
    slope, intercept, r, pvalue, see = linregress(means, vars)
    mx = np.mean(means)
    sx2 = ((means-mx)**2).sum()
    sd_intercept = see*np.sqrt(1./len(means) + 1.*mx*mx/sx2)
    sd_slope = see*np.sqrt(1./sx2)
    print 1./slope,intercept,sd_slope,sd_intercept
