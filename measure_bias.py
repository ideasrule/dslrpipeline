#!/usr/bin/python
import pyfits
import numpy as np
import sys
import os
import re
import matplotlib.pyplot as plt
from scipy import polyfit, polyval
from scipy.stats import linregress

def reject_outliers(data, m=10.):
    data = data.flatten()
    to_ret = data[abs(data - np.mean(data)) < m*np.std(data)]
    return to_ret

def get_biases(filelist, ch):
    sum_sqr = []
    sum = []
    variance = []
    hdu = None
    for filename in filelist:
        hdu = pyfits.open(filename)[0]
        data = hdu.data[ch]
        data = data[0:1728,73:]

        if len(sum) == 0:
            sum = data
            sum_sqr = data**2
        else:
            sum += data
            sum_sqr += data**2

    N = float(len(filelist))
    mean = sum/N
    variance = (sum_sqr - N*mean**2)/(N-1)
    return (np.mean(mean.flatten()), np.mean(variance.flatten()))

filelist = sys.argv[1:]

for ch in range(4):
    (mean,var) = get_biases(filelist, ch)
    print mean,var


