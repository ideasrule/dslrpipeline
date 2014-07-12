import pyfits
import numpy as np
import sys

#def reject_outliers(data, m=5.):
 #   data = data.flatten()
  #  return data[abs(data - np.mean(data)) < m*np.std(data)]

hdu=pyfits.open(sys.argv[1])[0]
#shape = hdu.data.shape

data = hdu.data[3]

#data = reject_outliers(hdu.data[0][0:1728,73:]).reshape(shape)
result = np.zeros(data.shape)
print data.shape
for i in range(data.shape[0]):
    print i
    for j in range(data.shape[1]):
        result[i][j] = np.var(data[i:i+25,j:j+25])

hdu.data = result
hdu.writeto(sys.argv[2])

