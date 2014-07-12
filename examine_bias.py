import pyfits
import numpy as np
import sys
import matplotlib.pyplot as plt

def reject_outliers(data, m=5.):
    data = data.flatten()
    return data[abs(data - np.mean(data)) < m*np.std(data)]

def plot(filelist, ch, bias_data):
    means=[]
    variances=[]
    JDs=[]
    temps=[]
    print "filename,JD,temps,mean_dark,im_mean,dark_var,imvar,exptime"
    for filename in filelist:
        hdu = pyfits.open(filename)[0]
#        hdu.data -= bias_data
        JDs.append(float(hdu.header["JD"]))
        temps.append(float(hdu.header["AMBTEMP"]))
        imregion = reject_outliers(hdu.data[ch][0:1728,73:])
#        imregion = reject_outliers(imregion)
#        darkregion = reject_outliers(hdu.data[ch][1734:1741,73:])
        darkregion = reject_outliers(hdu.data[ch][0:1728,32:65])
 #       darkregion = reject_outliers(darkregion)

        means.append(np.mean(imregion))
        variances.append(np.var(imregion))
        dark_variance = np.var(darkregion)
        print filename, JDs[-1],temps[-1],np.mean(darkregion),means[-1],dark_variance,variances[-1],hdu.header["EXPTIME"]
#    plt.plot(JDs, variances)
#    plt.savefig("JD_variances.png")
 #   plt.show()

f = open("file_objs")
bias_data = pyfits.open("master_bias.fits")[0].data
filelist = []
for line in f:
    elements = line.split()
    if len(elements) != 2: continue
    if elements[1]=="nightdark" or elements[1]=="mondark" or \
            elements[1]=="dark": filelist.append(elements[0])

ch = int(sys.argv[1])
plot(filelist, ch, bias_data)

