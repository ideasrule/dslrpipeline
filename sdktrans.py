import sys
import math

def transform(coeffs, x, y):
    order = 5
    i = 0
    result = 0
    width = 2601.
    height = 1734.
    for totpow in range(order+1):
        for ypow in range(totpow+1):
            xpow = totpow - ypow
  #          print i
            xscale = (x-width/2.)/(width/2.)
            yscale = (y-height/2.)/(height/2.)
            result += coeffs[i]*xscale**xpow*yscale**ypow/math.factorial(xpow)/math.factorial(ypow)
            i += 1
    return result

#f = open(sys.argv[1])
#line = f.readline()
#elements = line.split()[2:]
#elements = [ float(x) for x in elements ]
#print len(elements)/3

psffile = open(sys.argv[1])
line = psffile.readline()
coeffs = line.split()[2:]
coeffs = [float(x) for x in coeffs]

photfile = open(sys.argv[2])
print "S D K"
for line in photfile:
    if line[0]=='#': continue
    elements=line.split()
    x=float(elements[1])
    y=float(elements[2])
#    print x,y
    s=transform(coeffs[0:21],x,y)
    d=transform(coeffs[21:42],x,y)
    k=transform(coeffs[42:63],x,y)
    print s,d,k


