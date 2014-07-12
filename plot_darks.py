import matplotlib.pyplot as plt
import sys

f = open(sys.argv[1])

time=[]
vals=[]

time_high=[]
time_low=[]
vals_low=[]
vals_high=[]
#vals_dark=[]

for line in f:
    elements=line.split()
    if len(elements)==1: continue
#    if len(elements) != 7: continue
    if float(elements[7]) > 179: #exposure time is 180
        time.append(float(elements[1]))
        if float(elements[3]) > 2048.8:
            vals_high.append(float(elements[6]))
            time_high.append(float(elements[1]))
        else:
            vals_low.append(float(elements[6]))
            time_low.append(float(elements[1]))
#        vals_dark.append(float(elements[6]))

#plt.plot(time, 
plt.plot(time_low, vals_low, ".")
plt.plot(time_high, vals_high, ".")

plt.legend()
plt.show()
