import matplotlib.pyplot as plt
import sys

f = open(sys.argv[1])
time=dict()
vals=dict()
vals_dark=dict()

for line in f:
    elements=line.split()
    if len(elements) != 7: continue
#    imtype = elements[3]
    imtype="dark"
    if imtype not in time:
        time[imtype]=[]
        vals[imtype]=[]
        vals_dark[imtype]=[]
    if float(elements[6]) > 179: #exposure time is 180
        time[imtype].append(float(elements[1]))
        vals[imtype].append(float(elements[4]))
        vals_dark[imtype].append(float(elements[5]))

#colors=['r','g','b','y','B']
#i=0
for key in time:
    plt.plot(time[key], vals[key], ".", label=key)
    plt.plot(time[key], vals_dark[key], ".", label=key)

plt.legend()
#plt.scatter(time,vals)
plt.show()
