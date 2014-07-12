import sys
import os

f=open(sys.argv[1])
for line in f:
    elements=line.split()
    call_str = "grep " + elements[0] + " 2masscat_108 | awk '{print $15}'"
    os.system(call_str)
#    print call_str
