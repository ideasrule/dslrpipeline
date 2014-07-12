import pyfits
import sys
import os.path

for filename in sys.argv[1:]:
    f=open(filename, "r+")
    seqnum = os.path.basename(filename)[3:-7]

    all_lines=[]
    for line in f:        
        line = line.strip()
        if line[0]=='T' and line[1]=='-': line = "HA" + line
        elements=line.split()
#        print elements
        all_lines.append(' '.join(elements[0:55]))
        print all_lines[-1]
#        print line
#    print all_lines
#    f.seek(0)
 #   for line in all_lines:
  #      f.write(line)
   # f.close()
