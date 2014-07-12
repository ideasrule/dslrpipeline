import pyfits
import sys
import os.path

for filename in sys.argv[1:]:
    f=open(filename, "r+")
    seqnum = os.path.basename(filename)[3:-7]

    all_lines=[]
    skip = False
    for line in f:
        line = line.strip()
        if line[0]=="#":             
            #already correct? skip
            if line.find("seqnum") != -1:
                skip=True
                break
            line += "\t seqnum\n"
        else:
            line += "\t" + str(seqnum) + "\n"
        all_lines.append(line)
    if skip: continue
    f.seek(0)
    for line in all_lines:
        f.write(line)
    f.close()
