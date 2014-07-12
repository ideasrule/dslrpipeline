import pyfits
import sys
import os.path

for filename in sys.argv[1:]:
    f=open(filename, "r")
    fitsfile=filename[:-11] + ".fits"
    print fitsfile

    hdu=pyfits.open(fitsfile)[0]
    JD=hdu.header["JD"]
    Z = hdu.header["Z"]
    seqnum = os.path.basename(filename)[3:-7]
    outf = open(filename[:-11] + ".asubpixphot", "w")

    all_lines=[]
    skip = False
    for line in f:
     #   print line
        line=line.strip()
        pos = line.find("nan")
        if pos != -1 and line[0] != "#": continue
        if line[0]=="#": 
            #already correct? skip
            if line.find("JD") != -1 and line.find("Z") != -1:
                skip=True
                break
            line += "\t JD \t\t Z \t seqnum\n"
        else:
            line += "\t" + str(JD) + "\t" + str(Z) + "\t" + str(seqnum) + "\n"
        all_lines.append(line)
    if skip: continue
    for line in all_lines:
        outf.write(line)
    outf.close()
    f.close()
