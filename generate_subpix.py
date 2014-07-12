import pyfits

pattern = [[2,0],[3,1]]
for ch in range(4):
    sensitivity = [[0,0], [0,0]]
    for y in range(2):
        for x in range(2):
            if ch == pattern[y][x]: sensitivity[y][x] = 4

    hdu = pyfits.PrimaryHDU(sensitivity)
    hdu.writeto("subpix" + str(ch) + ".fits", clobber=True)
