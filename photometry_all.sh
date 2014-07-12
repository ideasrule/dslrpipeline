#!/bin/bash

for dir in $1
do
	files=$dir*_[0123].fits
	for f in $files
	do
	    if [ ! -f ${f%.fits}.solved ]; then continue; fi
	    catfile=${f%.fits}.catalog
	    photfile=${f%.fits}.phot
	    tmpfile=${f%.fits}.tmp
	    echo $photfile
	    awk '{print $1" "$2" "$3-0.5" "$4-0.5}' $catfile > $tmpfile
	    fiphot --format IXY,FfBbs --apertures 1.2:6:7,1.4:6:7,1.6:6:7,1.8:6:7,2:6:7,2.2:6:7,2.4:6:7,2.6:6:7,2.8:6:7,3:6:7 -il $tmpfile --output $photfile --input $f --gain 2.1 --aperture-mask-ignore saturated --nan-string NaN --single-background 5 -k --sky-fit mode,sigma=3,iterations=2 --disjoint-radius 2 --col-id 1 --col-xy 3,4
	    sed -i '1i#ID               x              y         BG            BGerr    Flux[0]  FluxErr[0] S[0]  Flux[1]  FluxErr[1] S[1] Flux[2]  FluxErr[2] S[2]   Flux[3]  FluxErr[3] S[3]  Flux[4]  FluxErr[4] S[4]  Flux[5]  FluxErr[5] S[5]   Flux[6]  FluxErr[6] S[6]   Flux[7]  FluxErr[7] S[7]   Flux[8]  FluxErr[8] S[8]   Flux[9]  FluxErr[9] S[9]' $photfile
	    grep -v NaN $photfile > $tmpfile
	    mv $tmpfile $photfile
	done
done