#!/bin/bash

source define_frame_number_lists.sh

CMPOS=$1
NPROC=$2

CMPOS_M1=$((CMPOS-1))
#FRAME_NUMBERS=(${fnum_lists[$CMPOS_M1]})
FRAME_NUMBERS=( $(for fnum in ${fnum_lists[$CMPOS_M1]}; do if [ "$(fiheader -g MNTSTATE ../RED/2-${fnum}_${CMPOS}.fits)" == "T" ] ; then echo $fnum; fi; done) )
FRAME_COUNT=${#FRAME_NUMBERS[@]}
CAT="../cat_${CMPOS}.2mass"
for ((offset=0; offset<NPROC; offset++)); do 
	for ((frame_ind=offset; frame_ind<FRAME_COUNT; frame_ind+=NPROC)) ; do
		fnum=${FRAME_NUMBERS[$frame_ind]}
		trans=../ASTROM/2-${fnum}_${CMPOS}.trans
		fits=../RED/2-${fnum}_${CMPOS}.fits
		outfname=../SUBPIX_SDK/2-${fnum}_${CMPOS}.spsdk
		ra=$(grep '^# 2MASS' $trans |awk '{print $3;}')
		dec=$(grep '^# 2MASS' $trans |awk '{print $4;}')
		grtrans --input ${CAT} --col-radec 2,3 --wcs \
			arc,ra=${ra},dec=${dec},degrees|\
		grtrans --input - --col-xy 2,3 --input-transformation $trans |\
		awk '{if($2>0 && $2<4096 && $3>0 && $3<4096) \
			printf "%-11s %16.5f %16.5f\n", $1, $2, $3;}'|\
		~/SVN/HATpipe/source/subpixel_sensitivity/FitPSF $fits \
			--input-columns id,x,y \
			--output-columns id,x,y,S,D,K,A,flux,flux_err,bg,bg_err,nbgpix,flag,sn,chi2,npix \
			--fit-order -1 --max-sources 10000 --max-sat-frac 0 --maxS 10\
			--source-assignment \
			${outfname%.spsdk}_source_assignment.fits > $outfname
		echo "done with $fnum sub-pixel S, D, K fit."
	done &
done
