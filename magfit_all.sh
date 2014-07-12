#!/bin/bash

D20files=$(awk '{if($2=="object" && $3=="G09363730_180" && $4=="180" && $5=="0" && $6=="D20")print $1}' database_file)

Tfiles=$(awk '{if($2=="object" && $3=="G09363730_180" && $4=="180" && $5=="0" && $6=="T")print $1}' database_file)

rm tmp[DT]$1

for f in $D20files
do
    fullname=RED/${f%_2.fits}_$1.fits
    if [ -f ${fullname%.fits}.solved ]
    then 
	echo ${fullname%.fits}.phot $fullname >> tmpD$1
    fi
done

for f in $Tfiles
do
    fullname=RED/${f%_2.fits}_$1.fits
    if [ -f ${fullname%.fits}.solved ]
    then 
	echo ${fullname%.fits}.phot $fullname >> tmpT$1
    fi
done

for mode in D T
do
    if [ $mode == "T" ]
    then
	refname="RED/10-20140226/10-361394"
    else
	refname="RED/10-20140226/10-361395"
    fi

    python ./MagnitudeFitting.py hatnet single $refname"_$1.fits" $refname"_$1.phot" --manual-frame-list tmp$mode$1 --log-config=logging.conf -p 30 --config-file=magfit$1.cfg

    awk '{print $1;}' tmp$mode$1 > tmp$mode$1_photonly
    
    python do_masterphotref.py hatnet $refname"_$1.fits" --log-config=logging.conf tmp$mode$1_photonly --config-file=magfit$1.cfg

    python ./MagnitudeFitting.py hatnet master $refname"_$1.fits" $refname"_$1.phot" --manual-frame-list tmp$mode$1 --log-config=logging.conf -p 30 --config-file=magfit$1.cfg

done