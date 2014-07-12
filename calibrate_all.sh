#!/bin/bash

#get first flat
for dir in RED/10-2014*/
do
    if [ -f $dir"masterflat.fits" ]
    then
	masterflat=$dir"masterflat.fits"
	break
    fi
done

for dir in 10-2014*/
do
    if [ -f "RED/"$dir"masterflat.fits" ]
    then
	masterflat="RED/"$dir"masterflat.fits"
    fi
    files=$(grep object database_file | grep $dir | awk '{if($4==180 && $5==0)print $1}')
    echo $dir
    echo $files
    python -u calibrate.py masterbias_3d.fits $masterflat $files
done