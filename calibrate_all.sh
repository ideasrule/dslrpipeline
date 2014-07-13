#!/bin/bash

#get first flat
for dir in RED/10-2014*/
do
    if [ -f $dir"masterflat0.fits" ]
    then
	masterflatdir=$dir
	break
    fi
done

for dir in 10-2014*/
do
    if [ -f "RED/"$dir"masterflat0.fits" ]
    then
	masterflatdir="RED/"$dir
    fi
    files=$(grep object database_file | grep $dir | awk '{if($4==180 && $5==0)print $1}')
    python -u calibrate.py $masterflatdir $files &
done

for job in $(jobs -p)
do
    wait $job
done