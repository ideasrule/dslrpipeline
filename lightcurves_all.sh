#!/bin/bash
for mode in D T
do
    for i in 0 1 2 3
    do
	files=$(sed 's/.phot/.aphot/g' tmp$mode"$i"_photonly)
	grcollect -b lightcurves/$mode$i"_%b" --col-base 1 -m 2g $files &
    done
done

for job in $(jobs -p)
do
    echo $job
    wait $job
done