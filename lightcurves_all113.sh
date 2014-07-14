#!/bin/bash
for i in 0 1 2 3
do
    files=$(sed 's/.phot/.aphot/g' tmp"$i"_photonly)
    grcollect -b lightcurves/D$i"_%b" --col-base 1 -m 2g $files &
done


for job in $(jobs -p)
do
    echo $job
    wait $job
done