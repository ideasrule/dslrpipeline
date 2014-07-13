#!/bin/bash
for mode in D T
do
    for i in 0 1 2 3 
    do
	for dir in RED/10-2014*/
	do
	    files=$(grep $dir tmp$mode"$i"_photonly)
	    python fitpsf_fistar.py $files &
	done

    done
done

for job in $(jobs -p)
do
    wait $job
done