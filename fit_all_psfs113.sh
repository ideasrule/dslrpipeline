#!/bin/bash
for i in 0 1 2 3 
do
    for dir in RED/10-2014*/
    do
	files=$(grep $dir tmp"$i"_photonly)
	python fitpsf_fistar.py $files &
    done    
done


for job in $(jobs -p)
do
    wait $job
done