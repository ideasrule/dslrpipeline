#!/bin/bash

for f in 10-2014*/
do
    mkdir "RED/"$f
    flats=$(grep flat database_file | grep $f | awk '{print $1;}')
    python do_3d_masterflat.py "RED/"$f $flats &
done

for job in $(jobs -p)
do
    wait $job
done