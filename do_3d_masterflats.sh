#!/bin/bash

for f in 10-2014*/
do
    mkdir "RED/"$f
    flats=$(grep flat database_file | grep $f | awk '{print $1;}')
    mastername="RED/"$f"masterflat.fits"
    python do_3d_masterflat.py masterbias_3d.fits $mastername $flats
done