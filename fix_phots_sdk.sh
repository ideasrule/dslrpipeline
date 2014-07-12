#!/bin/bash

for f in "$@"
do
    echo $f
    
    JD=$(fiheader --get jd ${f%.phot}.fits)
    Z=$(fiheader --get z ${f%.phot}.fits)    
    HA=$(fiheader --get ha ${f%.phot}.fits)
    python sdktrans.py ${f%.phot}.psf $f | paste $f - | awk -v jd=$JD -v z=$Z -v ha=$HA -v sdk=$sdk '{if (NR != 1) print $0, jd, z,ha; else print $0,"JD Z HA"}' > ${f%.phot}.aphot

done