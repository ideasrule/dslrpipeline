#!/bin/bash

files=10-2014*/10-*.fits

for f in $files
do
    object=$(fiheader --get object $f)
    imagetyp=$(fiheader --get imagetyp $f)
    aborted=$(fiheader --get aborted $f)
    exptime=$(fiheader --get exptime $f)
    mntstate=$(fiheader --get mntstate $f)
    jd=$(fiheader --get jd $f)
    echo $f $imagetyp $object $exptime $aborted $mntstate $jd
done