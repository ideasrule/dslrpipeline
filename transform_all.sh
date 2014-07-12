#!/bin/bash

for dir in $1
do
    echo $dir
    files=$dir*_[0123].fits
    for f in $files
    do
	if [ ! -f ${f%.fits}.solved ]; then continue; fi
	tmpfile=${f%.fits}.tmp
	obj=$(fiheader --get object $f)
	if [ $obj == "G09363730_180" ]
	then
	    wcs-rd2xy -w ${f%.fits}.wcs -i catalog_180.fits -o $tmpfile
	    python read_corr.py $tmpfile | paste 2masscat_180_starlist - | awk '{if ($3 > 0 && $3 < 2601 && $4 > 0 && $4 < 1731) print;}' > ${f%.fits}.catalog
	else
	    wcs-rd2xy -w ${f%.fits}.wcs -i catalog_113.fits -o $tmpfile
	    python read_corr.py $tmpfile | paste 2masscat_113_starlist - | awk '{if ($3 > 0 && $3 < 2601 && $4 > 0 && $4 < 1731) print;}' > ${f%.fits}.catalog
	fi
	rm $tmpfile
    done
done