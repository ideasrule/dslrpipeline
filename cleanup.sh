#!/bin/bash
for f in *_HAT*
do
    grep -v nan $f > tmpfile
    numlines=$(wc -l < tmpfile)
    if [ $numlines -lt 1000 ]
    then 
	rm $f
    else
	awk '{if(NF==61)print;}' tmpfile > $f
#	mv tmpfile $f
    fi
done
