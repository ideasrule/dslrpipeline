#!/bin/bash
for f in $1
do
    numlines=$(grep -v nan $f | wc -l)
    if [ $numlines -lt 1000 ]
    then 
	rm $f
    else
	grep -v nan $f | awk '{if(NF==61)print;}' > $f
    fi
done
