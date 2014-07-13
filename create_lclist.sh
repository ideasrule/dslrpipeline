#!/bin/bash

rm lclist lclist_tfa_master
for f in *_HAT*
do
    starname=${f:3:200}
    ch=${f:1:1}
    chindex=0
    if [ "$ch" == 0 ] || [ "$ch" == 3 ]
    then 
	chindex=15; 
    elif [ "$ch" == 2 ]
    then
	chindex=16
    else
	chindex=14
    fi
    mag=$(grep $starname ../2masscat_180 | awk -v chindex=$chindex '{print $chindex}')
    echo $f $mag >> lclist
    echo epd/$f.epdlc $(head -1 $f | awk '{print $2,$3}') >> lclist_tfa_master
done