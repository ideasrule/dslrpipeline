#!/bin/bash

find epd -name '*.epdlc' -exec sort -g -k 5,5 '{}' -o '{}' \;

for mode in D T
do
    for i in 0 1 2 3
    do
	grep $mode"$i"_HAT lclist_tfa_master > lclist_tfa_$mode$i
	sed 's/D0/'$mode''$i'/g' trendlist_master > trendlist$mode$i

	~jhartman/SVN/HATpipe/bin/vartools -l lclist_tfa_$mode$i -inputlcformat starid:1,t:5,mag:8,err:9 -TFA trendlist$mode$i "readformat" 0 5 8 datesfile 25.0 1 0 0 -o tfa columnformat starid,t,mag,err -redirectstats tfa/tfa_stats$mode$i -parallel 40

    done
done