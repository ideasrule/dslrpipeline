#!/bin/bash
~jhartman/SVN/HATpipe/bin/vartools -l lclist -inputlcformat \
starid:1,x:2,y:3,mag0:49,s:56,t:59,Z:60,h:61 \
-inlistvars avgmag0:2 \
-changevariable mag mag0 \
-expr 'err=1' -rescalesig -rms -clip 5 0 -clip 5 0 -clip 5 0 \
-rms -linfit \
'ax*sin(2*pi*x)+bx*cos(2*pi*x)+ay*sin(2*pi*y)+by*cos(2*pi*y) + ay2*(y-floor(y)) + as*s + az*Z + ah*h + const' ax,bx,ay,by,ay2,as,az,ah,const correctlc \
-expr 'mag0=mag0+avgmag0' -rms -header \
-numbercolumns -o epd nameformat "%s.epdlc" columnformat \
starid,x,y,s,t,Z,h,mag0,err -redirectstats epd/epd_stats -parallel 40