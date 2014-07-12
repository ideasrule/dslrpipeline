#!/bin/bash

#files=10-201405*/*.cr2
#for f in $files
#do 
#    ccdtemp=$(~/Image-ExifTool-9.62/exiftool -CameraTemperature $f | awk '{print $4;}')
#    fiheader --set ccdtemp=$ccdtemp ${f%.cr2}.fits
#    echo $f
#done

#first, produce database file
#./produce_datfile.sh > database_file

#master bias
#rm masterbias_3d*.fits
#biases=$(grep bias database_file | awk '{print $1}')
#python do_3d_masterbias.py $biases

#master darks
#rm masterdark_3d*.fits
#darks=$(grep dark database_file | grep 180 | awk '{print $1}')
#python do_3d_masterdark.py masterbias_3d.fits $darks

#master flats for every day
#./do_3d_masterflats.sh

#reduce frames and put them into RED/
#./calibrate_all.sh

#multiprocessing (40 cores) astrometry
python run_astrometry.py RED/10-2014*/*_[0123].fits

./multip.py ./transform_all.sh RED/10-2014*/

./multip.py ./photometry_all.sh RED/10-2014*/

#don't forget to set FILTERS for template frames!
fiheader --set FILTERS=V RED/10-20140226/10-36139[45]_[03].fits
fiheader --set FILTERS=R RED/10-20140226/10-36139[45]_2.fits
fiheader --set FILTERS=B RED/10-20140226/10-36139[45]_1.fits

./multip.py ./magfit_all.sh 0 1 2 3

./fit_all_psfs.sh

./multip.py ./fix_phots_sdk.sh RED/10-2014*/*.phot

./lightcurves_all.sh

cd lightcurves
#remove all nan's, remove files with less than 1000 points
./cleanup.sh

./create_lclist.sh
./run_epd.sh
#create datesfile
awk '{print $1,$7}' ../database_file > datesfile
./run_tfa.sh

cd tfa
for f in D0_HAT*; do filename=${f:3:15}; python print_devs.py $filename D; done > tfa_stats
