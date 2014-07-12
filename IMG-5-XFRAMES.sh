#!/bin/bash

Cap=18		# use column 18 for rejection for ap
Cism=28         # use column 28 for rejection for ism
SIG=3		# rejection: 3 sigma
Bap=2		# basename column, with apphot
Bism=1        # basename column, with ism
#B=1;		# With ISM
LCEXTap=rlc;	# LC extension
LCEXTism=ilc;	# With ISM
STARLISTFILEap=BASE/2MASS.txt;		# With apphot
STARLISTFILEism=BASE/photref.cat;			# With ISM

LCEXT=
STARLISTFILE=
C=
B=

TYPE=ap

SELF=$0

star=

while (( "$#" )) ; do
	if [ "$1" = "--ism" ] ; then
		TYPE=ism
	elif [ "$1" = "--ap" ] ; then
		TYPE=ap
	elif [ "$1" = "-c" ] ; then
		shift 1
		C=$1
	elif [ "$1" = "--sig" ] ; then
		shift 1
		SIG=$1
	elif [ "$1" = "-b" ] ; then
		shift 1
		B=$1
	elif [ "$1" = "--lcext" ] ; then
		shift 1
		LCEXT=$1
	elif [ "$1" = "--starlistfile" ] ; then
		shift 1
		STARLISTFILE=$1
	elif [ "$1" = "--help" ] ; then
		cat > /dev/stderr <<EOF
Usage:	$SELF [--help]
	[--ism | --ap] [-c colref] [--sig rejectionsigma]
	[-b basenamecolumn] [--lcext LCextension]
	[--starlistfile starlistfilename] [star1 star2 ...]
EOF
		exit 1
	else
		star=$1
	fi
	shift 1
done

if [ "$TYPE" = "ap" ] ; then
	for vartocheck in LCEXT STARLISTFILE B C ; do
		if [ -z "${!vartocheck}" ] ; then
			vartosub=${vartocheck}ap
			eval ${vartocheck}=${!vartosub}
		fi
	done
elif [ "$TYPE" = "ism" ] ; then
	for vartocheck in LCEXT STARLISTFILE B C ; do
		if [ -z "${!vartocheck}" ] ; then
			vartosub=${vartocheck}ism
			eval ${vartocheck}=${!vartosub}
		fi
	done
fi

if [ -n "$star" ] ; then

	lc=LC/$star.${LCEXT}

	if [ ! -f $lc ] ; then exit 0 ; fi

	# By GB
	gstater -i $lc -F -l $SIG -g -t -a -c $C -y -L - | \
		gawk -v B=$B '!/^#/ {print $B,1}'

else

	if [ ! -f photref.starlist ] ; then
		gawk '{print $1}' $STARLISTFILE | sort > photref.starlist
	fi	

	pexec	-n 4 -f photref.starlist -o - -a '%l' -u - -e star \
		-c "$SELF --lcext $LCEXT -b $B -c $C --sig $SIG \$star" | \
	grcollect - \
		--col-base 1 --col-stat 2 --stat count \
		--output - | \
	sort	-k 2,2n > fscount.dat
fi

