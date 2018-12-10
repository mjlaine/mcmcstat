#!/bin/sh

for f in ../stat/*.m; do
     ff=`basename $f`
#     n=`basename $f .eps`.png
     if [ ! -f $ff ];  then
	 echo $ff in stat, not here
     fi
done

# exit

for f in ../mcmc/*.m; do
     ff=`basename $f`
     if [ ! -f $ff ];  then
	 echo $ff in mcmc, not here
     fi
done
