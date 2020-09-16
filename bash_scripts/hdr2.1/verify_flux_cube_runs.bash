#!/bin/bash

#
# Verify that the flux limit cubes
# ran and prodcued output of non-zero
# size


for run in $(seq 1 17)
#for run in 17
do
   pushd run$run
   cat rflim | while read line
   do
     ifu=`echo $line | awk '{print $7}'` 
     datevshot=`echo $line | awk '{print $8}'`
     fn="${datevshot}_${ifu}_flim.fits"

     # Find file size and save return value
     duo=`du $fn`
     rval=$?
     size=`echo $duo | awk '{print $1}'` 

     if [ "$rval" != "0" ] || [ "$size" == "0" ]; then
         echo "pushd run${run} ; ${line} ; popd" >> ../torerun3
         echo $fn >> ../failures3
     fi
   done
   popd
done
