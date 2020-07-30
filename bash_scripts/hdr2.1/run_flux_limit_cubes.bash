#!/bin/bash

#
#
#  Split the flux limit cube run
#  file into sub files and sub-directories
#  to make running it easier
#
# example rflim0 file
# /work/00115/gebhardt/maverick/scripts/rsp/rfitsp1 150.043976 1.96127379 35 4505 50 004_106_033 20170107v013 1.7 3 3.5 2 1 105

# number of shots per folder
NPERFOLDER=200

# do the split or just the job submission
splitup=false
submit=true
run=0
i=$NPERFOLDER
saved_datevobs=none
pushd .

if [ "$splitup" == "true" ]
then
   cat /data/00115/gebhardt/getcen/rflim0 | while read line
    do
        datevobs=`echo $line | awk '{print $8}'`
        if [ "$datevobs" != "$saved_datevobs" ] 
        then
            saved_datevobs=$datevobs
            if [ $i -ge  $NPERFOLDER ] 
            then
                i=0
                run=$(( run + 1 ))
                popd
                mkdir run$run 
                pushd run$run
            fi
            i=$((i + 1))
            echo $datevobs >> shotlist
        fi
        echo $line >> rflim
    done
fi

if [ "$submit" == "true" ] 
then
   for run in $(seq 10 17)
   do
       pushd run$run
       nlines=`wc rflim | awk '{print $1}'`
   
       # number per line for 20 tasks
       nperline=`python -c "print(int(${nlines}/19))"`
       echo $nperline
       ~gebhardt/jobsplitter2w rflim 20 $nperline 24:00:00
       #sbatch rflim_1.slurm 
       popd
   done
fi



