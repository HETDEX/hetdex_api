#!/bin/bash 
#
# Loop over a list of VIRUS shots
# and submit slurm jobs to produce their data
# cubes. Wait for every NRUN to stop before
# submitting again

INLIST=$1
NRUN=10 # Number to submit at once
counter=1

while read line 
do

   date=`echo ${line} | awk '{print $2}'`
   shot=`echo ${line} | awk '{print $3}'`
   outfolder=`echo ${line} | awk '{print $2"v"$3}'`
   echo $date $shot

   # Run the setup
   vdrp_setup_flim $date $shot

   # cd into this dir
   pushd $outfolder/flim/

   # Find the name of the slurm file
   slurm_name=`ls flim*.slurm`
   counter=`expr ${counter} + 1`

   if [ "$counter" -eq "$NRUN" ]; then
       counter=1
       echo "Waiting for "$slurm_name
       sbatch --wait $slurm_name 
   else
       # Run without waiting
       echo "Run "$slurm_name
       sbatch $slurm_name
   fi

   popd

done<$INLIST
