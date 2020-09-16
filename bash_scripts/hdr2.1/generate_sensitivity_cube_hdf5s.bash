#!/bin/bash
#
#
# Generate the HDF5 sensitivity cubes from
# the flux limit FITs files. Also compute
# the averages across a shot
# 
#  Arguments:
#  $1  = run number to consider
#
# Daniel Farrow (MPE)
#

# Script to compute averages
AV_FLIM_SCRIPT=/home/04120/dfarrow/HETDEX_API_devel/hetdex_tools/compute_average_one_sigma.py 

# List of datevshots
DATE_FN=shotlist
redo="false"

# input and output paths
INPATH=/data/04120/dfarrow/hdr2.1_fits_cubes_rerun/run$1
OUTPATH=/data/04120/dfarrow/hdr2.1_hdfs_rerun

# File in to which to write the biweight location of the flux limits at 4540AA for each shot
AVERAGE_FLIM_FILE=$OUTPATH/average_one_sigma_run$1.txt

cd $INPATH

for datevobs in `cat $DATE_FN`
do
 
    echo $datevobs

    if [ -f ${OUTPATH}/${datevobs}_sensitivity_cube.h5 ]; then
        if [ "$redo" == "true" ]
        then
            mv ${OUTPATH}/${datevobs}_sensitivity_cube.h5 ${OUTPATH}/${datevobs}_sensitivity_cube.h5.old
        else
            echo "${OUTPATH}/${datevobs}_sensitivity_cube.h5 exists, skipping"
            continue
        fi
    fi
  
    files=`ls ${datevobs}_???_???_???_flim.fits`
    echo $files
    python $AV_FLIM_SCRIPT --fn-shot-average $AVERAGE_FLIM_FILE $files
    add_sensitivity_cube_to_hdf5 $files "${OUTPATH}/${datevobs}_sensitivity_cube.h5"
done
