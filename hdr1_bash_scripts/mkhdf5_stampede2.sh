# bash script to run the HDF5 shot creation scripts on stampede2 
workdir=${pwd}
cd /scratch/03946/hetdex/
python /work/03946/hetdex/hdr1/software/HETDEX_API/create_shot_hdf5.py -d $1 -o $2 -of $1v$2.h5
python /work/03946/hetdex/hdr1/software/HETDEX_API/create_astrometry_hdf5.py -d $1 -o $2 -of $1v$2.h5 --append
python /work/03946/hetdex/hdr1/software/HETDEX_API/create_cal_hdf5.py -d $1 -o $2 -of $1v$2.h5 --append
mv $1v$2.h5 $workdir
