#!/bin/bash
#
#script to collect info to populate HDF5 file in create_survey_hdf5.py
#
#Created 2018/01/23
#Author: Erin Mentuch Cooper
#

path_karldb='/work/00115/gebhardt/maverick/gettar'

cat $path_karldb/2017*sci > tmp_DR1
cat $path_karldb/2018*sci >> tmp_DR1

sed -i s/'\/work\/03946\/hetdex\/maverick\/'// tmp_DR1
sed -i s/'\/virus\/'/' '/ tmp_DR1
sed -i s/'virus0000'/' '/ tmp_DR1
sed -i s/'\/'/' '/ tmp_DR1
sed -i s/'\/virus\/'/' '/ tmp_DR1
sed -i s/'_105LL_sci.fits'// tmp_DR1
sed -i s/'_093LL_sci.fits'// tmp_DR1
sed -i s/'Parallel'/'parallel'/ tmp_DR1
grep -v virus000 tmp_DR1 > tmp2_DR1 #this removes test virus files 

awk '{print $1$2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' tmp2_DR1 > mastersci_DR1
sed -i '1s/^/SHOT DATE EXP DITHER TIMESTAMP OBJECT EXPTIME DARKTIME MJD TRAJCRA TRAJCDEC PARANGLE\n/' mastersci_DR1 

rm *tmp*
