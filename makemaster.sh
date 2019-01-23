#!/bin/bash
#
#script to collect info to populate HDF5 file in create_survey_hdf5.py
#
#Created 2018/01/23
#Author: Erin Mentuch Cooper
#

path_gpinfo='/work/00115/gebhardt/maverick/getgp/fwhm_and_fluxes_better.txt'
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

cp $path_gpinfo fwhm.tmp
sed -i s/'v'/' '/ fwhm.tmp
awk 'NR>1 {print $1$2,$1,$2,$5}' fwhm.tmp > master_fwhm_DR1
#awk '{print "awk '\''{if($1=="$1"&&$2=="$2") print $0}'\'' fwhm.tmp2 "}' mastersci_DR1 >tmp
#chmod a+x tmp
#tmp > fwhm_master
sed -i '1s/^/SHOT DATE EXP FWHM\n/' master_fwhm_DR1

grep 4540 /work/00115/gebhardt/maverick/detect/2*/res/2*sedtp.dat > tp.tmp
sed -i s/'\/work\/00115\/gebhardt\/maverick\/detect\/'// tp.tmp
sed -i s/'v'/' '/ tp.tmp
sed -i s/'\/res\/'/' '/  tp.tmp
sed -i s/':'/' '/ tp.tmp
sed -i '1s/^/DATE OBSID FILEID WAVE TP TP2 TP3\n/' tp.tmp

cp tp.tmp master_tp

rm *tmp*
