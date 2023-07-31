import sys
import os.path as op
import glob
import time
import numpy as np
from astropy.table import Table, vstack

from hetdex_api.detections import Detections
from hetdex_api.shot import Fibers

#awk '{print "python3 /work/05350/ecooper/stampede2/hdr3/get_det_flags.py", $1$2}' /scratch/03946/hetdex/hdr3/survey/hdr3.shotlist > run_get_det_flags

survey = "hdr3"
#survey = "hdr4"

merge = True
if merge:
    flag_cols = ['detectid','flag_pixmask', 'flag_badamp', 'flag_badpix', 'flag_badfib', 'flag_meteor', 'flag_largegal', 'flag_chi2fib' ]
    
    detflagfiles = glob.glob('/scratch/05350/ecooper/det_flags/det*')
    detflagstab = Table(names=flag_cols, dtype=[int, int, int, int, int, int,int, int])

    for file in detflagfiles:
        t = Table.read(file, format='ascii')
        detflagstab = vstack([detflagstab, t])
    detflagstab.write('det_flags_{}.fits'.format(survey))
    sys.exit()

shotid_use = int( sys.argv[1])

if op.exists('det_flags/det_flags_{}.txt'.format(shotid_use)):
    sys.exit("File exists for {}".format(shotid_use))

D = Detections(survey)
detlist = D.hdfile.root.Detections.read_where("shotid == shotid_use")['detectid']
curated_list = np.loadtxt('/work2/05350/ecooper/stampede2/hdr3/catalogs/line_hdr3.0.3.dets', dtype=int)
#curated_list = np.loadtxt('/work2/05350/ecooper/stampede2/hdr4/line_hdr4.0.0.dets', dtype=int)

common_list_det = set( detlist).intersection(curated_list)

flag_cols = ['detectid','flag_pixmask', 'flag_badamp', 'flag_badpix', 'flag_badfib', 'flag_meteor', 'flag_largegal', 'flag_chi2fib' ]
flag_table = Table(names=flag_cols, dtype=[np.int64, np.int32, np.int32, np.int32, np.int32, np.int32, np.int32, np.int32])

t0 = time.time()

if len(common_list_det) > 0:

    F = Fibers(shotid_use, survey=survey)

    for det in common_list_det:

        try:
            r = D.get_detection_flags(det)
            flag_table.add_row([det, r['flag_pixmask'], r['flag_badamp'], r['flag_badpix'], r['flag_badfib'], r['flag_meteor'], r['flag_largegal'], r['flag_chi2fib']])
        except:
            print("Problem getting flags for {}".format(det))


C = Detections(survey, catalog_type='continuum')
detlist = C.hdfile.root.Detections.read_where("shotid == shotid_use")['detectid']
curated_list = np.loadtxt('/work2/05350/ecooper/stampede2/hdr3/catalogs/cont_hdr3.0.3.dets', dtype=int)
#curated_list = np.loadtxt('/work2/05350/ecooper/stampede2/hdr4/cont_hdr4.0.0.dets', dtype=int)

common_list_cont = set( detlist).intersection(curated_list)

if len(common_list_cont) > 0:
    for det in common_list_cont:
        try:
            r = C.get_detection_flags(det)
            flag_table.add_row([det, r['flag_pixmask'], r['flag_badamp'], r['flag_badpix'], r['flag_badfib'], r['flag_meteor'], r['flag_largegal'], r['flag_chi2fib']])
        except:
            print("Problem getting flags for {}".format(det))
            
t1 = time.time()

print("Time to run for {} line and {} continuum detections: {:4.2f} min".format( len(common_list_det), len(common_list_cont), (t1-t0)/60) )

if len(flag_table) > 0:
    flag_table.write('det_flags/det_flags_{}.txt'.format(shotid_use), format='ascii')
