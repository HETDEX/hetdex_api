import sys
import os.path as op
import glob
import time
import numpy as np
from astropy.table import Table, vstack

from hetdex_api.detections import Detections
from hetdex_api.shot import Fibers

merge = True

if merge:
    version='4.0.0'
#    version='3.0.3'
    flag_cols = ['detectid','flag_pixmask', 'flag_badamp', 'flag_badpix', 'flag_badfib', 'flag_meteor', 'flag_largegal', 'flag_chi2fib' ]
    
    detflagfiles = glob.glob('/scratch/05350/ecooper/det_flags/det*')
    detflagfiles.sort()
    detflagstab = Table(names=flag_cols, dtype=[int, int, int, int, int, int,int, int])

    for file in detflagfiles:
        print(file)
        t = Table.read(file, format='ascii')
        detflagstab = vstack([detflagstab, t])
    detflagstab.write('det_flags_{}.fits'.format(version))
    sys.exit()

    
shotid_use = int( sys.argv[1])


#if op.exists('det_flags_hdr3/det_flags_{}.txt'.format(shotid_use)):
#    sys.exit("File exists for {}".format(shotid_use))


if False:#for hdr3 run shotid_use < 20210901000:
    survey = "hdr3"
else:
    survey = "hdr4"

D = Detections(survey=survey)

detlist = D.hdfile.root.Detections.read_where("shotid == shotid_use")['detectid']

if survey == 'hdr3':
    curated_list = np.loadtxt('/work2/05350/ecooper/stampede2/hdr3/catalogs/line_hdr3.0.3.dets', dtype=int)
else:
    curated_list = np.loadtxt('/scratch/projects/hetdex/hdr4/catalogs/line_hdr4.0.0.dets', dtype=int)

common_list_det = set( detlist).intersection(curated_list)

flag_cols = ['detectid','flag_pixmask', 'flag_badamp', 'flag_badpix', 'flag_badfib', 'flag_meteor', 'flag_largegal', 'flag_chi2fib' ]
flag_table = Table(names=flag_cols, dtype=[np.int64, np.int32, np.int32, np.int32, np.int32, np.int32, np.int32, np.int32])

t0 = time.time()

F = Fibers(shotid_use, survey=survey)

if len(common_list_det) > 0:

    for det in list( common_list_det):

        try:
            r = D.get_detection_flags(det, F=F)
            flag_table.add_row([det, r['flag_pixmask'], r['flag_badamp'], r['flag_badpix'], r['flag_badfib'], r['flag_meteor'], r['flag_largegal'], r['flag_chi2fib']])
        except:
            print("Problem getting flags for {}".format(det))

    D.close()
            
C = Detections(survey, catalog_type='continuum')
detlist = C.hdfile.root.Detections.read_where("shotid == shotid_use")['detectid']

if survey == 'hdr3':
    curated_list = np.loadtxt('/work2/05350/ecooper/stampede2/hdr3/catalogs/cont_hdr3.0.3.dets', dtype=int)
else:
    curated_list = np.loadtxt('/scratch/projects/hetdex/hdr4/catalogs/cont_hdr4.0.0.dets', dtype=int)

common_list_cont = set( detlist).intersection(curated_list)

if len(common_list_cont) > 0:
    for det in list(common_list_cont):
        try:
            r = C.get_detection_flags(det, F=F)
            flag_table.add_row([det, r['flag_pixmask'], r['flag_badamp'], r['flag_badpix'], r['flag_badfib'], r['flag_meteor'], r['flag_largegal'], r['flag_chi2fib']])
        except:
            print("Problem getting flags for {}".format(det))

F.close()
C.close()            
t1 = time.time()

print("Time to run for {} line and {} continuum detections: {:4.2f} min".format( len(common_list_det), len(common_list_cont), (t1-t0)/60) )

if len(flag_table) > 0:
    flag_table.write('det_flags/det_flags_{}.txt'.format(shotid_use), format='ascii', overwrite=True)
