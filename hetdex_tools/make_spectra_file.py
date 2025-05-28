import numpy as np
import os.path as op
import os
import sys

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, unique, vstack, join
from astropy.io import fits

from hetdex_api.detections import Detections

shotid = int(sys.argv[1])
shotid_use = shotid

print("Working on {}".format(shotid))
version = "4.0.0"

if op.exists("/scratch/05350/ecooper/spec_tables/spectra_{}_{}.fits".format(version, shotid)):
    sys.exit()

line_dets = np.loadtxt('/scratch/projects/hetdex/hdr4/catalogs/line_hdr{}.dets'.format(version), dtype=int)
cont_dets = np.loadtxt('/scratch/projects/hetdex/hdr4/catalogs/cont_hdr{}.dets'.format(version), dtype=int)

D = Detections(survey="hdr4", loadtable=False)
C = Detections(survey="hdr4", loadtable=False, catalog_type='continuum')

detlist = D.hdfile.root.Detections.read_where("shotid == shotid_use")['detectid']
common_line_det = list( set( detlist).intersection(line_dets))

contdetlist = C.hdfile.root.Detections.read_where("shotid == shotid_use")['detectid']
common_cont_det = list( set( contdetlist).intersection(cont_dets))

def get_spec_for_detlist(det_list, DetectsObject=None):

    ra = []
    dec = []
    spec1d_obs = []
    spec1d_obs_err = []
    spec1d = []
    spec1d_err = []
    apcor = []
    gmag = []

    for det in det_list:
        det_info = DetectsObject.get_detection_info(det)
 
        ra.append(det_info['ra'][0])
        dec.append(det_info['dec'][0])
        spec_tab = DetectsObject.get_spectrum(det, add_apcor=True)
        
        spec1d_obs.append(spec_tab["spec1d"])
        spec1d_obs_err.append(spec_tab["spec1d_err"])
        apcor.append(spec_tab["apcor"])

        spec_tab = DetectsObject.get_spectrum(det, add_apcor=True, deredden=True)
        spec1d.append(spec_tab["spec1d"])
        spec1d_err.append(spec_tab["spec1d_err"])
        gmag.append(DetectsObject.get_hetdex_mag(det))
        
    source_spectra = Table(
        [det_list, ra, dec, gmag, spec1d, spec1d_err, spec1d_obs, spec1d_obs_err, apcor], 
        names=[
            "detectid",
            "ra",
            "dec",
            "gmag",
            "spec",
            "spec_err",
            "spec_obs",
            "spec_obs_err",
            "apcor",
        ],
    )

    intensity_unit = (10 ** -17) * u.erg / u.s / u.cm ** 2 / u.AA

    source_spectra["spec_obs"].unit = intensity_unit
    source_spectra["spec_obs_err"].unit = intensity_unit
    source_spectra["spec"].unit = intensity_unit
    source_spectra["spec_err"].unit = intensity_unit
    source_spectra["ra"].unit = u.deg
    source_spectra["dec"].unit = u.deg

    return source_spectra

# run function
line_spectra = get_spec_for_detlist(common_line_det, DetectsObject=D)
cont_spectra = get_spec_for_detlist(common_cont_det, DetectsObject=C)

source_spectra = vstack( [line_spectra, cont_spectra]) 
source_spectra.write(
    "/scratch/05350/ecooper/spec_tables/spectra_{}_{}.fits".format(version, shotid), overwrite=True
)
