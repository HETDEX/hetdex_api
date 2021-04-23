# python3 make_curated_catalog.py 2.1.2
#
import sys
import os.path as op
import tables as tb
from astropy.table import Table, join, Column
import numpy as np

from multiprocessing import Pool

from hetdex_api.detections import Detections
from hetdex_api.config import HDRconfig
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import SensitivityCubeHDF5Container


def return_fiber_ratio(det):
    fiber_row = fiber_table.read_where("detectid == det")
    weights = np.sort(fiber_row["weight"])
    fiber_ratio = weights[-1] / weights[-2]
    return fiber_ratio


def get_flux_noise_1sigma(detid, mask=False):

    """ For a detid get the the f50_1sigma value

    No need to mask the flux limits if you are using the
    curated catalog which is already masked. This avoid a
    few good sources on the edges of masks that get flagged
    with a high flux limit value.
    """

    global config, det_table

    sncut = 1

    sel_det = det_table["detectid"] == detid
    shotid = det_table["shotid"][sel_det][0]
    ifuslot = det_table["ifuslot"][sel_det][0]

    det_table_here = det_table[sel_det]

    datevobs = str(shotid)[0:8] + "v" + str(shotid)[8:11]

    fn = op.join(config.flim_dir, datevobs + "_sensitivity_cube.h5")
    if mask:
        mask_fn = op.join(config.flimmask, datevobs + "_mask.h5")
        sscube = SensitivityCubeHDF5Container(fn, mask_filename=mask_fn, flim_model="hdr1")
    else:
        sscube = SensitivityCubeHDF5Container(fn, flim_model="hdr1")
    scube = sscube.extract_ifu_sensitivity_cube("ifuslot_{}".format(ifuslot))
    flim = scube.get_f50(
        det_table_here["ra"], det_table_here["dec"], det_table_here["wave"], sncut
    )
    sscube.close()

    return flim[0]


version = str(sys.argv[1])

config = HDRconfig()

# Note because refine is constantly updated, it isn't possible to
# truly replicate older catalogs. TODO for HDR3

detects = Detections(survey="hdr2.1", catalog_type="lines").refine()

sel_field = (
    (detects.field == "cosmos")
    | (detects.field == "dex-fall")
    | (detects.field == "dex-spring")
    | (detects.field == "egs")
    | (detects.field == "goods-n")
)

if version == "2.1.1":
    sel_chi2 = detects.chi2 < 1.2
    sel_wave = (detects.wave >= 3510) * (detects.wave <= 5490)
    sel_lw = detects.linewidth <= 6
    sel_cont = detects.continuum > -3
    sel_sn = detects.sn >= 4.8
    sel_chi2fib = detects.chi2fib < 4.5

    sel_cat = sel_field * sel_chi2 * sel_wave * sel_lw * sel_cont * sel_sn * sel_chi2fib

elif version == "2.1.2":

    sel_cut1 = (detects.sn >= 7) * (detects.chi2 <= 2.5)
    sel_cut2 = (detects.sn >= 4.8) * (detects.sn < 7) * (detects.chi2 <= 1.2)

    sel_cont = detects.continuum > -3
    sel_chi2fib = detects.chi2fib < 4.5
    sel_tp = detects.throughput >= 0.08

    sel = sel_field * sel_cont * sel_chi2fib * sel_tp * (sel_cut1 | sel_cut2)

    sel_wave = (detects.wave >= 3550) * (detects.wave <= 5470)
    sel_lw = (detects.linewidth <= 14) * (detects.linewidth > 6)

    sel1 = sel * sel_wave * sel_lw

    sel_wave = (detects.wave >= 3510) * (detects.wave <= 5490)
    sel_lw = detects.linewidth <= 6

    sel2 = sel * sel_wave * sel_lw

    sel_cat = sel1 | sel2

    det_table = detects[sel_cat].return_astropy_table()

    elixer_file = op.join(config.detect_dir, "catalogs", "elixer.2.1.2.h5")
    elixer_cat = tb.open_file(elixer_file, "r")

    cls = []
    mlname = []
    mlz = []
    mlprob = []

    for row in det_table:
        detectid_obj = row["detectid"]
        try:
            elix_row = elixer_cat.root.Detections.read_where("detectid == detectid_obj")
            row["plae_classification"] = elix_row["plae_classification"]
            row["combined_plae"] = elix_row["combined_plae"]
            row["combined_plae_err"] = elix_row["combined_plae_err"]
            mlname.append(elix_row["multiline_name"][0].decode())
            cls.append(elix_row["classification_labels"][0].decode())
            mlz.append(elix_row["multiline_z"][0])
            mlprob.append(elix_row["multiline_prob"][0])
        except Exception:
            mlname.append("")
            cls.append("")
            mlz.append(False)
            mlprob.append(0.0)

    det_table.add_column(mlname, name="multiline_name")
    det_table.add_column(cls, name="classification_labels")

elif version == "2.1.3":

    # get updated chi2fib values (see work in stampede2/notebooks/fiber_chi2_check.ipynb:
    chi2table = Table.read('chi2fib_all.tab', format='ascii.no_header', names=['detectid', 'chi2fib'])

    det = Table([detects.detectid],names=['detectid'])
    det_join = join(det, chi2table)

    detectid2 = np.array(det_join['detectid'])
    detects.chi2fib = np.array(det_join['chi2fib'])

    if np.sum(detects.detectid - detectid2) != 0:
        print('Something went wrong with appending updated chi2fib')

    sel_cut1 = np.logical_not(detects.gmag < 19) * (detects.sn >= 7) * (detects.chi2 <= 2.5) *  (detects.chi2 > 1.2) * (detects.continuum < 7.5)
    sel_cut2 = (detects.sn >= 4.8) * (detects.chi2 <= 1.2)

    sel_cont = detects.continuum > -3

    selchi2fib1 = detects.chi2fib < 4.5
    selchi2fib2 = np.invert((detects.chi2fib > 3) * (detects.continuum < 0.5))
    sel_chi2fib = selchi2fib1 * selchi2fib2

    sel_tp = detects.throughput >= 0.08

    sel = sel_field * sel_cont * sel_chi2fib * sel_tp * (sel_cut1 | sel_cut2)

    sel_wave = (detects.wave >= 3550) * (detects.wave <= 5470)
    sel_lw = (detects.linewidth <= 14) * (detects.linewidth > 6)

    sel1 = sel * sel_wave * sel_lw

    sel_wave = (detects.wave >= 3510) * (detects.wave <= 5490)
    sel_lw = detects.linewidth <= 6

    sel2 = sel * sel_wave * sel_lw

    sel_cat = sel1 | sel2

    det_table = detects[sel_cat].return_astropy_table()

    elixer_file = op.join(config.detect_dir, "catalogs", "elixer.2.1.2.h5")
    elixer_cat = tb.open_file(elixer_file, "r")

    cls = []
    mlname = []
    mlz = []
    mlprob = []

    counterpart_mag = []
    counterpart_mag_err = []
    counterpart_dist = []
    counterpart_catalog_name = []
    counterpart_filter_name = []

    fixed_mag = []
    fixed_mag_err = []
    fixed_catalog_name = []
    fixed_filter_name = []
    fixed_radius = []

    for row in det_table:
        detectid_obj = row["detectid"]
        try:
            elix_row = elixer_cat.root.Detections.read_where("detectid == detectid_obj")
            row["plae_classification"] = elix_row["plae_classification"]
            row["combined_plae"] = elix_row["combined_plae"]
            row["combined_plae_err"] = elix_row["combined_plae_err"]
            mlname.append(elix_row["multiline_name"][0].decode())
            cls.append(elix_row["classification_labels"][0].decode())
            mlz.append(elix_row["multiline_z"][0])
            mlprob.append(elix_row["multiline_prob"][0])
        except Exception:
            mlname.append("")
            cls.append("")
            mlz.append(False)
            mlprob.append(0.0)

        # append nearest source extracted neighbour match
        try:

            elix_row = elixer_cat.root.ExtractedObjects.read_where(
                "(detectid == detectid_obj) & (selected == True)"
            )

            if np.size(elix_row) == 0:
                elix_row = elixer_cat.root.ExtractedObjects.read_where(
                    "detectid == detectid_obj"
                )
            if np.size(elix_row) > 1:
                sel_r = elix_row["filter_name"] == b"r"
                if np.sum(sel_r) == 1:
                    counterpart_mag.append(elix_row["mag"][sel_r][0])
                    counterpart_mag_err.append(elix_row["mag_err"][sel_r][0])
                    counterpart_dist.append(elix_row["dist_baryctr"][sel_r][0])
                    counterpart_catalog_name.append(
                        elix_row["catalog_name"][sel_r][0].decode()
                    )
                    counterpart_filter_name.append(
                        elix_row["filter_name"][sel_r][0].decode()
                    )
                else:
                    counterpart_mag.append(elix_row["mag"][0])
                    counterpart_mag_err.append(elix_row["mag_err"][0])
                    counterpart_dist.append(elix_row["dist_baryctr"][0])
                    counterpart_catalog_name.append(
                        elix_row["catalog_name"][0].decode()
                    )
                    counterpart_filter_name.append(elix_row["filter_name"][0].decode())
            elif np.size(elix_row) == 1:
                counterpart_mag.append(elix_row["mag"][0])
                counterpart_mag_err.append(elix_row["mag_err"][0])
                counterpart_dist.append(elix_row["dist_baryctr"][0])
                counterpart_catalog_name.append(elix_row["catalog_name"][0].decode())
                counterpart_filter_name.append(elix_row["filter_name"][0].decode())
            else:
                counterpart_mag.append(np.nan)
                counterpart_mag_err.append(np.nan)
                counterpart_dist.append(np.nan)
                counterpart_catalog_name.append("")
                counterpart_filter_name.append("")
        except:
            counterpart_mag.append(np.nan)
            counterpart_mag_err.append(np.nan)
            counterpart_dist.append(np.nan)
            counterpart_catalog_name.append("")
            counterpart_filter_name.append("")

        # append fixed aperture mag
        try:
            elix_tab = elixer_cat.root.ElixerApertures.read_where(
                ("detectid == detectid_obj")
            )
            sel_r = elix_tab["filter_name"] == b"r"
            sel_g = elix_tab["filter_name"] == b"g"

            if np.any(sel_r):
                elix_r = elix_tab[sel_r]
                fixed_mag.append(elix_r["mag"][-1])
                fixed_mag_err.append(elix_r["mag_err"][-1])
                fixed_catalog_name.append(elix_r["catalog_name"][-1].decode())
                fixed_filter_name.append(elix_r["filter_name"][-1].decode())
                fixed_radius.append(elix_r["radius"][-1])
            elif np.any(sel_g):
                elix_g = elix_tab[sel_g]
                fixed_mag.append(elix_g["mag"][-1])
                fixed_mag_err.append(elix_g["mag_err"][-1])
                fixed_catalog_name.append(elix_g["catalog_name"][-1].decode())
                fixed_filter_name.append(elix_g["filter_name"][-1].decode())
                fixed_radius.append(elix_g["radius"][-1])
            else:
                sel = elix_tab["radius"] < 3
                elix_sel = elix_tab[sel]
                fixed_mag.append(elix_sel["mag"][-1])
                fixed_mag_err.append(elix_sel["mag_err"][-1])
                fixed_catalog_name.append(elix_sel["catalog_name"][-1].decode())
                fixed_filter_name.append(elix_sel["filter_name"][-1].decode())
                fixed_radius.append(elix_sel["radius"][-1])
        except:
            fixed_mag.append(np.nan)
            fixed_mag_err.append(np.nan)
            fixed_catalog_name.append("")
            fixed_filter_name.append("")
            fixed_radius.append(np.nan)

    det_table.add_column(mlname, name="multiline_name")
    det_table.add_column(cls, name="classification_labels")
    det_table.add_column(counterpart_mag, name="counterpart_mag")
    det_table.add_column(counterpart_mag_err, name="counterpart_mag_err")
    det_table.add_column(counterpart_dist, name="counterpart_dist")
    det_table.add_column(counterpart_catalog_name, name="counterpart_catalog_name")
    det_table.add_column(counterpart_filter_name, name="counterpart_filter_name")
    det_table.add_column(fixed_mag, name="forced_mag")
    det_table.add_column(fixed_mag_err, name="forced_mag_err")
    det_table.add_column(fixed_catalog_name, name="forced_catalog_name")
    det_table.add_column(fixed_filter_name, name="forced_filter_name")
    det_table.add_column(fixed_radius, name="forced_radius")
else:
    print("Provide a version : eg. 2.1.2")
    sys.exit()

fiber_table = detects.hdfile.root.Fibers

fiber_ratio = []
for det in det_table["detectid"]:
    fiber_row = fiber_table.read_where("detectid == det")
    weights = np.sort(fiber_row["weight"])
    fiber_ratio.append(weights[-1] / weights[-2])

det_table.add_column(fiber_ratio, name="fiber_ratio")

# get f50 values from Donghui:
#print('Adding f50 values from Donghui')
#fileh = tb.open_file("hdr213_fullfield_source_info_210323.h5", "r")
#detectid = fileh.root.hdr213_fullfield.detectid[:]
#f50 = fileh.root.hdr213_fullfield.f50[:]
#f50_tab = Table(
#    [detectid, f50[0], f50[1], f50[2], f50[3], f50[4], f50[5]],
#    names=["detectid", "f50_4pt8", "f50_5", "f50_5pt5", "f50_6", "f50_6pt5", "f50_7"],
#)
#fileh.close()

# add 1sigma flim value for each detection from HDR1 flux limit model
print('Adding f50_1sigma from flux limits')
p = Pool(24)
flim = p.map(get_flux_noise_1sigma, det_table['detectid'])
p.close()

det_table['flux_noise_1sigma'] = flim

det_table_join = join(det_table, f50_tab, join_type="left")

det_table_join.write("detect_hdr{}.fits".format(version), overwrite=True)
det_table_join.write("detect_hdr{}.tab".format(version), format="ascii", overwrite=True)
