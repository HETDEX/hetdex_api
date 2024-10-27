import sys
import os.path as op
import glob
import time
import numpy as np
from astropy.table import Table, vstack

from hetdex_api.survey import FiberIndex
from hetdex_api.detections import Detections
from hetdex_api.shot import Fibers

merge = True
version = "4.0.1"

flag_cols = [
    "detectid",
    "flag_pixmask",
    "flag_badamp",
    "flag_badpix",
    "flag_badfib",
    "flag_meteor",
    "flag_largegal",
    "flag_chi2fib",
    "flag_satellite",
    "flag_cal",
    "flag_throughput",
    "flag_shot",
]

if merge:
    detflagfiles = glob.glob("/scratch/05350/ecooper/det_flags/det_flags*")
    detflagfiles.sort()
    detflagstab = Table(names=flag_cols, dtype=[int, int, int, int, int, int, int, int, int, int, int, int])

    for file in detflagfiles:
        print(file)
        t = Table.read(file, format="ascii")
        detflagstab = vstack([detflagstab, t])
    detflagstab.sort('detectid')
    detflagstab.write("det_flags_{}.fits".format(version))
    sys.exit()

shotid_use = int(sys.argv[1])

if shotid_use in np.loadtxt(
    "/scratch/projects/hetdex/hdr4/catalogs/shots_hdr3_{}.txt".format(version),
    dtype=int,
):
    survey = "hdr3"
elif shotid_use in np.loadtxt(
    "/scratch/projects/hetdex/hdr4/catalogs/shots_hdr4_{}.txt".format(version),
    dtype=int,
):
    survey = "hdr4"
else:
    print("Something's wrong. {} not in shot list".format(shotid_use))
    sys.exit()

D = Detections(survey=survey)

detlist = D.hdfile.root.Detections.read_where("shotid == shotid_use")["detectid"]

if survey == "hdr3":
    curated_list = np.loadtxt(
        "/scratch/projects/hetdex/hdr4/catalogs/line_hdr3_{}.dets".format(version),
        dtype=int,
    )
elif survey == "hdr4":
    curated_list = np.loadtxt(
        "/scratch/projects/hetdex/hdr4/catalogs/line_hdr4_{}.dets".format(version),
        dtype=int,
    )
elif survey == "hdr5":
    curated_list = np.loadtxt(
        "/scratch/projects/hetdex/hdr5/catalogs/line_hdr5.0.0.dets", dtype=int
    )
else:
    print(f"Unsupported survey {survey}")

common_list_det = set(detlist).intersection(curated_list)

flag_table = Table(
    names=flag_cols,
    dtype=[
        np.int64,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
    ],
)

t0 = time.time()

FI = FiberIndex()
F = Fibers(shotid_use, survey=survey)

if len(common_list_det) > 0:
    for det in list(common_list_det):
        #print(det)
        try:
            r = D.get_detection_flags(det, F=F, FI=FI)
            flag_table.add_row(
                [
                    det,
                    r["flag_pixmask"],
                    r["flag_badamp"],
                    r["flag_badpix"],
                    r["flag_badfib"],
                    r["flag_meteor"],
                    r["flag_largegal"],
                    r["flag_chi2fib"],
                    r["flag_satellite"],
                    r["flag_badcal"],
                    r["flag_throughput"],
                    r["flag_shot"],
                ]
            )
        except:
            print("Problem getting flags for {}".format(det))

    D.close()

C = Detections(survey, catalog_type="continuum")
detlist = C.hdfile.root.Detections.read_where("shotid == shotid_use")["detectid"]

if survey == "hdr3":
    curated_list = np.loadtxt(
        "/scratch/projects/hetdex/hdr4/catalogs/cont_hdr3_{}.dets".format(version),
        dtype=int,
    )
elif survey == "hdr4":
    curated_list = np.loadtxt(
        "/scratch/projects/hetdex/hdr4/catalogs/cont_hdr4_{}.dets".format(version),
        dtype=int,
    )
elif survey == "hdr5":
    curated_list = np.loadtxt(
        "/scratch/projects/hetdex/hdr5/catalogs/cont_hdr5.0.0.dets", dtype=int
    )
else:
    print(f"Unsupported survey {survey}")

common_list_cont = set(detlist).intersection(curated_list)
print('common_list_cont size ', len(common_list_cont))
if len(common_list_cont) > 0:
    for det in list(common_list_cont):
        try:
            r = C.get_detection_flags(det, F=F, FI=FI)
            flag_table.add_row(
                [
                    det,
                    r["flag_pixmask"],
                    r["flag_badamp"],
                    r["flag_badpix"],
                    r["flag_badfib"],
                    r["flag_meteor"],
                    r["flag_largegal"],
                    r["flag_chi2fib"],
                    r["flag_satellite"],
                    r["flag_badcal"],
                    r["flag_throughput"],
                    r["flag_shot"],
                ]
            )

        except:
            print("Problem getting flags for {}".format(det))
FI.close()
F.close()
C.close()
t1 = time.time()

print(
    "Time to run for {} line and {} continuum detections: {:4.2f} min".format(
        len(common_list_det), len(common_list_cont), (t1 - t0) / 60
    )
)

if len(flag_table) > 0:
    flag_table.write(
        "det_flags_{}.txt".format(shotid_use), format="ascii", overwrite=True
    )
