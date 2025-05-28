import sys

import tables as tb
import numpy as np

from astropy.table import Table, join, hstack, unique
import astropy.units as u
from astropy.coordinates import SkyCoord

from hetdex_api.survey import Survey, FiberIndex
from hetdex_api.config import HDRconfig

from hetdex_api.shot import Fibers
import healpy as hp

import time

"""
Load and process fiber data from the HETDEX survey, generating flags and saving them to an HDF5 file.

This script reads survey information and a mask version as command-line arguments, loads the fiber data, 
applies various flags, and saves the results in an HDF5 format.

    sys.argv[1] (str): The survey name. This should be a string representing the specific HETDEX survey 
                       for which you want to load and process fiber data (e.g., 'hdr2.1', 'hdr5').

    sys.argv[2] (str): The mask version. This should be a string identifier for the mask version 
                       being applied to the fiber data (e.g., '4.0.1'). It is used to label the output 
                       HDF5 file and distinguish between different versions of mask processing.

Notes: best to run as a slurm job for about 3 hours (based on hdr4)


Example:

python create_fiber_index_mask.h5 hdr4 4.0.1

"""

survey = sys.argv[1]
mask_version = sys.argv[2]

config = HDRconfig(survey)

make_ascii_tables = True # flag in code whether to make ascii tables or not

t0 = time.time()
FibIndex = FiberIndex(
    load_fiber_table=True, survey=survey, keep_mask=False
)  # don't load previous mask
t1 = time.time()
print("Loaded fibers table in {:3.2f} min.".format((t1 - t0) / 60))

sat_flag = FibIndex.get_satellite_flag()
shot_flag = FibIndex.get_shot_flag()
throughput_flag = FibIndex.get_throughput_flag()
amp_flag = FibIndex.get_amp_flag()
gal_flag = FibIndex.get_gal_flag()
meteor_flag = FibIndex.get_meteor_flag()
badfib_flag = FibIndex.get_badfiber_flag()

FibIndex.fiber_table["flag_shot"] = shot_flag
FibIndex.fiber_table["flag_throughput"] = throughput_flag
FibIndex.fiber_table["flag_badamp"] = amp_flag
FibIndex.fiber_table["flag_largegal"] = gal_flag
FibIndex.fiber_table["flag_meteor"] = meteor_flag
FibIndex.fiber_table["flag_satellite"] = sat_flag
FibIndex.fiber_table["flag_badfib"] = badfib_flag

FibIndex.fiber_table["flag"] = (
    amp_flag
    * gal_flag
    * meteor_flag
    * shot_flag
    * throughput_flag
    * badfib_flag
    * sat_flag
)

flag_table = FibIndex.fiber_table[
    "fiber_id",
    "flag",
    "flag_badamp",
    "flag_badfib",
    "flag_meteor",
    "flag_satellite",
    "flag_largegal",
    "flag_shot",
    "flag_throughput",
]

print(
    "Fraction of flagged fibers is: {:5.3f}".format(
        np.sum(flag_table["flag"]) / len(flag_table)
    )
)

# save bool table to H5
print('Saving to output file.')
fileh = tb.open_file('fiber_mask_{}.h5'.format(mask_version), 'w')
tableflags = fileh.create_table(fileh.root, 'Flags', flag_table.as_array() )
fileh.close()

print('Done')
