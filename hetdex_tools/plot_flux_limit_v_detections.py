"""

Add the flux limit at the position of the detection
to the detections catalogue

AUTHOR: Daniel Farrow (MPE) 2020

"""

from __future__ import print_function
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from os.path import basename
from argparse import ArgumentParser
import tables as tb
from astropy.table import vstack, Table
from hetdex_api.detections import Detections
from hetdex_api.survey import Survey
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import SensitivityCubeHDF5Container

parser = ArgumentParser(description="Produce test plots of the flux limit versus the flux error")
parser.add_argument("sensitivity_cubes", nargs="+", help="Sensitivity cube(s) HDFs to use")
parser.add_argument("--output_fn", help="Output file to store catalogue with flux limit", default='detect_with_flim.fits')
parser.add_argument("--input_fn",  help="Input file with flux limits to plot", default=None)
parser.add_argument("--sncut", type=float, default=4.5)

args = parser.parse_args()

sncut = args.sncut
detections = Detections("hdr2.1")
survey = Survey("hdr2.1")

tables_out = []

if not args.input_fn:
    for fn in args.sensitivity_cubes:
    
        datevshot = basename(fn)[:12]
        dateshot = int(datevshot.replace("v", ""))
    
        print("Considering {:s}. Parsed dateshot to be {:d}".format(fn, dateshot))
    
        sscube = SensitivityCubeHDF5Container(fn)  
        survey_ttable = survey.hdfile.root.Survey.read_where('datevobs == datevshot')
    
        for ifuslot, scube in sscube.itercubes():
     
            ifuslot_ = ifuslot.replace("ifuslot_", "")
            det_ttable = detections.hdfile.root.Detections.read_where('(shotid == dateshot) & (sn > sncut) & (wave > 3530) & (wave < 5450) & (ifuslot == ifuslot_)')
            det_table = Table(det_ttable)
    
            flim = scube.get_f50(det_ttable["ra"], det_ttable["dec"], det_ttable["wave"], args.sncut)
            det_table["flux_limit_true"] = flim 
    
            tables_out.append(det_table)
    
    table_out = vstack(tables_out)

else:
    table_out = Table.read(args.input_fn)

plt.figure()

cm = plt.cm.get_cmap('RdYlBu')

sc = plt.scatter(table_out["flux"], table_out["flux_limit_true"]*1e17, marker=".",
                 c=table_out["sn"], cmap=cm, vmin=4.5, vmax=9.0)

plt.colorbar(sc, label="Catalogue S/N")

plt.plot([0.0, 1000], [0.0, 1000], "r-")
plt.xlim(0, 50)
plt.ylim(0, 50)
plt.xlabel("Flux $10^{-17}$ erg/s/cm$^{-2}$")
plt.ylabel("50% Flux Limit $10^{-17}$ erg/s/cm$^{-2}$ ")

plt.savefig("flim_test.pdf")

if args.output_fn:
    table_out.write(args.output_fn)
