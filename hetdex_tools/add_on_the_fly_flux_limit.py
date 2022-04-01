"""

Add on the fly sigma values to a catalogue

AUTHOR: Daniel Farrow (MPE; 2021)

"""

from argparse import ArgumentParser
from astropy.table import Table
from hetdex_api.flux_limits.shot_sensitivity import ShotSensitivity

parser = ArgumentParser(description="Add on the fly sigmas to a catalogue")
parser.add_argument("--wavenpix", default=3, type=int, 
                    help="Take this number of wave pixels either side")
parser.add_argument("--aprad", default=3.5, type=float, 
                    help="Aperture radius in arcseconds")
parser.add_argument("--ascii", action="store_true")
parser.add_argument("datevshot")
parser.add_argument("input")
parser.add_argument("output")
opts = parser.parse_args()


if opts.ascii:
    table = Table.read(opts.input, format="ascii")
else:
    table = Table.read(opts.input)

# Cut to this datevobs
try:
    shotid = int(opts.datevshot.replace("v", ""))
    table = table[table["shotid"] == shotid]
except KeyError as e:
    print("Couldn't find datevobs in catalogue, running on all")

ra = table["ra"]
dec = table["dec"]
wave = table["wave"]

s = ShotSensitivity(opts.datevshot, wavenpix=opts.wavenpix, 
                    rad=opts.aprad)
sigmas, norm = s.get_f50(ra, dec, wave, 1.0, 
                         direct_sigmas=True)

#table_out = Table()
#table_out["detectid"] = table["detectid"]
#table_out["sigma_{:d}_{:2.1f}".format(opts.wavenpix, opts.aprad)] = sigmas
#table_out.write(opts.output)

# use this if no detecid
table["sigma_{:d}_{:2.1f}".format(opts.wavenpix, opts.aprad)] = sigmas
table["norm_{:d}_{:2.1f}".format(opts.wavenpix, opts.aprad)] = norm
table["datevshot"] = opts.datevshot
table.write(opts.output)



