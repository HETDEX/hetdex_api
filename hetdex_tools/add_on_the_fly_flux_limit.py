"""

Add on the fly sigma values to a catalogue

AUTHOR: Daniel Farrow (MPE; 2021)
Updated: Daniel Farrow (Hull; 2025)

Some AI generated code included as indicated in the function
headers.

"""
import h5py
import cProfile
from argparse import ArgumentParser
from astropy.table import Table
from hetdex_api.flux_limits.shot_sensitivity import ShotSensitivity

def h5_to_astropy_table(file_path, verbose=True):
    """
    Converts datasets in an HDF5 file to an Astropy Table.
    
    Parameters:
        file_path (str): Path to the .h5 file.
    
    Returns:
        astropy.table.Table: Table containing the datasets.
        
    Reference:
        Original: AI-generated code, Microsoft Copilot (May 2025), prompted by
        Daniel Farrow (2025), prompt: 'A python function to convert an 
        h5 file to a astropy table'
        
        This version: Adapted to make it run by Daniel Farrow (2025)
    """
    try:
        h5_file = h5py.File(file_path, 'r')
        data_dict = {}
        dataset = h5_file["simdata"]
        with h5_file:
            for name, data in dataset.items():
                if verbose:
                    print(name, type(data), data[:10])
                data_dict[name] = data
                table = Table(data_dict)      
    except Exception as e:
        h5_file.close()
        raise(e)

    # Remove column with troublesome data type
    table.remove_column("ampin")
    
    return table

# Following https://stackoverflow.com/questions/
#           54465048/profiling-a-python-3-6-module-from-the-command-line
pr = cProfile.Profile()
pr.enable()

parser = ArgumentParser(description="Add on the fly sigmas to a catalogue")
parser.add_argument("--wavenpix", default=3, type=int, 
                    help="Take this number of wave pixels either side")
parser.add_argument("--aprad", default=3.5, type=float, 
                    help="Aperture radius in arcseconds")
parser.add_argument("--ascii", action="store_true")
parser.add_argument("--hdf5", action="store_true")
parser.add_argument("--header", help="Header for the ascii file. "
                                     "col1: anything col2: column " 
                                     "name")
parser.add_argument("datevshot")
parser.add_argument("input")
parser.add_argument("output")
opts = parser.parse_args()

if opts.header:
    names = [] 
    with open(opts.header, 'r') as fp: 
        for line in fp: 
            names.append(line.strip().split()[1]) 

    table = Table.read(opts.input, format="ascii", 
                       names=names)

elif opts.ascii:
    table = Table.read(opts.input, format="ascii")
elif opts.hdf5:
    table = h5_to_astropy_table(opts.input)
else:
    table = Table.read(opts.input)

#table = table[:50000]

# Cut to this datevobs
try:
    shotid = int(opts.datevshot.replace("v", ""))
    table = table[table["shotid"] == shotid]
except KeyError as e:
    try:
        table = table[table["datevobs"] == opts.datevshot]
    except KeyError as e:
        print("Couldn't find datevobs in catalogue, running on all")

# Remember to add '_in' for the input positions
ra = table["RAin"]
dec = table["Decin"]
wave = table["wavein"]
flux = 1e-17*table['fluxin']

# Unfinished !!!
s = ShotSensitivity(opts.datevshot, wavenpix=opts.wavenpix, 
                    rad=opts.aprad)

sigmas, norm = s.get_f50(ra, dec, wave, 1.0, 
                         direct_sigmas=True)

f50s = s.f50_from_noise(sigmas, wave, sncut, 
                        linewidth = linewidth)

c = s.return_completeness(flux, ra, dec, wave, 
                          sncut, f50s=f50s,
                          linewidth = linewidth)

#table_out = Table()
#table_out["detectid"] = table["detectid"]
#table_out["sigma_{:d}_{:2.1f}".format(opts.wavenpix, opts.aprad)] = sigmas
#table_out.write(opts.output)

table["sigma_{:d}_{:2.1f}".format(opts.wavenpix, opts.aprad)] = sigmas
table["norm_{:d}_{:2.1f}".format(opts.wavenpix, opts.aprad)] = norm
table["datevshot"] = opts.datevshot

table.write(opts.output)

pr.disable()
pr.print_stats(sort="cumtime")

