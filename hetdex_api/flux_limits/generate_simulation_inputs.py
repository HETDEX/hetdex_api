"""

Generate input to the source
simulations using flux limit cubes

AUTHOR: Daniel Farrow (MPE; 2020)

"""
from os.path import basename
from argparse import ArgumentParser
from numpy import log10, power, unique, isfinite, logical_not
from numpy.random import uniform, seed
from astropy.table import Table, vstack
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import SensitivityCubeHDF5Container


def rdz_flux_from_hdf_cubes(filename, minfrac=0.1, maxfrac=20.0, nperifu=1000):
    """
    Generate ra, dec, redshift and
    flux from

    Parameters
    ----------
    filename : str
        the filename of the HDF5
    minfrac, maxfrac : float
        the minimum and maximum fraction
        of the flux limit used to define 
        a range of fluxes
    """

    hdfcont = SensitivityCubeHDF5Container(filename)

    ra = []
    dec = []
    wave = []
    flux = []
    flims = []
    ifuslots = []
    tables = []
    ifu_ra = []
    ifu_dec = []

    for ifuslot, scube in hdfcont.itercubes():

        # Generate randoms in pixel space then transform
        xsize = scube.f50vals.shape[2] 
        ysize = scube.f50vals.shape[1]
        zsize = scube.f50vals.shape[0]
        # Pixels centers start at 0.0 end at xsize - 1
        x = uniform(-0.5, xsize + 0.5, size=nperifu)
        y = uniform(-0.5, ysize + 0.5, size=nperifu)
        # Not redshift!
        z = uniform(0, zsize - 1, size=nperifu)
        r, d, l = scube.wcs.all_pix2world(x, y, z, 0)
        flim = scube.get_f50(r, d, l)

        flim[flim < 1e-30] = 100
        flim[logical_not(isfinite(flim))] = 100

        logf = uniform(log10(flim*minfrac), log10(flim*maxfrac), size=nperifu)
        f = power(10, logf)
       
        ra.extend(r)
        dec.extend(d)
        wave.extend(l)
        flux.extend(f)
        flims.extend(flim)
        ifuslots.extend([ifuslot]*nperifu)
        ifra, ifdec, l = scube.wcs.all_pix2world(xsize/2 - 0.5, ysize/2. - 0.5, 0, 0)
        ifu_ra.extend([ifra]*nperifu)
        ifu_dec.extend([ifdec]*nperifu)
        print("{:s} {:f} {:f}".format(ifuslot, ifra, ifdec))

    table = Table([ra, dec, wave, flux, flims, ifuslots, ifu_ra, ifu_dec], 
                  names=["ra", "dec", "wave", "flux", "flim", "ifuslot", "ifu_ra", "ifu_dec"])

    print("Generated {:d} sources from {:s}".format(len(ra), filename))

    return table

def produce_rmksource_runfile_old(fnout, table, nsplit):
    """
    Produce a file full of rmksource commands

    Example output
    --------------
    rmksource 1165 150.194 2.348 5400 6e-17 20190204v014

    New Example output
    ------------------
    rmksource2 217.788452 51.7916641 008_093_054 20190207v024

    """

    n = 0
    with open(fnout, "w") as fout:
        for line in table:
            fout.write("rmksource {:d} {:8.5f} {:8.5f} {:4.1f} {:e} {:s}".format(
                       line["id"], line["ra"], line["dec"],
                       line["wave"], line["flux"], line["field"]))
            n = n + 1
            if n == nsplit:
                fout.write("\n") 
                n = 0
            else:
                fout.write("; ")
                

def produce_rmksource_runfile(fnout, ifuslots, ras, decs, outroot, nsplit): 
    """
    Produce a file full of rmksource commands

    New Example output
    ------------------
    rmksource 217.788452 51.7916641 008_093_054 20190207v024

    """
    boutroot = basename(outroot)
    n = 0
    with open(fnout, "w") as fout:
        for ifu, r, d in zip(ifuslots, ras, decs):
            fout.write("rmksource {:8.5f} {:8.5f} {:s} {:s}".format(r, d, ifu, boutroot)) 

            n = n + 1
            if n == nsplit:
                fout.write("\n") 
                n = 0
            else:
                fout.write("; ")
 
 
if __name__ == "__main__":
  
    parser = ArgumentParser(description="Generate inputs for source simulations")
    parser.add_argument("--flimcut", help="Don't simulate source above this", type=float, default=1e-15)
    parser.add_argument("--nsplit", default=1, help="Number of jobs per line in .run file")
    parser.add_argument("filenames", nargs="+", help="Sensitivity HDF5 files to generate sources for")
    parser.add_argument("outroot", help="Root for output files")
    opts = parser.parse_args() 

    tables = []
    for filename in opts.filenames:    

        ttable = rdz_flux_from_hdf_cubes(filename)
        field = filename.strip("_sensitivity_cube.h5")
        ttable["field"] = field
        tables.append(ttable)
 

    table = vstack(tables)
    table["id"] = range(len(table))
    table = table[table["flim"] < opts.flimcut]
    print("After flimcut left with {:d} sources!".format(len(table)))

    table.write("{:s}_full_input.txt".format(opts.outroot), format="ascii.commented_header")

    # Split into IFUS
    unique_ifus = unique(table["ifuslot"])

    ras = []
    decs = []
    for ifuslot in unique_ifus:
        table_ifu = table[table["ifuslot"] == ifuslot]
        table_ifu.write("{:s}_{:s}.input".format(opts.outroot, ifuslot), include_names=["ra", "dec", "wave", "flux"],
                        format="ascii.no_header")
        ras.append(table_ifu["ifu_ra"][0]) 
        decs.append(table_ifu["ifu_dec"][0])

    produce_rmksource_runfile("{:s}.run".format(opts.outroot), unique_ifus, ras, decs, opts.outroot, opts.nsplit)

