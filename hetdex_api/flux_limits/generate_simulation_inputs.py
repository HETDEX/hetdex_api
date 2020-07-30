"""

Tools to generate inputs to the source
simulations using flux limit cubes. Note the 
ra, dec, redshifts and fluxes for the simulations 
with luminosity function weights weren't generated using 
the tools here, but using a separate python library. The 
wrappers for rs1 etc. where however generated with this code.

AUTHOR: Daniel Farrow (MPE; 2020)

"""
from os.path import basename, dirname, join, isfile
from six import iteritems
from argparse import ArgumentParser
from numpy import (log10, power, unique, isfinite, logical_not, cos, 
                   sin, deg2rad, zeros, ceil, sqrt, square, transpose, savetxt)
from numpy.random import uniform, seed
import astropy.units as u
from astropy.io.fits import Header
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from pyhetdex.coordinates.tangent_projection import TangentPlane
from pyhetdex.het.fplane import FPlane
from hetdex_api.survey import Survey 
from hetdex_api.mask_tools.generate_sky_masks import get_fplane    
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import SensitivityCubeHDF5Container
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube

class HDF5MockContainer(object):
    """
    A mocked SensitivityCubeHDF5Container, so you
    don't actually need to write anything to disk.
    Just used for the simulations.
    """
    def __init__(self):
        self.ifu_slot_dict = {}

    def add_sensitivity_cube(self, datevshot, ifuslot, scube, flush=False):
        """
        Add a sensitivity cube to the mock container.

        Parameters
        ----------
        datevshot : str
            not used, here for
            compatibility
        ifuslot : str
            the IFU identifier
        scube : SensitivityCube
            the sensitivity cube object to
            add
        """
        self.ifu_slot_dict[ifuslot] = scube

    def itercubes(self, datevshot=None):
        """
        Yield ifuslot, sensitivity cube

        Parameters
        ----------
        datevshot : str (optional)
            This is ignored and
            only included to
            match the API of the
            real version

        """
        if datevshot:
            raise Exception("This is only a mock container, datevshot ignored!")

        for ifuslot, scube in iteritems(self.ifu_slot_dict):
            yield ifuslot, scube

    def extract_ifu_sensitivity_cube(self, ifuslot):
        return self.ifu_slot_dict[ifuslot]


def generate_sencube_hdf(datevshot, ra, dec, pa, fplane_output_dir, 
                         nx, ny, nz, ifusize, skip_ifus=["000", 
                         "600", "555", "601", "602", "603", "604", 
                         "610", "611", "612", "613", "614", "615", 
                         "616"], hdf_filename=None):
    """
    Generate an empty real or mock sensitivity HDF5 container, 
    with the proper astrometry in the cubes. Real containters
    are of the SensitivityCubeHDF5Container class and are 
    written to a file. The mock containers are of
    HDF5MockContainer class and do not have a real HDF5
    file - useful for simulations. 

    Parameters
    ----------
    datevshot : str
        the 8 digit YYYYMMDDvSSS date
        of the shot and the shot, used to get the
        correct focal plane file
    ra, dec, pa : float
        the astrometry of the shot
    fplane_output_dir : str
        directory to output fplane files to
    hdf_filename : str (optional)
        if passed, generate a real
        SensitivityCubeHDF5Container
        with this filename. If None
        generate a mock container.
    ifusize : float
        size of x,y of IFU in arcsec
    skip_ifus : list (optional)
        the IFUSLOTS to skip

    Returns
    -------
    hdfcont : SensitivityCubeHDF5Container or HDF5MockContainer
       a real or mock sensivity cube container depending on
       the ``hdf_filename`` parameter
    """ 
    if hdf_filename:
        hdfcont = SensitivityCubeHDF5Container(hdf_filename, mode="w") 
    else:
        hdfcont = HDF5MockContainer()

    # Generate the shot astrometry
    rot = 360.0 - (pa + 90.)
    tp = TangentPlane(ra, dec, rot)

    date = datevshot[:8]
    fplane_bn = "{:s}_fplane.txt".format(date)
    fplane_fn = join(fplane_output_dir, fplane_bn)

    if not isfile(fplane_fn):
       get_fplane(fplane_fn, datestr=str(date))
       fplane = FPlane(fplane_fn)
    else:
       fplane = FPlane(fplane_fn) 

    for ifuslot, ifu in iteritems(fplane.difus_ifuslot):

        if ifuslot in skip_ifus:
            continue

        ifuslot_str = "ifuslot_" + ifuslot
        # Note x, y swapped in focal fplane 
        ra_ifu, dec_ifu = tp.xy2raDec(ifu.y, ifu.x)
        scube = create_sensitivity_cube_from_astrom(ra_ifu.item(), dec_ifu.item(), pa, 
                                                    nx, ny, nz, ifusize)
        hdfcont.add_sensitivity_cube(datevshot, ifuslot_str, scube)

    return hdfcont

def create_sensitivity_cube_from_astrom(racen, deccen, pa, nx, ny, nz, ifusize, 
                                        wrange=[3470.0, 5542.0]): 
    """
    Return an (empty) sensitivity cube object to fill
    with data from simulations later

    Parameters
    ----------
    racen, deccen : float
        the central coordinates of the IFU
    pa : float
        the IFU rotation
    nx, ny, nz : int
        the dimensions of the cube in ra, dec, wave
    ifusize : float
        the length of an IFU side in arcsec
    wrange : array (optional)
        the lower and upper wavelength
        limits in Angstrom
    """

    cards = {}
    cards["NAXIS"] = 3  
    cards["NAXIS1"] = nx
    cards["NAXIS2"] = ny
    cards["NAXIS3"] = nz
    cards["CTYPE1"] = "RA---TAN" 
    cards["CTYPE2"] = "DEC--TAN"
    cards["CTYPE3"] = "Wave "
    cards["CUNIT1"] = "deg " 
    cards["CUNIT2"] = "deg "

    cards["CRPIX1"] = nx/2. + 0.5
    cards["CRPIX2"] = ny/2. + 0.5
    cards["CRPIX3"] = 1.0
   
    coord = SkyCoord(racen*u.deg, deccen*u.deg)
    cards["CRVAL1"] = racen #deg
    cards["CRVAL2"] = deccen #deg
    cards["CRVAL3"] = wrange[0] #AA

    deltapix = (float(ifusize)/nx/3600.0)
 
    # this is rotation in focal plane, maybe not the IFU
    rot = deg2rad(pa)
    cards["CROTA2"] = pa
    cards["CD1_1"] = deltapix*cos(rot)
    cards["CD1_2"] = deltapix*sin(rot)
    cards["CD1_3"] = 0.0
    cards["CD2_1"] = -1.0*deltapix*sin(rot)
    cards["CD2_2"] = deltapix*cos(rot)
    cards["CD2_3"] = 0.0
    cards["CD3_1"] = 0.0
    cards["CD3_2"] = 0.0
    cards["CD3_3"] = (wrange[1] - wrange[0])/nz

    header = Header(cards=cards)
    sigmas = zeros((nz, ny, nx))
    alphas = zeros((nz, ny, nx))

    return SensitivityCube(sigmas, header, None, alphas, aper_corr=1.0, nsigma=1.0)



def rdz_flux_from_hdf_cubes(hdfcont, minfrac=0.2, maxfrac=2.0, nperifu=1000, 
                            add_flim=False, sncut=6.0, logspaced=False):
    """
    Generate ra, dec, redshift and
    flux from a sensivity cube, distribute
    uniformly in x,y,z cube pixels

    Parameters
    ----------
    hdfcont :  SensitivityCubeHDF5Container
        sensitivity cubes for the shot to simulate
    minfrac, maxfrac : float
        the minimum and maximum fraction
        of the flux limit used to define 
        a range of fluxes
    nperifu : int (optional)
        number of sources per IFU
    sncut : float
        the cut on SNR to pass to the flux
        limits API if add_flim=True
    """

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
 
        #if not "086" in ifuslot:
        #    continue

        # Generate randoms in pixel space then transform
        xsize = scube.sigmas.shape[2] 
        ysize = scube.sigmas.shape[1]
        zsize = scube.sigmas.shape[0]
        # Pixels centers start at 0.0 end at xsize - 1
        x = uniform(-0.5, xsize + 0.5, size=nperifu)
        y = uniform(-0.5, ysize + 0.5, size=nperifu)
        # Not redshift!
        z = uniform(0, zsize - 1, size=nperifu)
        r, d, l = scube.wcs.all_pix2world(x, y, z, 0)

        if add_flim:
            flim = scube.get_f50(r, d, l, sncut)
            flim[flim < 1e-30] = 100
            flim[logical_not(isfinite(flim))] = 100
            flims.extend(flim)
        else:
            # fixed value to generate fluxes around
            flim = 5e-17
        
        if logspaced: 
            logf = uniform(log10(flim*minfrac), log10(flim*maxfrac), size=nperifu)
            f = power(10, logf)
        else:
            sqrtf = uniform(sqrt(flim*minfrac), sqrt(flim*maxfrac), size=nperifu)
            f = square(sqrtf)

        ra.extend(r)
        dec.extend(d)
        wave.extend(l)
        flux.extend(f)
        ifuslots.extend([ifuslot]*nperifu)
        ifra, ifdec, l = scube.wcs.all_pix2world(xsize/2 - 0.5, ysize/2. - 0.5, 0, 0)
        ifu_ra.extend([ifra]*nperifu)
        ifu_dec.extend([ifdec]*nperifu)
        print("{:s} {:f} {:f}".format(ifuslot, ifra, ifdec))


    table = Table([ra, dec, wave, flux, ifuslots, ifu_ra, ifu_dec], 
                  names=["ra", "dec", "wave", "flux", "ifuslot", "ifu_ra", "ifu_dec"])

    if add_flim:
        table["flim"] = flim

    print("Generated {:d} sources".format(len(ra)))

    return table

               

def produce_rs1_runfile(fnout, ifuslots, ras, decs, outroot, nsplit):
    """
    Produce a file full of ``rs1`` commands, ready to 
    run as a simulation

    Parameters
    ----------
    fnout : str
        file name to write commands
        to. Can be fed into e.g. the
        jobsplitter command to generate
        SLURM inputs
    ifuslots, ras, decs : list or array 
        The IFUSLOTS, right ascensions (deg.)
        and declinations (deg.) of the IFUs
        you plan on adding simulated sources
        to
    outroot : str
        runs basename(outroot) on this to
        generate the root of the output to
        pass to rs1
    nsplit : int
        indicates which of the multiple
        realisations this rs1 corresponds
        to

    Example output
    -------------
    rs1 25.51297  0.56949 35 4505 50 ifuslot_086_0 20190105v013 1.7 3 3.5 0. 1 110

    """
    fitradecspmode = 110
    boutroot = basename(outroot)
    n = 0
    with open(fnout, "w") as fout:
        for ifu, r, d in zip(ifuslots, ras, decs):
            fout.write("rs1 {:8.5f} {:8.5f} 35 4505 50 ".format(r, d) +
                       "{:s} {:s} 1.7 3 3.5 0.2 7 {:d}".format(ifu, 
                                                              boutroot, 
                                                              fitradecspmode)) 
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

def split_into_ifus(table, unique_fields, unique_ifus, outdir,
                    NMAX=1000, nsplit_start=0, nsplit=1):
    """
    Split a catalogue into input catalogues
    and commands on single IFUs. Outputs
    catalogues and files full of calls
    to Karl's HETDEX simulation code

    Parameters
    ----------
    unique_fields, unique_ifus : array
        fields and IFUs to split up
    outdir : str
        output directory
    NMAX : int (optional)
        maximum source per IFU,
        further splits if bigger than this
 
    """

    # Maximum sources per ifu in one sim run
    # Second part of the filename
    fn2s = []
    ras = []
    decs = []
    for field in unique_fields:
        table_field = table[table["field"] == field]
        for ifuslot in unique_ifus:

            table_ifu = table_field[table_field["ifuslot"] == ifuslot]
            # Split into sub sets if too big to simulate at once
            nsplit = int(ceil(len(table_ifu)/NMAX))
            for isplit in range(nsplit):
               ttable = table_ifu[isplit*NMAX : (isplit + 1)*NMAX]

               date = field[:-4]
               fplane_fn = "{:s}_fplane.txt".format(date)
               if not isfile(fplane_fn):
                   get_fplane(fplane_fn, datestr=date)
                   fplane = FPlane(fplane_fn)
               else:
                   fplane = FPlane(fplane_fn) 


               ifu = fplane.by_id(ifuslot[-3:], "ifuslot")               
               

               if isplit + nsplit_start > 0: 
                   fn2 = ifuslot + "_{:d}".format(isplit + nsplit_start)
               else:
                   # ingestion script compatible fn
                   fn2 = "{:03d}_{:s}_{:s}".format(ifu.specid, ifu.ifuslot, ifu.ifuid)

               outfn = join(outdir, "{:s}_{:s}.input".format(field, fn2))
               savetxt(outfn, transpose([ttable["ra"], ttable["dec"],
                                         ttable["wave"], ttable["flux_extincted"]
                                        ]))
               ras.append(ttable["ifu_ra"][0]) 
               decs.append(ttable["ifu_dec"][0])
               fn2s.append(fn2)
 
        produce_rs1_runfile(join(outdir, "{:s}_{:d}.run".format(field, nsplit_start)), 
                            fn2s, ras, decs, field, nsplit)



def generate_sources_to_simulate(args = None):
    """
    Generate source simulation inputs
    spread evenly over the sensitivity cube
    pixels, from an input sensitivity cube. Not
    suitable as a random catalogue, as not uniform
    in redshift - good for probing the detection
    efficiency etc. over redshift though

    """  

    parser = ArgumentParser(description="Generate inputs for source simulations, uniform in x,y,z datacube coords")
    parser.add_argument("--nsplit-start", default=0, type=int, help="What number to start the file split labeling at")
    parser.add_argument("--ifusize", help="Size of IFU, in arcsec", default=52.)
    parser.add_argument("--nperifu", help="Number of sources to add per IFU", default=1000, type=int)
    parser.add_argument("--flimcut", help="Don't simulate source above this", type=float, default=1e-15)
    parser.add_argument("--fix-flux", help="If not None, set all sources to this flux value", default=None, type=float)
    parser.add_argument("--frac-range", default=[0.2, 2.0], nargs=2, type=float, help="The range in flux as a fraction of flim")
    parser.add_argument("--nsplit", default=1, help="Number of jobs per line in .run file")
    parser.add_argument("--nmax", default=1000, help="Maximum number of sources in one sim", type=int)
    parser.add_argument("filelist", help="Ascii file with list of sensitivity HDF5 files or date shots")
    parser.add_argument("outdir", help="Directory for output files")
    opts = parser.parse_args(args=args) 

    survey_obj = None
    tables = []
    with open(opts.filelist, "r") as fp:
        for line in fp:
            input_ = line.strip()
            if ".h5" in input_:
                print("Assuming {:s} is an HDF5 file".format(input_))
                field = filename.strip("_sensitivity_cube.h5")
                sencube_hdf = SensitivityCubeHDF5Container(filename)
                add_flim = True
            else:
                print("Assuming {:s} is date shot".format(input_))
                field = input_        
                add_flim = False

                if type(survey_obj) == type(None):
                    survey_obj = Survey("hdr2")

                shot = survey_obj[survey_obj.datevobs == field]
                sencube_hdf = generate_sencube_hdf(shot.datevobs[0], shot.ra[0], 
                                                   shot.dec[0], shot.pa[0], opts.outdir,
                                                   31, 31, 1036, 
                                                   opts.ifusize)

        ttable = rdz_flux_from_hdf_cubes(sencube_hdf, add_flim=add_flim, minfrac=opts.frac_range[0], 
                                         maxfrac=opts.frac_range[1], nperifu=opts.nperifu)		 
        ttable["field"] = field
        tables.append(ttable)
 

    table = vstack(tables)
    table["id"] = range(len(table))

    if opts.fix_flux:
        table["flux"] = opts.fix_flux


    if add_flim:
        table = table[table["flim"] < opts.flimcut]
        print("After flimcut left with {:d} sources!".format(len(table)))

    #table.write("full_input_{:d}.txt".format(opts.nsplit_start), format="ascii.commented_header")

    # Split into IFUS and fields
    unique_ifus = unique(table["ifuslot"])
    unique_fields = unique(table["field"])
 
    # Second part of the filename
    fn2s = []
    split_into_ifus(table, unique_fields, unique_ifus,
                    opts.outdir, NMAX=opts.nmax, 
                    nsplit_start=opts.nsplit_start,
                    nsplit=opts.nsplit)


if __name__ == "__main__":
    generate_sources_to_simulate()
