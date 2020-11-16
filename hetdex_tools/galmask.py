"""
Tools for creating elliptical galaxy masks, with parameters from the RC3 galaxy catalog

Created on 2020/10/28
Last Updated: 2020/11/9

Updated to use config.py path to RC3 cat
Updated to use hetdex_api.detections to grab detections

@author: John Feldmeier (YSU)
"""

from __future__ import print_function

import numpy as np

import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
import os, time
import pylab as pl

from os.path import expanduser
from numpy.random import random

from astropy import units as u
from astropy import coordinates
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astroquery.skyview import SkyView

from regions import EllipseSkyRegion, EllipsePixelRegion
from regions import CircleSkyRegion, CirclePixelRegion
from regions import LineSkyRegion, LinePixelRegion
from regions import PixCoord

from hetdex_api.config import HDRconfig
from hetdex_api.detections import Detections

config = HDRconfig()

# # Usage cases for galmask.py
#
# 1. Reading in needed tables:
#
# To use the code at all, you will need to read in the RC3 catalog of sources
# within the HETDEX fields.  This can be done using the following methods:
#
#         tboth = read_rc3_tables_old() or tboth = read_rc3_tables()
#
# You will also need to read in a table of HETDEX detections:
#
#         tdetect = read_detect_old() or tdetect = read_detect()
#
# 2. How the elliptical regions work:
#
# The regions use the *astropy.regions* formalism
# (see https://astropy-regions.readthedocs.io/en/latest/),
# and use an elliptical isophote derived from a galaxy catalog.  Currently,
# the only catalog used is the RC3 (Third Reference Catalog of Bright Galaxies),
# but other galaxy catalogs can be added in the future, since all that is
# needed are the parameters of the ellipse.
#
# The table is a comma separated variable file that contains the galaxy
# shape information.  The key parameters needed are:
# the **Coords* string, that gives the center of the galaxy in J2000 ICRS
# coordinates, the **SemiMajorAxis** and **SemiMinorAxis** parameters
# of the ellipse, and a **PositionAngle**, defined in the standard way
# (counter-clockwise from North through East on the sky).  The ellipse can
# be scaled by the **d25scale parameter** in almost all of the methods to make
# the ellipse larger (or smaller).  d25scale = 1.0 is the isophote that
# containts 25 magnitudes per square arcsecond in the B band, as defined by
# the RC3 catalog.  In some cases, the RC3 values for the ellipses were absent
# or unreliable - these were fixed either by referring to the NED values, or
# by visual inspection.
#
# 3. Visualizing the galaxy ellipses:
#
# To visualize what the ellipses look like for a particular galaxy, you can use
# the following two methods:
#
# show_ellipse_ring_source(tboth, ind1, coord, rings, 'out.png', save_file=False, show_notebook=True)
#
# and
#
# show_ellipse_ring_detect(tboth, ind1, tdetect, rings, 'out.png', save_file=False, show_notebook=True)
#
#
# show_ellipse_ring_source currently downloads an image from the Skyview web site
# (https://skyview.gsfc.nasa.gov/current/cgi/titlepage.pl), using *astroquery* and
# using matplotlib, creates a plot that overlays the ellipses of a specific
# galaxy, as well as give the position for a single source. The input
# parameters are:
#
# *tboth* refers to the RC3 table loaded previously,
# *ind1* is the integer index referring to a specific galaxy,
# *coord* is a single *SkyCoord* object,
# *rings* refers to a numpy array showing the ellipses to be plotted, with the last value giving the size of
# an outer circle, called *rlimit*.  It's presumed that all real sources belonging to the galaxy are inside of
# rlimit, and it's better for the rlimit value to be larger than you would expect (since it's only used for plotting
# and intial selection.  *rings* is given in units of d25, and a useful example is
# *rings = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.5])*.
#
# The image can be saved to an output file (*save_file = True*) where *out.png*
# would be replaced by the filename you wanted, or it can be displayed to a Jupyter
# notebook (*show_notebook = True*).  You can display things both ways, if you wish.
#
# show_ellipse_ring_detect is almost identical to show_ellipse_ring_source, but
# it instead plots all data from *tdetect* that lie within the outermost
# circle on the image. This second method calls another useful method, called *numrings*
#
#         num_rings(t, index, tdetect, rings)
#
# which will allow you to determine the number of *tdetect* sources within
# a specific scaled elliptical isophote, making radial plots possible.  numrings
# has the following parameters
# where *t* is the RC3 galaxy table,
# *index* is the particular value in that table,
# *tdetect* is as before, and *rings* is as before.
#
# 4. Searching for objects within a single (known) galaxy ellipse
#
# To search whether a number of coordinates are within the ellipse of a
# *single* known galaxy, use the following method:
#
# isclose, name, zgal = onegal_flag_from_coords(coords, tboth, index, d25scale = 1.0)
#
# where *coords* are SkyCoords as before, *tboth* is the RC3 table, index
# refers to the specific galaxy, and *d25* scale is the scaling for the ellipses.
# The method returns three output numpy arrays:
#
# - *isclose* is a Boolean array that indicates whether a source was inside
#   the scaled ellipse or not.
# - *name* gives the name of the galaxy if it is inside an ellipse
# - *zgal* gives the NED redshift for that corresponding galaxy
#   (for comparison purposes)
#
# If there is no match, the values will be 'isclose = False, name = None, zgal = None'
#
# 5. Searching for objects within the entire galaxy table, checking all ellipses
#
# To search whether a number of coordinates are within the ellipse of all
# galaxies in the table, use the following method:
#
# isclose, name, zgal = gals_flag_from_coords(coords, tboth, d25scale=1.0, nmatches=3)
#
# Where the *isclose* *tboth* and *d25matches* are as before.  This method will
# loop through the entire catalog (defined by tboth), and will check the closest
# *nmatch* sources.  The *nmatches* parameter is due to the fact that ellipses
# might overlap (this does happen for interacting galaxy pairs).  If you don't
# want to check for overlaps, set *nmatches = 1* (the default)
#
# 6. Searching the entire table for objects inside each ellipse:
#
# outTest, outNames, outRedshifts = check_all_large_gal(tboth, tdetect, scale, verbose=False)
#
# where *outTest* is a Boolean array for each source, similar to *isclose*,
# *outNames* is similar to *name* above, and
# *outRedshifts* is similar to *zgal* above.  Since this checks every single source in
# a large catalog (over one million sources in HDR 2.10), against over 900 galaxy regions,
# this will take a while (about 250s on a moderate power laptop).
#
# 7. Usage examples:
#
# There are a number of test cases, ranging from simple to more complex given below as 'test1', through 'test10'


def read_rc3_tables():
    """
    Read in the RC3 table created by John Feldmeier for HETDEX.  
    
    2020.11.04 - not yet ready - will require some config parameter, probably
    in 
    The key parameters in the table are:
    'Coords'        - string representation of the galaxy center
    'SemiMajorAxis' - the semi-major axis of the galaxy, measured in arcminutes
    'SemiMinorAxis' - the semi-minor axis of the galaxy, measured in arcminutes
    'PositionAngle' - The angle of the major axis, measured from North through 
                      east, measured in degrees.
    'NEDRedshift'   - The redshift of the galaxy, taken from NED (if available)
    """
    tboth = Table.read(config.rc3cat, format="ascii.csv")
    return tboth


def check_gal_table():
    """Check to ensure that the galaxy table is not nonsensical.  Do later."""
    pass


def read_detect(version='2.1.2'):
    detects = Detections(curated_version=version)
    detect_table = detects.return_astropy_table()
    return detect_table


def create_dummy_wcs(c1):
    """
    Create a simple fake WCS for use in the regions subroutine.
    c1 is a center SkyCoord, so that the pixel values do not go to very large or small values.
    """

    # version from Erin Mentuch Cooper, with small changes by John Feldmeier

    gridsize = np.float64(
        3600.0
    )  # image size in arcsec.  Somewhat arbitrary, but make large enough to enclose
    gridstep = np.float64(
        1.0
    )  # make the fake pixel scale one arsec per pixel for convenience

    # Create coordinate center
    ra_cen = c1.ra.deg
    dec_cen = c1.dec.deg

    ndim = np.int(2 * gridsize / gridstep + 1)
    center = np.int(ndim / 2)
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [ra_cen, dec_cen]
    w.wcs.crpix = [center, center]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cdelt = [-gridstep / gridsize, gridstep / gridsize]
    return w


def create_ellreg(t, index, d25scale=1.0):
    """
    Creates an elliptical sky region from astropy.regions, with info from RC3 catalog

    t - a table of galaxy regions, similar to that found in read_rc3_tables
    index - the value of the table to be used
    scale - a scaling factor.  Leaving at 1.0 means to use the D25 isophote, 
    and it's likely we will want to scale this larger.
    """

    coords = SkyCoord(t["Coords"][index], frame="icrs")

    # The ellipse region uses the major and minor axes, so we have to multiply by
    # two first, before applying any user scaling.

    major = (t["SemiMajorAxis"][index]) * np.float64(2.0) * d25scale * u.arcmin
    minor = (t["SemiMinorAxis"][index]) * np.float64(2.0) * d25scale * u.arcmin
    pa = (t["PositionAngle"][index]) * u.deg
    ellipse_reg = EllipseSkyRegion(center=coords, height=major, width=minor, angle=pa)

    return ellipse_reg

def create_gal_region(galaxy_cat, row_index=None, pgcname=None, d25scale=1.):
    """
    Similar to ellreg above but can take a galaxy name as input.
    
    Create galaxy ellipse region using an input galaxy catalog_table (likely need
    to change this API as it heavily follows the RC3 catalog format)
    
    Parameters
    ----------
    galaxy_catalog_table: an astropy table
        table of nearby galaxy positions and sizes. Must have a central
        coordinate and the SemiMajorAxis, SemiMinorAxis, Position Angle
        info as table column names
    row_index: int
        a row index in the galaxy catalog to create the ellipse for
    pgcname: str
        the PGCNAME in the RC3 cat. This is a string. eg. "PGC 43255" for NGC 4707
    d25scale: float
        how many times D25 should to scale the region
    """
    if row_index is not None:
        index = row_index
    elif pgcname is not None:
        index = np.where(galaxy_cat['PGC'] == pgcname)[0][0]
        
    coords = SkyCoord(galaxy_cat['Coords'][index],frame='icrs')

    # The ellipse region uses the major and minor axes, so we have to multiply by
    # two first, before applying any user scaling.
    
    major = (galaxy_cat['SemiMajorAxis'][index]) * np.float64(2.) * d25scale * u.arcmin
    minor = (galaxy_cat['SemiMinorAxis'][index]) * np.float64(2.) * d25scale * u.arcmin
    pa    = (galaxy_cat['PositionAngle'][index]) * u.deg
    ellipse_reg = EllipseSkyRegion(center=coords, height=major, width=minor, angle=pa)
    
    return ellipse_reg
                                                                            

def onegal_flag_from_coords(coords, t, index, d25scale=1.0):
    """ Search a specific galaxy ellipse, and see if a single source lies within that region.
    If so, return True, the name of the galaxy, and the NED redshift.  Otherwise, return
    False, None, None."""

    # create fake WCS for regions use
    mywcs = create_dummy_wcs(coords)

    flag = False
    source_name = None
    source_redshift = None

    ellipse = create_ellreg(t, index, d25scale=d25scale)
    if ellipse.contains(coords, mywcs):
        flag = True
        source_name = t["PGC"][index]
        source_redshift = t["NEDRedshift"][index]

    return (flag, source_name, source_redshift)


def gals_flag_from_coords(coords, galaxy_cat, d25scale=1.0, nmatches=1):
    """
    Returns a boolean flag value to mask sources near large galaxies
    Parameters
    ----------
    coords
        an astropy.coordinates SkyCoord object
    galaxy_cat
        an astropy table containing the large galaxy parameters. This
        is catered for the RC3 catalog stored in config.rc3cat
    d25
        The scaling of ellipses.  1.0 means use the ellipse for D25.
        Experimentation showed a value of 1.75 might be more appropriate
    nmatches
        the closest nmatches are searched for.  nmatches = 1 means
        search the closest coordinate only.  nmatches = 3 is recommended

    Returns
    -------
    flag - boolean
        True if the source lies within the ellipse defined 
        False if the source is not within the scaling
    source_name - string
        The t['PGC'] value of the galaxy matched to
    source_redshift 
        The t['NEDRedshift'] value of the galaxy matched to

    """

    gal_coords = SkyCoord( galaxy_cat["Coords"])

    # calculate angular distances to all of the sources, and pick out the n closest ones
    d2d = coords.separation(gal_coords)
    ind1 = np.argsort(d2d)  # sort from closest to farthest away
    id_close = ind1[0:nmatches]

    # create fake WCS for regions use
    mywcs = create_dummy_wcs(coords)

    flag = False
    source_name = None
    source_redshift = None

    for idnum in id_close:
        ellipse = create_ellreg(galaxy_cat, idnum, d25scale=d25scale)
        if ellipse.contains(coords, mywcs):
            flag = True
            source_name = galaxy_cat["PGC"][idnum]
            source_redshift = galaxy_cat["NEDRedshift"][idnum]

    return (flag, source_name, source_redshift)


def show_ellipse_ring_source(
    t,
    index,
    coords,
    rings,
    outname,
    image_survey="SDSSg",
    save_file=True,
    show_notebook=False,
):

    """For a single HETDEX near a galaxy, overplot the ellipsoidal rings from the RC3, along
    with the source object given by coords onto an image for reference.  The image comes
    from SkyView, and so it requires connectivity.
    
    Input:
       
    t - the RC3 both survey, read in with read_rc3_tables_old(), or read_rc3_tables()
    
    index - the index of the galaxy to be plotted
    
    coords - the SkyCoord coordinate of the HETDEX source to be plotted
    
    rings - np.array with scaling factors, given in units of D25 that show the ellipses to
    be overplotted, which are shown in blue.  The last value gives the green outermost 
    circle, called r_limit, that is intended to be much larger than the mask we will 
    be using.  A recommended list is rings = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.5])
    
    image_survey - Which imaging survey is called from Skyview
    
    Output:
    
    The image is saved as a png file (if save_file = True) or it can be shown
    in a Jupyter Notebook, (if show_notebook = True).  
    """

    import warnings

    warnings.filterwarnings("ignore")  # get rid of annoying "color"
    # matplotlib warnings for right now

    # Parameters that might need to be changed:

    source_size = 1.5 * u.arcsec  # How large the source circle is in the plot.

    # Use the outermost scale as the radius.  Make this somewhat larger than the

    source_scale = rings[-1]

    # Read in galaxy parameters from the RC3 table.
    c1 = SkyCoord(t["Coords"][index], frame="icrs")
    pgc = t["PGC"][index]
    name = t["Name1"][index]
    major = t["SemiMajorAxis"][index] * u.arcmin
    minor = t["SemiMinorAxis"][index] * u.arcmin
    pa = t["PositionAngle"][index] * u.deg

    rlimit = major * source_scale
    radouter = (
        rlimit * 2.1
    )  # The size of the image displayed, slightly larger than 2 times
    nrings = rings.size

    # Create title information on ellipse parameters for the plot and the source coordinate
    s1 = str(pgc) + " " + str(name) + " " + t["Coords"][index] + "\n"
    s2 = "({0:0.03f} {1:0.03f} PA: {2:0.01f})".format(major, minor, pa) + "\n"
    s2a = "Source Coord: " + coords.to_string("hmsdms")
    s3 = s1 + s2 + s2a

    # Create a one arcminute scale bar on the East edge of the image, pointing North-South
    temp1 = c1.ra + (major * (source_scale * 0.925))
    temp2 = c1.dec - (0.5 * u.arcmin)
    start_sky = SkyCoord(temp1, temp2)
    temp3 = c1.ra + (major * (source_scale * 0.925))
    temp4 = c1.dec + (0.5 * u.arcmin)
    end_sky = SkyCoord(temp3, temp4)
    scale_bar = LineSkyRegion(start=start_sky, end=end_sky)

    # Now, create the ellipse regions list to be overplotted.  Loop over all but the last ring.
    regionList = []
    for i in range(0, nrings - 1):
        scale = rings[i]
        ellipse = create_ellreg(t, index, d25scale=scale)
        regionList.append(ellipse)

    # Create the outerRing circle, which shows the radius that we searched to.
    outerRing = []
    tempreg = CircleSkyRegion(center=c1, radius=rlimit)
    outerRing.append(tempreg)

    # Create a region from the source we want to overplot
    sourceCircles = []

    if np.size(coords) > 1:
        for coord in coords:
            tempreg = CircleSkyRegion(center=coord, radius=source_size)
            sourceCircles.append(tempreg)
    else:
        tempreg = CircleSkyRegion(center=coords, radius=source_size)
        sourceCircles.append(tempreg)

    try:
        imglist = SkyView.get_images(position=c1, radius=radouter, survey=image_survey)

    # If SDSS g is not found, try using the DSS2 image
    except:
        imglist = SkyView.get_images(position=c1, radius=radouter, survey="DSS2 Blue")

    # the returned value is a list of images, but there is only one
    img = imglist[0]

    # 'img' is now a fits.HDUList object; the 0th entry is the image
    mywcs = wcs.WCS(img[0].header)

    # Begin plotting here
    fig = plt.figure(figsize=(8.0, 8.0))
    ax = fig.add_subplot(111)
    s4 = "RA - Object: {0:d} ".format(index) + str(t["Notes"][index]) + "\n"
    s5 = "R limit: {0:0.03f}".format(rlimit)
    s7 = s4 + s5

    ax.set_xlabel(s7)
    ax.set_ylabel("DEC")
    ax.set_title(s3)

    # Plot the Skyview fits file
    ax.imshow(
        img[0].data,
        cmap="gray_r",
        interpolation="none",
        origin="lower",
        norm=pl.matplotlib.colors.LogNorm(),
    )

    # overplot the one arcminute scale bar
    pixel_region = scale_bar.to_pixel(mywcs)
    pixel_region.plot(ax=ax, color="black", linewidth=4)

    # overplot the RC3 ellipses that were chosen
    for reg in regionList:
        pixel_region = reg.to_pixel(mywcs)
        pixel_region.plot(ax=ax, color="blue", linewidth=2)

    # overplot source position as purple circle
    for reg in sourceCircles:
        pixel_region = reg.to_pixel(mywcs)
        pixel_region.plot(
            ax=ax, edgecolor="purple", facecolor="purple", linewidth=5, fill=True
        )

    # overplot search radius
    for reg in outerRing:
        pixel_region = reg.to_pixel(mywcs)
        pixel_region.plot(ax=ax, color="green", linewidth=4)

    # Show the output
    if save_file:
        fig.savefig(outname, dpi=400)

    if show_notebook:
        plt.show()

    plt.close()


def num_rings(t, index, tdetect, rings):
    """
    Find the total number of HETDEX sources from tdetect inside a set of scaled ellipsoidal rings.
    
    Input
    
    t - the RC3 catalog table
    index - the index corresponding to the specific galaxy
    tdetect - the HETDEX ctalog
    rings = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.5])
    """
    # Read in specific galaxy parameters from the table.

    c1 = SkyCoord(t["Coords"][index], frame="icrs")
    pgc = t["PGC"][index]
    name = t["Name1"][index]
    major = t["SemiMajorAxis"][index] * u.arcmin
    minor = t["SemiMinorAxis"][index] * u.arcmin
    pa = t["PositionAngle"][index] * u.deg

    source_scale = rings[-1]
    rlimit = major * source_scale
    nrings = len(rings)
    temprings = np.zeros(nrings)

    # Now, create a list of sources near the galaxy from the detect table.

    catalog = SkyCoord(tdetect["ra"], tdetect["dec"], unit=(u.deg, u.deg))

    d2d = c1.separation(catalog)  # Find all of separations
    catalogmsk = d2d < rlimit  # Create Boolean mask of objects inside of this
    near_sources = catalog[catalogmsk]

    # Create fake WCS, centered on the galaxy
    mywcs = create_dummy_wcs(c1)

    # Loop over all but the last ring, and find the number of sources within each scaled ellipse

    for i in range(0, nrings - 1):
        scale = rings[i]
        ellipse = create_ellreg(t, index, d25scale=scale)
        nellipse = 0
        for skycoord in near_sources:
            if ellipse.contains(skycoord, mywcs):
                nellipse += 1
        temprings[i] = nellipse

    # Finally, for the last ring, find the number of sources within the large outer circle
    # Create the outerRing circle, which shows the radius that we searched to.

    ntotal = 0
    outerRing = CircleSkyRegion(center=c1, radius=rlimit)
    for skycoord in near_sources:
        if outerRing.contains(skycoord, mywcs):
            ntotal += 1

    temprings[-1] = ntotal

    outrings = temprings

    return outrings


def show_ellipse_ring_detect(
    t,
    index,
    tdetect,
    rings,
    outname,
    image_survey="SDSSg",
    save_file=True,
    show_notebook=False,
):

    """For an object near a galaxy, overplot the ellipsoidal rings from the RC3, along
    
    Download the SDSS g image of each galaxy, and overlay the ellipses requested,
       with the source object given by coords.
       
    
    image_survey = The default SDSSg Which imaging survey is called from Skyview
    The image is saved as a png file (if save_file = True) or it can be shown
    in the Jupyter Notebook, if show_notebook = True.  The four ellipses 
    are D25, 2 * D25, and 3 * D25, 4 * D25
    rings = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.5])"""

    # rings = np.array with scaling factors 1.0, 2.0, 3.0, 4.0, 4.5
    # outrings = number of objects within that ellipse scaling, with the last entry being the circle

    import warnings

    warnings.filterwarnings("ignore")  # get rid of annoying "color"
    # matplotlib warnings for right now

    # Parameters that might need to be changed:

    source_size = 1.5 * u.arcsec  # How large the source circle is in the plot.

    source_scale = rings[
        -1
    ]  # Use the outermost scale as the radius.  Make this somewhat larger than the
    # other ones, for reference

    # Read in galaxy parameters from the table.
    c1 = SkyCoord(t["Coords"][index], frame="icrs")
    pgc = t["PGC"][index]
    name = t["Name1"][index]
    major = t["SemiMajorAxis"][index] * u.arcmin
    minor = t["SemiMinorAxis"][index] * u.arcmin
    pa = t["PositionAngle"][index] * u.deg

    rlimit = major * source_scale
    radouter = rlimit * 2.1  # The size of the image displayed, slightly larger than 2.
    nrings = rings.size

    # Now, create a list of sources near the galaxy from the detect table.

    catalog = SkyCoord(tdetect["ra"], tdetect["dec"], unit=(u.deg, u.deg))

    d2d = c1.separation(catalog)  # Find all of separations
    catalogmsk = d2d < rlimit  # Create Boolean mask of objects inside of this
    near_sources = catalog[catalogmsk]

    # Create title information on ellipse for the plot
    s1 = str(pgc) + " " + str(name) + " " + t["Coords"][index] + "\n"
    s2 = "({0:0.03f} {1:0.03f} PA: {2:0.01f})".format(major, minor, pa)

    s3 = s1 + s2

    # Create a one arcminute scale bar on the East edge of the image, pointing North-South
    temp1 = c1.ra + (major * (source_scale * 0.925))
    temp2 = c1.dec - (0.5 * u.arcmin)
    start_sky = SkyCoord(temp1, temp2)
    temp3 = c1.ra + (major * (source_scale * 0.925))
    temp4 = c1.dec + (0.5 * u.arcmin)
    end_sky = SkyCoord(temp3, temp4)
    scale_bar = LineSkyRegion(start=start_sky, end=end_sky)

    # Now, create the ellipse regions list to be overplotted.  Loop over all but the last ring.
    regionList = []
    for i in range(0, nrings - 1):
        scale = rings[i]
        ellipse = create_ellreg(t, index, d25scale=scale)
        regionList.append(ellipse)

    # Create the outerRing circle, which shows the radius that we searched to.
    outerRing = []
    tempreg = CircleSkyRegion(center=c1, radius=rlimit)
    outerRing.append(tempreg)

    # Create a region from the sources we want to overplot
    sourceCircles = []
    for coord in near_sources:
        tempreg = CircleSkyRegion(center=coord, radius=source_size)
        sourceCircles.append(tempreg)

    # Find the total number of sources within each ring
    out_rings = num_rings(t, index, tdetect, rings)

    s4 = "RA - Object: {0:d} ".format(index) + str(t["Notes"][index]) + "\n"
    s5 = "R limit: {0:0.03f}".format(rlimit) + "\n"
    s6 = "Ring Counts: "
    for val in out_rings:
        s6 = s6 + "{0:d} ".format(np.int(val))
    s6 = s6 + "\n"
    s7 = s4 + s5 + s6

    try:
        imglist = SkyView.get_images(position=c1, radius=radouter, survey=image_survey)
    # If SDSS g is not found, try using the DSS2 image
    except:
        imglist = SkyView.get_images(position=c1, radius=radouter, survey="DSS2 Blue")

    # the returned value is a list of images, but there is only one
    img = imglist[0]

    # 'img' is now a fits.HDUList object; the 0th entry is the image
    mywcs = wcs.WCS(img[0].header)

    # Begin plotting here
    fig = plt.figure(figsize=(8.0, 8.0))
    ax = fig.add_subplot(111)

    ax.set_xlabel(s7)
    ax.set_ylabel("DEC")
    ax.set_title(s3)

    # Plot the Skyview fits file
    ax.imshow(
        img[0].data,
        cmap="gray_r",
        interpolation="none",
        origin="lower",
        norm=pl.matplotlib.colors.LogNorm(),
    )

    # overplot the one arcminute scale bar
    pixel_region = scale_bar.to_pixel(mywcs)
    pixel_region.plot(ax=ax, color="black", linewidth=4)

    # overplot the RC3 ellipses that were chosen
    for reg in regionList:
        pixel_region = reg.to_pixel(mywcs)
        pixel_region.plot(ax=ax, color="blue", linewidth=2)

    # overplot source position as purple circle
    for reg in sourceCircles:
        pixel_region = reg.to_pixel(mywcs)
        pixel_region.plot(
            ax=ax, edgecolor="purple", facecolor="purple", linewidth=5, fill=True
        )

    # overplot search radius
    for reg in outerRing:
        pixel_region = reg.to_pixel(mywcs)
        pixel_region.plot(ax=ax, color="green", linewidth=4)

    # Show the output
    if save_file:
        fig.savefig(outname, dpi=400)

    if show_notebook:
        plt.show()

    plt.close()


def check_all_large_gal(tboth, tdetect, scale, verbose=False):

    # Set up output arrays
    ndetect = len(tdetect)
    outTest = np.zeros(ndetect, dtype=bool)  # Set all to False intially
    outRedshift = np.zeros(ndetect)
    tempNames = []
    for i in range(0, ndetect):
        tempNames.append("NoGalaxy")

    # Safety margin - this can probably be reduced further, once we check things
    margin = 1.1
    nprint = 50  # Number to print progress

    # Create skyCoord of detections
    catalog = SkyCoord(tdetect["ra"], tdetect["dec"], unit=(u.deg, u.deg))

    # Loop over each galaxy.  Create a subset of objects that are close to the galaxy, and then
    # see if they fall within the ellipse. Close is defined as a fraction larger than the semi-major axis.

    ngals = len(tboth)
    for i in range(0, ngals):
        if verbose:
            if (i % nprint) == 0:
                print("Checking Galaxy: ", i, tboth["PGC"][i], tboth["Name1"][i])
        c1 = SkyCoord(tboth["Coords"][i], frame="icrs")
        rlimit = tboth["SemiMajorAxis"][i] * scale * margin * u.arcmin

        d2d = c1.separation(catalog)  # Find all of separations to the objects
        catalogmsk = (
            d2d < rlimit
        )  # Create Boolean mask of objects inside of this rlimit value
        if (
            len(catalogmsk) == 0
        ):  # If there are no close sources, go on to the next galaxy
            break
        else:  # Otherwise, loop over the sources, and look at the ellipses
            near_sources = catalog[catalogmsk]
            id_near_sources = tdetect["detectid"][catalogmsk]
            for j in range(0, len(near_sources)):
                coords = near_sources[j]
                idtest = id_near_sources[j]
                isclose, name, zgal = onegal_flag_from_coords(
                    coords, tboth, i, d25scale=scale
                )

                if isclose:
                    ind2 = np.int(np.where(tdetect["detectid"] == idtest)[0])
                    outTest[ind2] = True
                    tempNames[ind2] = name
                    outRedshift[ind2] = zgal

    outNames = np.array(tempNames, dtype=object)

    return (outTest, outNames, outRedshift)


# Test cases here- change to True, if you want to run the tests.
test1 = False
test2 = False  # test2 - Check to see if RC3 catalog loaded properly, from config location - not yet
test3 = False  # test3 - load in the detections file
test4 = False  # test4 - double check that the make wcs task worked
test5 = False  # test5 - create an elliptical region from the galaxy table, and make sure it makes sense
test6 = False  # test6 - test to see if we find objects near a specific galaxy, using the coordinates alone.
test7 = False  # test7 - Look over the entire catalog to see if we find objects near a galaxy, using the coordinates alone.
test8 = False  # test8 - Plot a single source, overlaid on the ellipsoidal rings of a specific galaxy
test9 = False  # test9 - Plot all the HETDEX sources around a particular galaxy
test10 = False  # test10 - search the entire detect list against the entire galaxy list

if test2:
    print(
        "\nTest 2 - Checking to see if RC3 catalog loaded from configuration location:"
    )
    print(70 * "-")
    tboth = read_rc3_tables()
    print(tboth.info)
    print("\nFirst and last entries: \n")
    print(tboth[0])
    print()
    print(tboth[-1])

if test4:
    print("\nTest 4 - create dummy WCS, and verify it:")
    print(70 * "-")

    # Example - the location of NGC 4707
    c1 = SkyCoord("12h48m22.87s +51d09m52.9s", frame="icrs")

    mywcs = create_dummy_wcs(c1)
    print(mywcs)
    print("For coordinate: ", c1.to_string("hmsdms"))
    pixcoord = PixCoord.from_sky(c1, mywcs, origin=1)
    print(pixcoord)
    print("Two arcsecond offsets in the cardinal directions")

    # create offsets, and make sure they are what we would expect
    sep = 2.0 * u.arcsec
    pas = np.float64([0.0, 180.0, 90.0, 270.0]) * u.deg

    for pa in pas:
        ctemp = c1.copy()
        skycoord = ctemp.directional_offset_by(pa, sep)
        print(skycoord.to_string("hmsdms"))
        print(skycoord.separation(ctemp).arcsec)
        pixcoord = PixCoord.from_sky(skycoord, mywcs, origin=1)
        print(pixcoord)
        print()


if test5:
    print(
        "\nTest 5 - create a random elliptical region from the RC3 galaxy table, and print it out:"
    )
    print(70 * "-")
    print("Loading RC3 catalog:")
    tboth = read_rc3_tables()

    ntot = len(tboth)
    index = np.random.randint(ntot)
    print("For random index: ", index)
    print(tboth["PGC"][index])
    ellipse1 = create_ellreg(tboth, index, d25scale=1.0)

    ellipse2 = create_ellreg(tboth, index, d25scale=2.0)
    print(ellipse1)
    print("Scaled up by a factor of 2:")
    print(ellipse2)


if test6:
    print("\nTest 6 - Search for sources by coordinates only for a single galaxy:")
    print(70 * "-")

    tboth = read_rc3_tables()

    tdetect = read_detect()

    # Find 200 objects in the tdetect list that are closest to NGC 4707
    nclose = 200
    index = np.int(np.where(tboth["Name1"] == "NGC 4707")[0])
    c1 = SkyCoord(tboth["Coords"][index], frame="icrs")

    # Read in coordinates from detect table, and sort by order of distance
    c2 = SkyCoord(tdetect["ra"], tdetect["dec"], unit=(u.deg, u.deg))
    d2d = c1.separation(c2)  # Find angular separations
    ind1 = np.argsort(d2d)  # sort by order of distance

    numd25_1 = np.zeros(nclose)
    numd25_2 = np.zeros(nclose)

    for i in range(0, nclose):
        coords = c2[ind1[i]]
        idstring = tdetect["detectid"][ind1[i]]
        dist = np.around(d2d[ind1[i]].arcmin, decimals=2)
        isclose, name, zgal = onegal_flag_from_coords(
            coords, tboth, index, d25scale=1.0
        )
        isclose2, name2, zgal2 = onegal_flag_from_coords(
            coords, tboth, index, d25scale=2.0
        )
        print(
            idstring,
            coords.to_string(),
            dist,
            isclose,
            name,
            zgal,
            isclose2,
            name2,
            zgal2,
        )
        numd25_1[i] = isclose
        numd25_2[i] = isclose2

    print(
        "Out of {0:d} objects closest, the number within D25, 2D25: {1:0.01f} {2:0.01f}".format(
            nclose, np.sum(numd25_1), np.sum(numd25_2)
        )
    )


if test7:
    print("\nTest 7 - Search for sources by coordinates only for a range of galaxies:")
    print(70 * "-")

    print("Loading RC3 catalog:")
    tboth = read_rc3_tables()
    
    print("Loading HDR catalog")
    tdetect = read_detect()

    ntot = len(tboth)

    # find a random galaxy that has objects near it.  Loop over random indices until we find one that has sources nearby

    while True:
        index = np.random.randint(ntot)
        c1 = SkyCoord(tboth["Coords"][index], frame="icrs")
        # Read in coordinates from detect table, and sort by order of distance
        c2 = SkyCoord(tdetect["ra"], tdetect["dec"], unit=(u.deg, u.deg))
        d2d = c1.separation(c2)  # Find angular separations
        if np.min(d2d).arcmin < 1.0:
            break

    # Find 20 objects in the tdetect list that are closest to a particular galaxy
    nclose = 20

    ind1 = np.argsort(d2d)  # sort by order of distance
    print("For random index: ", index)
    print(tboth["PGC"][index])
    numd25_1 = np.zeros(nclose)
    numd25_2 = np.zeros(nclose)

    for i in range(0, nclose):
        coords = c2[ind1[i]]
        idstring = tdetect["detectid"][ind1[i]]
        dist = np.around(d2d[ind1[i]].arcmin, decimals=2)
        isclose, name, zgal = gals_flag_from_coords(
            coords, tboth, d25scale=1.0, nmatches=3
        )
        isclose2, name2, zgal2 = gals_flag_from_coords(
            coords, tboth, d25scale=2.0, nmatches=3
        )
        print(
            idstring,
            coords.to_string(),
            dist,
            isclose,
            name,
            zgal,
            isclose2,
            name2,
            zgal2,
        )
        numd25_1[i] = isclose
        numd25_2[i] = isclose2

    print(
        "Out of {0:d} objects closest, the number within D25, 2D25: {1:0.01f} {2:0.01f}".format(
            nclose, np.sum(numd25_1), np.sum(numd25_2)
        )
    )


if test8:
    print("\nTest 8: Plot a single source near an RC3 galaxy: ")
    print(70 * "-")

    print("Loading RC3 catalog:")
    tboth = read_rc3_tables()

    print("Loading HDR catalog")
    tdetect = read_detect()

    # examples of near and far sources
    detectid_near = 2100710607  # source very near NGC 4707 nucleus
    detectid_far = 2100176247  # source farther away from NGC 4707 nucleus

    ind1 = np.int(np.where(tboth["Name1"] == "NGC 4707")[0])
    ind2 = np.int(np.where(tdetect["detectid"] == detectid_near)[0])
    coords = SkyCoord(tdetect["ra"][ind2], tdetect["dec"][ind2], unit=(u.deg, u.deg))

    # Define ellipsoidal rings, in units of D25
    rings = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.5])

    show_ellipse_ring_source(
        tboth, ind1, coords, rings, "", save_file=False, show_notebook=True
    )

    ind3 = np.int(np.where(tdetect["detectid"] == detectid_far)[0])
    coords = SkyCoord(tdetect["ra"][ind3], tdetect["dec"][ind3], unit=(u.deg, u.deg))

    show_ellipse_ring_source(
        tboth, ind1, coords, rings, "", save_file=False, show_notebook=True
    )


if test9:
    print("\nTest 9: Overplot all tdetect sources near a single RC3 galaxy:")
    print(70 * "-")

    print("Loading RC3 catalog:")
    tboth = read_rc3_tables()

    print("Loading HDR catalog")
    tdetect = read_detect()

    # Example - plot all of the tdetect sources located around NGC 4707:
    ind1 = np.int(np.where(tboth["Name1"] == "NGC 4707")[0])
    rings = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.5])
    show_ellipse_ring_detect(
        tboth, ind1, tdetect, rings, "", save_file=False, show_notebook=True
    )


if test10:
    print(
        "\nTest 10: Search the entire tdetect list, and find out how many sources would be within the masks:"
    )
    print("This test takes a significant amount of time, (5 minutes) per scale factor!")
    print(70 * "-")

    tboth = read_rc3_tables()
    tdetect = read_detect()

    d25scales = np.array([1.0, 1.25, 1.5, 1.75, 2.0])  # Different values of d25 scales
    ntot = len(tdetect)

    # If you need to, create a smaller sub-sample of the catalog and call it
    # "tsearch"
    # npieces = 4
    # tsearch = f[:ntot //npieces]

    for scale in d25scales:
        startTime = time.time()
        print("Beginning scale search: ", scale)
        outTest, outNames, outRedshifts = check_all_large_gal(
            tboth, tdetect, scale, verbose=False
        )

        print("\nFive random sources from tdetect:")
        for j in range(0, 5):
            index = np.random.randint(ntot)
            print(
                tdetect["detectid"][index],
                outTest[index],
                outNames[index],
                outRedshifts[index],
            )

        print("\nFive random matches from tdetect:")
        for j in range(0, 5):
            while True:
                index = np.random.randint(ntot)
                if outTest[index]:
                    break
            print(
                tdetect["detectid"][index],
                outTest[index],
                outNames[index],
                outRedshifts[index],
            )

        nfound = np.sum(outTest)
        print()
        frac = np.around((nfound / ntot), decimals=4)
        print("Scale, nfound, fraction: ", scale, nfound, frac)
        endTime = time.time()
        deltaTime = np.around((endTime - startTime), decimals=2)
        print("\nTime elapsed for entire run (s): ", deltaTime)

print("\nTests are complete.")
