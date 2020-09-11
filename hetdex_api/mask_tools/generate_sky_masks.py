"""

Generate a masks for the data and the bad amplifiers in Mangle vertices format. 
Assume dither positions are fiducial. Also, usually extend IFU corners out
to 30 arcseconds to account for possible detections on the IFU edge. 
Word of warning: this makes the masks harder to use to compute survey area.

.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>

"""
try:
    # Python 3
    from urllib.request import urlopen
    from urllib.error import HTTPError
except ImportError:
    # Python 2
    from urllib2 import urlopen, HTTPError

from os.path import isfile
from six import iteritems
import re
from numpy import array
import tables as tb
from astropy.table import Table
from pyhetdex.het.fplane import FPlane, NoIFUError
from pyhetdex.coordinates.tangent_projection import TangentPlane
from hetdex_api.mask_tools.amplifier_positions import amp_corners, swapped_around_amps


def get_fplane(filename, datestr='', actpos=False, full=True):
    """
    Function from 

    https://fplaneserver.readthedocs.io/en/latest/user.html#fplane-retrieval
    Accessed 08/10/2019

    """

    url = 'https://luna.mpe.mpg.de/fplane/' + datestr

    if actpos:
        url += '?actual_pos=1'
    else:
        url += '?actual_pos=0'

    if full:
        url += '&full_fplane=1'
    else:
        url += '&full_fplane=0'

    try:
        resp = urlopen(url)
    except HTTPError as e:
        print(url)
        raise Exception(' Failed to retrieve fplane file, server '
                        'responded with %d %s' % (e.getcode(), e.reason))

    with open(filename, 'w') as f:
        f.write(resp.read().decode())


def generate_ifu_mask(output_fn, survey_hdf, badshots_file, ramin, ramax, decmin, decmax, specific_shot=None,
                      xsize=25.0, ysize=25.0, other_cuts={}, specific_field=None):
    """
    Generate a mask of IFU corners from the survey HDF 

    Parameters
    ----------
    output_fn : str
        a file to output the ra/dec corners to
    survey_hdf : str
        path to the survey HDF
    badshots_file : str
        path to the bad shots file
    ramin, ramax, decmin, decmax : float
        restrict the mask to a subregion
    specific_shot : str (Optional)
        overides the ra and dec range
        and instead only outputs a bad
        mask only for the given shotid
    xsize, ysize : float
        half of the size of the IFU in x and y
        in arcseconds. Vertices are produced
        a +/-xsize and +/-ysize. Optional,
        default xsize=ysize=25.0
    other_cuts : dict
        dictionary of shot property and a
        2 element list of minimum and maximum
        allowed value
    """
    # Read in the survey file
    survey_hdf = tb.open_file(survey_hdf)
       
    # Read bad shots file
    table_bad_shots = Table.read(badshots_file, format="ascii", names = ["shotid"])
    bad_shots = array(table_bad_shots["shotid"])

    if specific_shot is not None:
        survey_ttable = survey_hdf.root.Survey.read_where('shotid == specific_shot')
    elif specific_field is not None:
        survey_ttable = survey_hdf.root.Survey.read_where('field == specific_field')
    else:   
        query = '(ra < ramax) & (ra > ramin) & (dec < decmax) & (dec > decmin)'

        for param, lims in iteritems(other_cuts):
            query += ' & ({:s} > {:f}) & ({:s} < {:f})'.format(param, lims[0], param, lims[1])

        print(query) 
        survey_ttable = survey_hdf.root.Survey.read_where(query)

    # Loop over the datevshots and see if there are bad amps
    polys = []
    shids = []
    for line in survey_ttable:
        print(line)
        # skip bad shots
        if line["shotid"] in bad_shots:
            continue

        date = line["date"]

        rot = 360.0 - (line["pa"] + 90.)
        tp = TangentPlane(line["ra"], line["dec"], rot)
        fplane_fn = "fplanes/{:d}_fplane.txt".format(date)
   
        if not isfile(fplane_fn):
           get_fplane(fplane_fn, datestr=str(date))
           fplane = FPlane(fplane_fn)
        else:
           fplane = FPlane(fplane_fn) 

        rect = [[-1.0*xsize, -1.0*xsize, xsize, xsize],
                [-1.0*ysize, ysize, ysize, -1.0*ysize]]
    
        for ifu in fplane.ifus: 
            x = array(rect[0]) + ifu.y  
            y = array(rect[1]) + ifu.x
            ra, dec = tp.xy2raDec(x, y)
            polys.append([ra[0], dec[0], ra[1], dec[1], 
                          ra[2], dec[2], ra[3], dec[3]])  
            shids.append(line["shotid"])
 
    # Should now have a list of polygons to output
    with open(output_fn, "w") as fp:
        for poly, shid in zip(polys, shids):
            fp.write("{:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:d}\n".format(*poly, shid))




def generate_bad_amp_mask_fits(output_fn, survey_hdf, badamps_fits, ramin, ramax, decmin, 
                               decmax, specific_shot = None):
    """
    Generate a Mangle-compatible list of ra/dec pairs corresponding
    to the corners of the bad amplifiers on the sky. The amplifiers
    are split into squares and rectangles to better follow their 
    shape.

    Parameters
    ----------
    output_fn : str
        a file to output the ra/dec corners to
    survey_hdf : str
        path to the survey HDF
    badamps_file : str
        path to the bad amps file
    ramin, ramax, decmin, decmax : float
        restrict the mask to a subregion
    specific_shot : str (Optional)
        overides the ra and dec range
        and instead only outputs a bad
        mask only for te given shotid
    """
    # Read in the survey file
    survey_hdf = tb.open_file(survey_hdf)
   
    if specific_shot is not None:
        # 20190209027
        survey_ttable = survey_hdf.root.Survey.read_where('shotid == specific_shot')
    else:    
        survey_ttable = survey_hdf.root.Survey.read_where('(ra < ramax) & (ra > ramin) & (dec < decmax) & (dec > decmin)')

    # Read in the bad amps
    table_bad_amps = Table.read(badamps_fits)
    table_bad_amps = table_bad_amps[table_bad_amps["flag"] == 0]
    pattern = re.compile("multi_[0-9]{3}_([0-9]{3})_[0-9]{3}_([RLU]{2})") 
    table_bad_amps["AMP"] = [pattern.findall(x)[0][1] for x in table_bad_amps["multiframe"]]
    table_bad_amps["IFUSLOT"] = [pattern.findall(x)[0][0] for x in table_bad_amps["multiframe"]]

    # Loop over the datevshots and see if there are bad amps
    polys = []
    amps = []
    for line in survey_ttable:
        bad_amps_here = table_bad_amps[table_bad_amps["shotid"] == line["shotid"]]

        # If any, grab the focal plane and generate 
        # a tangent plane for the astrometry
        if len(bad_amps_here) > 0:
            #print("{:d} has bad amps. Adding to mask".format(line["shotid"]))
    
            date = line["date"]
            fplane_fn = "fplanes/{:d}_fplane.txt".format(date)

            if not isfile(fplane_fn):
               get_fplane(fplane_fn, datestr=str(date))
            else:
               fplane = FPlane(fplane_fn) 
    

            rot = 360.0 - (line["pa"] + 90.)
            tp = TangentPlane(line["ra"], line["dec"], rot)
                   
            for bad_amp in bad_amps_here:
                try:
                    ifu = fplane.by_ifuslot("{:s}".format(bad_amp["IFUSLOT"]))
                except NoIFUError:
                    print("Warning. IFU {:s} not found for dateobs {:d}".format(bad_amp["IFUSLOT"], line["shotid"]))
                    continue
    
                # Check if the amps in this IFU are swapped around
                ampkey = "{:s}{:s}".format(bad_amp["IFUSLOT"], bad_amp["AMP"])
                if ampkey in swapped_around_amps:
                    amp = swapped_around_amps[ampkey]
                else:
                    amp = bad_amp["AMP"]
                 
                # coordinates of amplifier for default dither and IFU cen
                rects_to_mask = amp_corners[amp]
    
                for rect in rects_to_mask:
                    # Flip is correct
                    x = array(rect[0]) + ifu.y  
                    y = array(rect[1]) + ifu.x
                    ra, dec = tp.xy2raDec(x, y)
                    polys.append([ra[0], dec[0], ra[1], dec[1], 
                                  ra[2], dec[2], ra[3], dec[3]])
                   
                    amps.append(bad_amp["AMP"])
    
    # Should now have a list of polygons to output
    with open(output_fn, "w") as fp:
        for poly, amp in zip(polys, amps):
            fp.write("{:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:s}\n".format(*poly, amp))

   

def generate_bad_amp_mask(output_fn, survey_hdf, badamps_file, badshots_file, 
                          ramin, ramax, decmin, decmax, specific_shot=None):

    """
    Generate a Mangle-compatible list of ra/dec pairs corresponding
    to the corners of the bad amplifiers on the sky. The amplifiers
    are split into squares and rectangles to better follow their 
    shape.

    Parameters
    ----------
    output_fn : str
        a file to output the ra/dec corners to
    survey_hdf : str
        path to the survey HDF
    badamps_file : str
        path to the bad amps file
    ramin, ramax, decmin, decmax : float
        restrict the mask to a subregion
    specific_shot : str (Optional)
        overides the ra and dec range
        and instead only outputs a bad
        mask only for te given shotid

    """

    # Read in the survey file
    survey_hdf = tb.open_file(survey_hdf)
   
    if specific_shot is not None:
        # 20190209027
        survey_ttable = survey_hdf.root.Survey.read_where('shotid == specific_shot')
    else:    
        survey_ttable = survey_hdf.root.Survey.read_where('(ra < ramax) & (ra > ramin) & (dec < decmax) & (dec > decmin)')

    # Read in the bad amps
    table_bad_amps = Table.read(badamps_file, format="ascii", names = ["IFUSLOT", "AMP", "multiframe", 
                                                                       "start", "end"])
    # Read bad shots file
    table_bad_shots = Table.read(badshots_file, format="ascii", names = ["shotid"])
    bad_shots = array(table_bad_shots["shotid"])

    # Loop over the datevshots and see if there are bad amps
    polys = []
    amps = []
    for line in survey_ttable:

        date = line["date"]
        fplane_fn = "fplanes/{:d}_fplane.txt".format(date)

        if not isfile(fplane_fn):
           get_fplane(fplane_fn, datestr=str(date))
        else:
           fplane = FPlane(fplane_fn) 
 
        if line["shotid"] in bad_shots:
            # if shot bad mask all IFUS
            print("Masking whole bad shot found {:d}".format(line["shotid"]))
            bad_amps_here = Table()
            bad_amps_here["IFUSLOT"] = [int(x) for x in fplane.ifuslots]
            bad_amps_here["AMP"] = "AA"
        else:
            # otherwise just mask bad amps
            sel = (table_bad_amps["start"] <= date) & (table_bad_amps["end"] >= date) 
            bad_amps_here = table_bad_amps[sel]
    
        # If any, grab the focal plane and generate 
        # a tangent plane for the astrometry
        if len(bad_amps_here) > 0:
           print("{:d} has bad amps. Adding to mask".format(line["shotid"]))
    
           rot = 360.0 - (line["pa"] + 90.)
           tp = TangentPlane(line["ra"], line["dec"], rot)
                  
           for bad_amp in bad_amps_here:
               try:
                   ifu = fplane.by_ifuslot("{:03d}".format(bad_amp["IFUSLOT"]))
               except NoIFUError:
                   print("Warning. IFU {:d} not found for dateobs {:d}".format(bad_amp["IFUSLOT"], line["shotid"]))
                   continue
    
               if bad_amp["AMP"] == "AA":
                   # whole IFU with a little extra border
                   rects_to_mask = [[[-30.0, -30.0, 30.0, 30.0],
                                     [-30.0, 30.0, 30.0, -30.0]]]
               else:
                   # Check if the amps in this IFU are swapped around
                   ampkey = "{:03d}{:s}".format(bad_amp["IFUSLOT"], bad_amp["AMP"])
                   if ampkey in swapped_around_amps:
                       amp = swapped_around_amps[ampkey]
                   else:
                       amp = bad_amp["AMP"]
                    
                   # coordinates of amplifier for default dither and IFU cen
                   rects_to_mask = amp_corners[amp]
    
               for rect in rects_to_mask:
                   # Flip is correct
                   x = array(rect[0]) + ifu.y  
                   y = array(rect[1]) + ifu.x
                   ra, dec = tp.xy2raDec(x, y)
                   polys.append([ra[0], dec[0], ra[1], dec[1], 
                                 ra[2], dec[2], ra[3], dec[3]])
                  
                   amps.append(bad_amp["AMP"])
    
    # Should now have a list of polygons to output
    with open(output_fn, "w") as fp:
        for poly, amp in zip(polys, amps):
            fp.write("{:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:7.6f} {:s}\n".format(*poly, amp))



if __name__ == "__main__":

    # The HDR2.1 bias region
    # change this config
    survey_hdf = "survey/survey_hdr2.1.h5"
    badamps_file = "hdr2.1_issues/badamps.list"
    badshots_file = "hdr2.1_issues/badshots.list"
    badamps_fits="survey/amp_flag.fits"

    # whole sky
    #generate_ifu_mask("full_hdr2pt1_30as.vert", survey_hdf, badshots_file, 0.0, 360.0, -90.0, 90.0, 
    #                  xsize=30.0, ysize=30.0

    #generate_ifu_mask("fall_hdr2pt1.vert", survey_hdf, badshots_file, 0.0, 360.0, -90.0, 90.0, 
    #                  xsize=30.0, ysize=30.0, specific_field="dex-fall")

    #generate_bad_amp_mask("hdr2.1_bas_mask.vert", survey_hdf, badamps_file, badshots_file, 0, 360.0, -90, 90)
    generate_bad_amp_mask_fits("hdr2.1_bad_mask_fits.vert", survey_hdf, badamps_fits, 0, 360.0, -90, 90)

    #(response_4540 > 0.08) & (fwhm_moffat < 2.6)
    #generate_ifu_mask("hdr2.1_ifus.vert", survey_hdf, badshots_file, 195, 215, 50, 52, 
    #                   xsize=26, ysize=26, other_cuts={'fwhm_moffat' : [0.0, 2.6], 'response_4540' : [0.095, 999]})
