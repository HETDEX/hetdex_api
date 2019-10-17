# -*- coding: utf-8 -*-
"""
Created on October 8 2019
@author: Erin Mentuch Cooper

Extracts 1D spectrum at specified RA/DEC.

Provide either an ID/RA/DEC as input or a file
with a list of ID, RA, DEC. Can be in the format
of an astropy table provided columns are appropriately 
marked 'ID','ra', 'dec'. 

Option can be to provide a single
 shotid/datevobs to work on for batch processing, 
or can let program figure out which shots
to extract spectra on.


For a single coordinate, you can search and extract all shots in HDR1

python3 get_spec.py --ra 150.02548 --dec 2.087987 --ID cosmos_LAE --outfile cosmos_LAE

To speed up, you can use the --multiprocess option. This works for interactive use
on a node either through a jupyter notebook or an interactive node called by idev

python3 get_spec.py --multiprocess --ra 150.02548 --dec 2.087987 --ID cosmos_LAE --outfile cosmos_LAE  

For working on a single shot only you can use the -s option. 
It is recommended to batch process each 
SHOTID individually and provide the same source list in each call.

python3 get_spec.py -ra 8.86535 -dec 0.59352  -s 20190104008

This can also work on an input file. It should either have the format of
3 columns = ID/RA/DEC or be an astropy table with columns 
labeled 'ID', 'ra', 'dec'

python3 get_spec.py --multiprocess True -i '3dhst_input.cat' -o '3dhst' 

This for example can be broken into batch jobs:

python3 get_spec.py -i '3dhst_input.cat'  


"""

import sys
import argparse as ap
import os
import os.path as op
import glob

import numpy as np
import pickle
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")

from input_utils import setup_logging


from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, join

from hetdex_api.extract import Extract
from hetdex_api.survey import Survey
from hetdex_api.shot import *

from copy import deepcopy
from collections.abc import Mapping

from multiprocessing import Pool, Process, Manager
import time

def merge(dict1, dict2):
    """ Return a new dictionary by merging two dictionaries recursively. """
    result = deepcopy(dict1)

    for key, value in dict2.items():
        if isinstance(value, Mapping):
            result[key] = merge(result.get(key, {}), value)
        else:
            result[key] = deepcopy(dict2[key])

    return result

def get_source_spectra(shotid, args):
    E = Extract()
    source_dict = {}
    if len(args.matched_sources[shotid]) > 0:
        args.log.info('Working on shot: %s' % shotid)
        fwhm = args.survey.fwhm_moffat[ args.survey.shotid == shotid ][0]
        moffat = E.moffat_psf(fwhm, 10.5, 0.25)
        E.load_shot(shotid)
        
        for ind in args.matched_sources[shotid]:
            if np.size(args.coords) > 1:
                info_result = E.get_fiberinfo_for_coord(args.coords[ind], radius=7.)
            else:
                info_result = E.get_fiberinfo_for_coord(args.coords, radius=7.)
            if info_result is not None:
                if np.size(args.ID) > 1:
                    args.log.info('Extracting %s' % args.ID[ind])
                else:
                    args.log.info('Extracting %s' % args.ID)
                
                ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result
                weights = E.build_weights(xc, yc, ifux, ifuy, moffat)
                result = E.get_spectrum(data, error, mask, weights)
                spectrum_aper, spectrum_aper_error = [res for res in result]
                if np.size(args.ID) > 1:
                    if args.ID[ind] in source_dict:
                        source_dict[args.ID[ind]][shotid] = [spectrum_aper, spectrum_aper_error, weights.sum(axis=0)]
                    else:
                        source_dict[args.ID[ind]] = dict()  
                        source_dict[args.ID[ind]][shotid] = [spectrum_aper, spectrum_aper_error, weights.sum(axis=0)] 
                else:
                    if args.ID in source_dict:
                        source_dict[args.ID][shotid] = [spectrum_aper, spectrum_aper_error, weights.sum(axis=0)]
                    else:
                        source_dict[args.ID] = dict()
                        source_dict[args.ID][shotid] = [spectrum_aper, spectrum_aper_error, weights.sum(axis=0)]

        E.fibers.close()
        return source_dict


def get_source_spectra_mp(source_dict, shotid, manager, args):
    E = Extract()

    if len(args.matched_sources[shotid]) > 0:
        args.log.info('Working on shot: %s' % shotid)
        fwhm = args.survey.fwhm_moffat[ args.survey.shotid == shotid ][0]
        moffat = E.moffat_psf(fwhm, 10.5, 0.25)
        E.load_shot(shotid)

        for ind in args.matched_sources[shotid]:
            if np.size(args.coords) > 1:
                info_result = E.get_fiberinfo_for_coord(args.coords[ind], radius=7.)
            else:
                info_result = E.get_fiberinfo_for_coord(args.coords, radius=7.)
            if info_result is not None:
                if np.size(args.ID) > 1:
                    args.log.info('Extracting %s' % args.ID[ind])
                else:
                    args.log.info('Extracting %s' % args.ID)

                ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result
                weights = E.build_weights(xc, yc, ifux, ifuy, moffat)
                result = E.get_spectrum(data, error, mask, weights)
                spectrum_aper, spectrum_aper_error = [res for res in result]
                if np.size(args.ID) > 1:
                    if args.ID[ind] in source_dict:
                        source_dict[args.ID[ind]][shotid] = [spectrum_aper, spectrum_aper_error, weights.sum(axis=0)]
                    else:
                        source_dict[args.ID[ind]] = manager.dict()
                        source_dict[args.ID[ind]][shotid] = [spectrum_aper, spectrum_aper_error, weights.sum(axis=0)]
                else:
                    if args.ID in source_dict:
                        source_dict[args.ID][shotid] = [spectrum_aper, spectrum_aper_error, weights.sum(axis=0)]
                    else:
                        source_dict[args.ID] = manager.dict()
                        source_dict[args.ID][shotid] = [spectrum_aper, spectrum_aper_error, weights.sum(axis=0)]

        E.fibers.close()
        
def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create rsp tmpXXX.datfiles.""",
                               add_help=True)

    parser.add_argument("-s", "--shot",
                        help='''ShotID, e.g., 20170321v009, YYYYMMDDvOBS''',
                        type=str, default=None)
    
    parser.add_argument("-i", "--infile",
                        help='''File with table of ID/RA/DEC''', default=None)

    parser.add_argument("-o", "--outfile",
                        help='''File to store pickled dictionary output''', default='output', type=str)

    parser.add_argument("-ra", "--ra",
                        help='''ra, e.g., right ascension in degrees''',
                        type=float, default=None)

    parser.add_argument("-dec", "--dec",
                        help='''ra, e.g., right ascension in degrees''',
                        type=float, default=None)
    
    parser.add_argument("-rad", "--rad",
                        help='''radius, e.g., aperture radius in arcsec''',
                        type=float, default=3.0)

    parser.add_argument("-id", "--ID",
                        help='''source name''',
                        type=str, default=None)
    
    parser.add_argument("--multiprocess", "-mp", help='''Multiprocessing Flag''', 
                        type=bool, default=False)

    parser.add_argument("--merge", '-merge', help='''Boolean trigger to merger are pickle files in cwd''',
                        default=False, type=bool)

    parser.add_argument("--mergepath",'-mpath', help='''Path to location of pickle files to merge''',
                        default=os.getcwd(), type=str)
 

    args = parser.parse_args(argv)
    args.log = setup_logging()

    if args.merge:
        all_source_dict = {}
        files = glob.glob(op.join(args.mergepath,'*.pkl'))
        args.log.info('Merging all pickle files in ' + args.mergepath)
        for file in files:
            file_dict = pickle.load( open( file, 'rb'))
            all_source_dict = merge( all_source_dict, file_dict )
            
        outfile = args.outfile+'.pkl'
        
        pickle.dump( all_source_dict, open( outfile, "wb" ) )
        args.log.info('Saved output file to '+ outfile)
        sys.exit('Exiting script')

    if args.infile:

        args.log.info('Loading External File')

        try: 
            table_in = Table.read(args.infile, format='ascii')
        except:
            table_in = Table.read(args.infile, format='no_header', names=['ID','ra','dec'])

        args.ID = table_in['ID']
        args.ra = table_in['ra']
        args.dec = table_in['dec']

    else:
        if args.ID == None:
            args.ID = 'DEX_' + str(args.ra).zfill(4)+'_'+str(args.dec).zfill(4)

        args.log.info('Extracting for ID: %s' % args.ID)

    args.coords = SkyCoord(args.ra*u.deg, args.dec*u.deg)

    args.survey = Survey('hdr1')

    # if args.shot exists, only select that shot

    if args.shot:
        try:
            sel_shot = args.survey.shotid == int(args.shot)
        except:
            sel_shot = args.survey.datevobs == str(args.shot)

        args.survey = args.survey[sel_shot]

    else:
        args.log.info('Searching through all shots')
            
    args.matched_sources = {}
    shots_of_interest = []

    count = 0

    # this radius applies to the inital shot search and requires a large aperture for the wide FOV of VIRUS
    max_sep = 11.0 * u.arcminute
    
    args.log.info('Finding shots of interest')
 
    for i, coord in enumerate(args.survey.coords):
        dist = args.coords.separation(coord)
        sep_constraint = dist < max_sep
        shotid = args.survey.shotid[i]
        idx = np.where(sep_constraint)[0]
        if np.size(idx) > 0:
            args.matched_sources[shotid] = idx
            count += np.size(idx)
            if len(idx) > 0:
                shots_of_interest.append(shotid)
        
    args.log.info('Number of shots of interest: %i' % len(shots_of_interest))
    args.log.info('Extracting %i sources' % count)

    if args.multiprocess:

        manager = Manager()
        Sources = manager.dict()
        
        start = time.time()
        njobs = np.size(shots_of_interest)

        ntasks = 32

        for i in np.arange(0, njobs, ntasks):
            jobs = [ Process(target=get_source_spectra_mp, args=(Sources, shotid, manager, args)) 
                     for shotid in shots_of_interest[i: np.minimum(njobs, i+ntasks)]
                 ]
            
            for j in jobs:
                j.start()
            for j in jobs:
                j.join()
    
        end = time.time()
        args.log.info( 'Extraction of sources completed in %.2f minutes.' % ((end-start)/60.))
        
        #convert manager dict to a regular dictionary object to pickle
        Source_dict = {}
        count = 0
        for ID in Sources.keys():
            Source_dict[ID] = {}
            for shotid in Sources[ID].keys():
                count += 1
                Source_dict[ID][shotid] = list( Sources[ID][shotid])

    else:
        Source_dict = {}
        for shotid in shots_of_interest:
            shot_source_dict = get_source_spectra(shotid, args)

            if shot_source_dict:
                Source_dict = merge( Source_dict, shot_source_dict)

    args.survey.close()
    
    outfile = args.outfile+'.pkl'
    pickle.dump( Source_dict, open( outfile, "wb" ) )
        
if __name__ == '__main__':
    main()
