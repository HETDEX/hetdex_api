# -*- coding: utf-8 -*-
"""
Created on March 29 2019
@author: Erin Mentuch Cooper

Extracts 2D spectral image for a specific detection

You can specify the size in pixels..more options to come

python get_2D_spec.py -det 1000572696 -dx 40 -dy 100

"""

from __future__ import print_function

import glob
import re

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import ZScaleInterval

from input_utils import setup_logging
from hetdex_api.shot import *
from hetdex_api.detections import *


def get_2Dimage(detectid_obj, detects, fibers, width=100, height=50):
    
    fiber_table = Table(detects.hdfile.root.Fibers.read_where("detectid == detectid_obj") )
    
    weight_total = np.sum(fiber_table['weight'])

    #open blank image to store summed image
    im_sum = np.zeros((height, width))

    for fiber_i in fiber_table:
        wave_obj = fiber_i['wavein']
        multiframe_obj = fiber_i['multiframe']
        expnum_obj = fiber_i['expnum']
        fibnum_obj = fiber_i['fibnum']
        weight = fiber_i['weight']/weight_total
        im_fib = fibers.get_fib_image2D(wave_obj, fibnum_obj, multiframe_obj, expnum_obj, width=width, height=height)
        im_sum += weight * im_fib
        
    return im_sum

def get_2Dimage_wave(detectid_obj, detects, fibers, width=100, height=50):
    fiber_table = Table(detects.hdfile.root.Fibers.read_where("detectid == detectid_obj") )
    maxfib = np.argmax(fiber_table['weight'])    
    wave_obj = fiber_table['wavein'][maxfib]
    multiframe_obj = fiber_table['multiframe'][maxfib]
    expnum_obj = fiber_table['expnum'][maxfib]
    fibnum_obj = fiber_table['fibnum'][maxfib]
    idx = np.where((fibers.fibidx == (fibnum_obj - 1) ) * (fibers.multiframe == multiframe_obj) * (fibers.expnum == expnum_obj))[0][0]
    im_fib = fibers.wavelength[idx]

    sel = np.where(im_fib >= wave_obj)[0][0]
    x1 = np.maximum(0, sel - int(width/2))
    x2 = np.minimum(1031, sel + int(width/2))
    
    x1_slice = np.minimum(0, width - (x2-x1)) 
    x2_slice = x2-x1

    im_wave = np.zeros(width)
    im_wave[x1_slice:x2_slice] = im_fib[x1:x2]
    
    return im_wave


def save_2Dimage(detectid_i, detects, fibers, width=100, height=20, path=os.getcwd()):
    #try:
    im_sum = get_2Dimage(detectid_i, detects, fibers, width=width, height=height)
    im_wave = get_2Dimage_wave(detectid_i, detects, fibers, width=width, height=height)
    im_out = np.vstack([im_wave, im_sum])
    zscale = ZScaleInterval(contrast=0.5,krej=2.5) 
    vmin, vmax = zscale.get_limits(values=im_out)
    plt.imshow(im_out[1:-1,:],vmin=vmin, vmax=vmax,extent=[im_out[0,0], im_out[0,-1], -int(height/2.), int(height/2.)], origin="lower",cmap=plt.get_cmap('gray'),interpolation="none")
    plt.savefig(op.join(path, 'imsum2D_' + str(detectid_i) +  '.png'))
    np.savetxt(op.join(path, 'imsum2D_' + str(detectid_i) + '.txt'), im_out)
    #except:
     #   print("failed to create sum, dims prob don't match")

def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create rsp tmpXXX.datfiles.""",
                               add_help=True)

    parser.add_argument("-s", "--datevobs",
                        help='''ShotID, e.g., 20170321v009, YYYYMMDDvOBS''',
                        type=str, default=None)
    
    parser.add_argument("-dets" "-detectlist",
                        help='''filelist of detectids''', default=None)

    parser.add_argument("-ra", "--ra",
                        help='''ra, e.g., right ascension in degrees''',
                        type=float, default=None)

    parser.add_argument("-dec", "--dec",
                        help='''ra, e.g., right ascension in degrees''',
                        type=float, default=None)

    parser.add_argument("-rad", "--rad",
                        help='''radius, e.g., aperture radius in arcsec''',
                        type=float, default=3.0)

    parser.add_argument("-w", "--wave",
                        help='''wavelength in AA''',
                        type=float, default=None)

    parser.add_argument("-dw", "--dwave",
                        help='''delta wavelength in AA''',
                        type=float, default=50.0)

    parser.add_argument("-dx", help='''Width of image in wavelength dim in pix''')

    parser.add_argument("-dy", help='''Height of fiber image in pixels''')

    parser.add_argument("-p", "--path", help='''Path to save output''',
                        default=os.getcwd(), type=str)

    parser.add_argument("-i", "--infile",
                        help='''File with table of ID/RA/DEC''', default=None)



    args = parser.parse_args(argv)
    args.log = setup_logging()

    # initiate Fibers and Detections Classes
    print(args)

    shotid_i = args.datevobs

    detects = Detections('hdr1').refine()

    print("opening shot: "+str(shotid_i))
    fibers = Fibers(args.datevobs)

    
    if args.infile:
        try:
            catalog = Table.read(args.infile, format=ascii)
        except:
            catalog = Table.read('/work/05350/ecooper/hdr1/catalogs/hdr1_sngt6pt5_for.tab', format='ascii')
            
        selcat = (catalog['shotid'] == int(shotid_i))     
        detectlist = np.array(catalog['detectid'][selcat])

    elif args.dets:
        detectlist = np.loadtxt(args.dets, dtype=int)


    for detectid_i in detectlist:
        print(detectid_i)
        save_2Dimage(detectid_i, detects, fibers, width=100, height=20, path=args.path)
            
    if args.ra and args.dec : 

        obj_coords = SkyCoord(args.ra * u.deg, args.dec * u.deg, frame='icrs')
        idx = fibers.query_region_idx(obj_coords, radius=(args.rad/3600.))

        output = Table()
        output['ra'] = fibers.coords.ra[idx]*u.deg
        output['dec'] = fibers.coords.dec[idx]*u.deg
        filenames = []

        fileidx = 101
        for i in idx:
            filename = 'tmp' + str(fileidx) + '.dat'
            filenames.append(filename)
            save_rsp_spectrum(fibers, i, file=filename, )
            fileidx += 1

        output['filename'] = np.array(filenames)
        ascii.write(output, 'fib_coords.dat', overwrite=True)

    fibers.close()
    
tb.file._open_files.close_all() 
   
if __name__ == '__main__':
    main()
