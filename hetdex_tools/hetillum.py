# -*- coding: utf-8 -*-
"""

Function to calculate het pupil illumination depending
on tracker position and field location

author = Niv Drory, Hanshin Lee
"""

import matplotlib.pyplot as plt
from astropy.io import fits
import subprocess

illumlib = '/work/03946/hetdex/hdr1/software/illum_lib'

def hetillum_(tx,ty,tr,fx,fy,platform=False,plot=False):
    """
    Calculate and return HET pupil illumination relative to 
    center track, center field.

    Inputs: tx,ty,tr tracker x and y in mm, rho in degrees
    fx, fy: field position in degrees
    platform=: Include obscuration by work platform. The work platform was
    place prior to Nov 6, 2017.
    plot=: plot image of pupil with obstructions

    Returns: illumnination relative to 0,0,0 tracker and 0,0 field.
    Note that the value can be very slightly greater than 1 due to numerics.
    """
    if platform==True:
        f = '-p'
    else:
        f = ''
    if plot==True:
        p = '-i'
    else:
        p = ''
    s = illumlib+'/hetillum -x %s %s [%f,%f,%f] [%f,%f] 256'% (f,p,tx,ty,tr,fx,fy)
    # print(s)
    proc = subprocess.Popen(s.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o, e = proc.communicate()

    if plot==True:
        hdulist = fits.open('pupilmask.fits')
        my_image = hdulist[0].data[0,0:,0:]
        plt.imshow(my_image, cmap='gray')
        plt.show()

    v = float(o.decode('ascii'))
    return v
