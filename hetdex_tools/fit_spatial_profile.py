# -*- coding: utf-8 -*-

"""
Created on March 29 2019
@author: Erin Mentuch Cooper

Updated: June 4 2024 EMC - add coordinate/shot/wave/wave_range option

see file /work/05350/ecooper/stampedede2/cosmos/run_profile.py for how to execute

Fits 1D PSF model to line detections using the ML data product h5 file from HDR4.0.0

python3 fit_spatial_profile --start --end

"""

import tables as tb
import numpy as np
import os.path as op
import argparse as ap
import time

from astropy.table import Table, join, unique
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

#import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import ZScaleInterval

from hetdex_api.input_utils import setup_logging
from hetdex_api.config import HDRconfig
from hetdex_api.detections import Detections
from hetdex_api.survey import Survey
from hetdex_tools.interpolate import make_narrowband_image
from hetdex_tools.get_spec import get_spectra
from hetdex_tools import mcmc_gauss

from astropy import wcs
from astropy.visualization import ZScaleInterval
from astropy.stats import sigma_clipped_stats

from astropy import wcs
from astropy.visualization import ZScaleInterval
from astropy.stats import sigma_clipped_stats

from photutils.aperture import aperture_photometry, CircularAnnulus

import pyimfit

import logging
from hetdex_api.input_utils import setup_logging

#plotting preferences
plt.style.use("default")
plt.rcParams["axes.linewidth"] = 2
plt.rcParams.update({"font.size": 12})

plt.rcParams["lines.linewidth"] = 2

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.labelsize"] = 12.0
plt.rcParams["ytick.labelsize"] = 12.0

global config, imsize, pixscale, D_hdr3, D_hdr4

D_hdr3 = Detections('hdr3')
D_hdr4 = Detections('hdr4')
    
try:
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
except:
    print('Could not find data release')
imsize=20
pixscale=0.25

config = HDRconfig('hdr4')

def gaussian(x,u1,s1,A1=1.0,y=0.0):
    if (x is None) or (u1 is None) or (s1 is None):
        return None
    return A1 * (np.exp(-np.power((x - u1) / s1, 2.) / 2.) / np.sqrt(2 * np.pi * s1 ** 2)) + y

# This function converts the row data from the H5 file into an HDU FITS standard
def get_hdu(detectid_obj):
    im_data = fileh.root.LineImages.read_where('detectid == detectid_obj')[0]
    header = fits.Header.fromstring(im_data['header'])
    hdu = fits.PrimaryHDU(im_data['line_image'], header=header)
    hdu_error = fits.ImageHDU(im_data['line_image_err'], header=header)
    hdu_x = fits.ImageHDU(im_data['delta_ra'], header=header)
    hdu_y = fits.ImageHDU(im_data['delta_dec'], header=header)
    return fits.HDUList([hdu, hdu_error, hdu_x, hdu_y])


def fit_profile(detectid=None,
                coords=None,
                wave=None,
                linewidth=None,
                shotid=None,
                name=None,
                extract_class=None,
                apply_mask=False,
                plot=False,
                survey=LATEST_HDR_NAME):

    global config, fileh, D_hdr4, D_hdr3, imsize, pixscale

    if detectid is not None:
        if str( detectid)[0] == '3':
            D = D_hdr3
        if str( detectid)[0] == '4':
            D = D_hdr4
        if str( detectid)[0] == '2':
            D = Detections('hdr2.1')
        
        det_info = D.get_detection_info(detectid)[0]

        wave_obj = det_info['wave']
        flux = det_info['flux']
        linewidth = det_info['linewidth']
        cont = det_info['continuum']
        sn_line = det_info['sn']
        chi2_line = det_info['chi2']

        fwhm = D.get_survey_info(detectid)['fwhm_virus'][0]

        detectid_obj = detectid
        
        #if len(fileh.root.LineImages.read_where('detectid == detectid_obj')) == 0:
        #    return None

        #sn_im = fileh.root.LineImages.read_where('detectid == detectid_obj')[0]['sn_im']

        name = detectid
        hdu = get_hdu(detectid)

    else:
        wave_obj = wave
        
        hdu = make_narrowband_image(
            coords=coords,
            wave_range=[wave - 2*linewidth, wave + 2*linewidth],
            shotid=shotid,
            survey=survey,
            imsize= imsize * u.arcsec,
            include_error=True,
            ffsky=False,
            extract_class=extract_class,
            subcont=False,
            convolve_image=False,
            interp_kind='cubic',
            apply_mask=apply_mask,
            fill_value=0.0,
        )
        
        if name is None:
            name = "{:4.3f}_{:4.2f}_{:4.0f}_{} ".format( coords.ra.deg, coords.ra.deg, wave, shotid)
        S = Survey(survey)
        fwhm = S.fwhm_virus[ S.shotid == shotid][0]
        
    im = hdu[0].data
    error = hdu[1].data
    w = wcs.WCS(hdu[0].header)

    image_data = hdu[0].data

    if np.sum(image_data) == 0:
        return None

    mask = (image_data == 0) 
    mean, median, stddev = sigma_clipped_stats(hdu[0].data[~mask], sigma=2, maxiters=5)
    
    im_x = hdu[2].data
    im_y = hdu[3].data
    im_r = np.sqrt( im_x**2 + im_y**2)

    mask = (image_data == 0) | ( (im_r > 5) & (image_data >= 3*stddev) ) | (im_r>8)

    im = np.ma.masked_where( mask, image_data)

    xpix, ypix = np.shape(image_data)
    xcen = int(xpix / 2)
    ycen = int(ypix / 2)
    dx = int(1.5/pixscale)
    dy = int(1.5/pixscale)

    #create fiber psf                                                                                               
    boxsize=imsize/2
    fiber_pixscale=pixscale

    xl, xh = (0.0 - boxsize / 2.0, 0.0 + boxsize / 2.0 + fiber_pixscale)
    yl, yh = (0.0 - boxsize / 2.0, 0.0 + boxsize / 2.0 + fiber_pixscale)
    x, y = (np.arange(xl, xh, fiber_pixscale), np.arange(yl, yh, fiber_pixscale))
    xgrid, ygrid = np.meshgrid(x, y)

    fiber_psf = np.zeros_like(xgrid)
    r = np.sqrt( xgrid**2 + ygrid**2)
    sel_r = r <= 0.75
    fiber_psf[sel_r] = 1

    # create moffat model for both point source model and core+exp model

    moffat_model_desc = pyimfit.SimpleModelDescription()

    moffat_model_desc.x0.setValue(xcen, [xcen-dx,xcen+dx])
    moffat_model_desc.y0.setValue(ycen, [ycen-dy,ycen+dx])

    moffatmodel = pyimfit.make_imfit_function("Moffat")

    moffatmodel.PA.setValue(0, fixed=True)
    moffatmodel.ell.setValue(0, fixed=True)
    moffatmodel.fwhm.setValue( fwhm/pixscale, fixed=True)
    #moffatmodel.fwhm.setValue( fwhm/pixscale, [1.0, 4.0])
    moffatmodel.beta.setValue( 3.5, fixed=True)
    moffatmodel.I_0.setValue(10, [0.001, 10000])
    moffat_model_desc.addFunction(moffatmodel)


    imfit_fitter_moffat = pyimfit.Imfit(moffat_model_desc, psf=fiber_psf)
    try:
        fit = imfit_fitter_moffat.fit( image_data, mask=mask, error=error)
    except:
        return None
    chi2_moffat = np.float32( imfit_fitter_moffat.reducedFitStatistic)

    # get moffat intensity level
    moffat_fitparams = imfit_fitter_moffat.getRawParameters()
    moffat_x0 = np.float32( moffat_fitparams[0])
    moffat_y0 = np.float32( moffat_fitparams[1])
    moffat_I0 = np.float32( moffat_fitparams[4])

    moffat_im = imfit_fitter_moffat.getModelImage()

    # get RA/DEC of peak emission
    ra_fit, dec_fit = w.wcs_pix2world([[moffat_x0, moffat_y0]], 0)[0]

    dr=0.25
    #dr = 1
    sb = []
    sb_error = []
    sb_moffat_im = []

    r_in_array = np.arange(1, xcen, dr)
    r_in_array = np.logspace(-3, np.log10(2*xcen), num=50)

    for r_in in r_in_array:
        r_out = r_in + dr
        aper = CircularAnnulus((moffat_x0, moffat_y0), r_in, r_out)
        phot_table = aperture_photometry(
            image_data, aper, mask=mask, error=error,
        )
        sb.append(phot_table["aperture_sum"][0] / aper.area)
        sb_error.append(phot_table["aperture_sum_err"][0] / aper.area)

        phot_table = aperture_photometry(moffat_im, aper)
        sb_moffat_im.append(phot_table["aperture_sum"][0] / aper.area)

    sn_max = np.float32( np.max(sb) / (stddev))

    sb = 10 ** -17 * np.array(sb) / (pixscale ** 2)
    sb_error = 10 ** -17 * np.array(sb_error) / (pixscale ** 2)
    sb_moffat_im = 10 ** -17 * np.array(sb_moffat_im) / (pixscale ** 2)

    if plot==True:
        
        fig = plt.figure(figsize=(16,4), constrained_layout=True)
        sb_plus = sb + sb_error
        sb_plusMult = sb_plus/sb
        sb_minus = sb/sb_plusMult

        subplots = fig.subfigures(1,4)

        ax1 = subplots[0].subplots(1,1) #plt.subplot(141,projection=w)

        #ax1.remove()
        #ax1 = fig.add_subplot(141,projection=w)

        #ax1 = fig.add_subplot(1, 4, 1, projection=w)
        ax1.imshow(image_data, vmin= -1 * stddev, vmax=4 * stddev, alpha=0.3, origin='lower')
        ax1.imshow(im, vmin= -1 * stddev, vmax=4 * stddev, origin='lower')

        #plt.xlabel("RA")
        #plt.ylabel("Dec")

        # this helps to bring labels closer to plot borders
        #lon = ax1.coords[0]
        #lat = ax1.coords[1]
        #lon.set_axislabel('RA', minpad=0.5)
        #lat.set_axislabel('Dec', minpad=-0.4)
        #lon.set_ticklabel(exclude_overlapping=True)
        #plt.colorbar()

        if detectid is not None:
            plt.text(
                0.05,
                0.05,
                "S/N im={:3.2f}".format(sn_im),
                size=18,
                color="black",
                transform=ax1.transAxes
            )
        #plt.title(str(detectid))

        plt.text(
            0.05,
            0.9,
            str(name),
            size=16,
            color="black",
            transform=ax1.transAxes
        )

        plt.scatter( moffat_x0, moffat_y0, marker='x',  color='tab:orange')

        ax2 = subplots[1].subplots(1,1)# plt.subplot(142)

        if detectid is not None:
            detectid_obj = detectid
            spec_table = Table(fileh.root.Spec1D.read_where('detectid == detectid_obj'))
            spec_table.rename_column('spec1D','spec')
            spec_table.rename_column('spec1D_err', 'spec_err')
        else:
            # fit at peak position
            
            spec_table = get_spectra( coords=SkyCoord(ra=ra_fit, dec=dec_fit, unit='deg'),
                                      survey=survey,
                                      shotid=shotid,
                                      multiprocess=False,
                                      #apply_mask=True,
                                      loglevel='WARNING')
            
        wave_rect = 2.0 * np.arange(1036) + 3470.0
        #plt.plot(wave_rect, spec_table['spec1D'][0]*10**-17 * u.erg / (u.cm ** 2 * u.s * u.AA))
        plt.errorbar(wave_rect, 
                     spec_table['spec'][0], 
                     yerr=spec_table['spec_err'][0], elinewidth=1, ecolor='lightgrey')
        plt.xlabel('wavelength (AA)')
        plt.ylabel('spec ergs/s/cm^2/AA')

        if detectid is not None:
            model_gauss= gaussian(wave_rect,u1=wave_obj,s1=linewidth,A1=flux,y=cont)
            plt.plot(wave_rect, model_gauss, label='sn={:3.2f} lw={:3.2f} chi2={:3.2f}'.format(sn_line, linewidth, chi2_line))

            plt.text(
                0.05,
                0.05,
                "S/N={:3.2f}".format(sn_line),
                size=18,
                transform=ax2.transAxes,
                color="black",
            )
            
        else:
            #instatiate an MCMC_Gauss object and populate
            #you can optionally pass a logger instance to the constructor, otherwise it will make its own 
            logger = setup_logging()
            logger.setLevel(logging.WARNING)
            fit = mcmc_gauss.MCMC_Gauss(logger=logger)
            #set the initial guesss
            #(here you can see it is set wrong to show we converge on the correct answer)
            fit.initial_A = 25
            fit.initial_y = 0
            fit.initial_sigma = 5
            fit.initial_mu =wave_obj
            fit.initial_peak = None

            #set the data to fit
            fit.data_x = wave_rect
            fit.data_y = spec_table['spec'][0]
            fit.err_y = spec_table['spec_err'][0]
            fit.err_x = np.zeros(len(fit.err_y))

            
            #these are the defaults and don't have to be set
            fit.max_sigma = 10.0
            fit.min_sigma = 1.7
            fit.range_mu = 5.0
            fit.max_A_mult = 2.0
            fit.max_y_mult = 2.0
            fit.min_y = -10.0 

            fit.burn_in = 250
            fit.main_run = 1000
            fit.walkers = 100

            fit.run_mcmc()
            wave_fit = fit.mcmc_mu[0]
            linewidth_fit = fit.mcmc_sigma[0]
            flux_fit =fit.mcmc_A[0]
            cont_fit = fit.mcmc_y[0]
            sn_line = fit.mcmc_snr

            model_gauss = gaussian(wave_rect,u1=wave_fit,s1=linewidth_fit, A1=flux_fit,y=cont_fit)

            sel_wave = np.abs( wave_rect - wave_fit) < 50
            y = spec_table['spec'][0][sel_wave]
            fit = gaussian(wave_rect[sel_wave],u1=wave_fit,s1=linewidth_fit, A1=flux_fit,y=cont_fit)
            yerr = spec_table['spec_err'][0][sel_wave]
            n_free = 3
            N = np.sum(sel_wave)
            
            chi2_line = 1.0/(N-n_free) * np.sum(((fit - y)/yerr)**2)
            
            plt.plot(wave_rect, model_gauss, label='sn={:3.2f} lw={:3.2f} chi2={:3.2f}'.format(sn_line, linewidth_fit, chi2_line))
            
        ymax = np.max( spec_table['spec'][0][ np.where( np.abs( wave_rect-wave_obj) < 100)[0]])
        ymin = np.min( spec_table['spec'][0][ np.where( np.abs( wave_rect-wave_obj) < 100)[0]])
        plt.ylim(1.5*ymin, 1.5*ymax)
        
        plt.xlim(wave_obj - 100,wave_obj+100)


        plt.text(
            0.05,
            0.05,
            "S/N={:3.2f}".format(sn_line),
            size=18,
            transform=ax2.transAxes,
            color="black",
        )

        plt.legend()

        ax3 = subplots[2].subplots(1,1) #plt.subplot(143)
        plt.errorbar(r_in_array*pixscale, sb, yerr=sb_error, label='data',
                 linestyle="none", color='tab:blue', marker="o")

        plt.plot(r_in_array*pixscale, sb_moffat_im,
             label='Moffat Model ($\chi2$={:3.2f})'.format(chi2_moffat),
             color='tab:orange')
        plt.xlim(0,10)
        plt.xlabel('Semi-major axis (arcsec)')
        plt.axhline(y=(stddev/pixscale**2)*10**-17, color='grey', linestyle='dotted', label='stddev background')
        #plt.axhline(y=sb_bkg*10**-17, color='cyan', linestyle='dashed', label='stdev background')
        #plt.axhline(y=(mean)*10**-17/pixscale**2, color='grey', linestyle='dashed', label='subtracted mean background')

        plt.text(
        0.05,
        0.05,
        "S/N max={:3.2f}".format(sn_max),
        size=18,
        transform=ax3.transAxes,
        color="black",
        )
        plt.ylabel(r'SB (erg/s/cm$^2$/arcsec$^2$)')
        plt.legend(fontsize=10)
        plt.ylim(10**-21,10**-15)
        plt.yscale('log')


        if detectid is not None:
            detectid = detectid_obj
            # plot each fiber for 4th object in example table
            obj_data = fileh.root.FiberImages.read_where('detectid == detectid_obj')[0]
            height=9 # in pixels
            wave = obj_data['im_wave']
            im_sum = obj_data['im_sum'] # this is the 2D summed image, 1st dim is height in fiber dims, 2nd dim is wave dim
            im_array = obj_data['im_array'] # this is the 4 brightest fibers, 1st dim is fibers, 2nd dim is fiber dims, 3rd is wavelength
            zscale = ZScaleInterval(contrast=0.5,krej=2.5)
            vmin, vmax = zscale.get_limits(values=im_sum)

            ax4 = subplots[3].subplots(4,1, sharex=True) #plt.subplot(144)

            if wave[25] < 3470:
                start = np.where( wave>= 3470)[0][0]
                end = 75
            elif wave[75] > 5540:
                start=25
                end = np.where( wave >= 5540)[0][0]
            else:
                start = 25
                end = 75

            ax4[0].imshow(im_sum[:, start:end],
                          vmin=vmin, vmax=vmax,
                          extent=[wave[start], wave[end], -int(height/2.), int(height/2.)], 
                          origin="lower",cmap=plt.get_cmap('gray'),interpolation="none", aspect=2)
            axis_i = 1
            
            for im_i in np.arange(0,3):
                zscale = ZScaleInterval(contrast=0.5,krej=2.5)
                vmin, vmax = zscale.get_limits(values=im_array[im_i, :, start:end])
                ax4[axis_i].imshow(im_array[im_i, :, start:end],vmin=vmin, vmax=vmax,extent=[wave[start], wave[end], -int(height/2.), int(height/2.)], origin="lower",
                                   cmap=plt.get_cmap('gray'),interpolation="none", aspect=2)

                axis_i += 1
                
        plt.savefig('figures/fits_{}.png'.format(name))

    if detectid is not None:
        return detectid, sn_im, sn_max, moffat_x0, moffat_y0, chi2_moffat, sn_line, chi2_line, linewidth
    else:
        return name, sn_max, moffat_x0, moffat_y0, chi2_moffat, sn_line, chi2_line, linewidth_fit



def get_parser():
    """ Function that returns a parser"""

    parser = ap.ArgumentParser(
        description="""Create fiber cutout for a given position or detectid""",
        add_help=True,
    )

    parser.add_argument(
        "-s", "--start", help="""Index to start at""", type=int, default=0,
    )

    parser.add_argument(
        "--end",
        "-e",
        type=int,
        help="""Index to end at""",
        default=10,
    )

    return parser


def main(argv=None):
    """ Main Function """

    parser = get_parser()
    args = parser.parse_args(argv)
    args.log = setup_logging()

    args.log.info(args)

    filename = 'fit_profile/output_{}_{}.fits'.format(args.start, args.end)

    global config, D_hdr3, D_hdr4, fileh

    D_hdr3 = Detections('hdr3')
    D_hdr4 = Detections('hdr4')

    mlfile = op.join( config.hdr_dir['hdr4'], 'catalogs','ml','detect_ml_4.0.0.h5')
    fileh = tb.open_file(mlfile, 'r') 

    detlist = fileh.root.LineImages.read_coordinates( np.arange( args.start,  args.end), 'detectid')

    t0 = time.time()
    res = []
    for det in detlist:
        try:
            res.append( fit_profile(det))
        except:
            args.log.info('Failed for {}'.format(det))
            
    det_list = []
    sn_im_list = []
    sn_max_list = []
    moffat_x0_list = []
    moffat_y0_list = []
    chi2_moffat_list = []
    sn_line_list = []
    chi2_line_list = []
    linewidth_list = []
    for r in res:
        if r is None:
            continue
        
        detectid, sn_im, sn_max, moffat_x0, moffat_y0, chi2_moffat, sn_line, chi2_line, linewidth = r
        det_list.append(detectid)
        sn_im_list.append(sn_im)
        sn_max_list.append(sn_max)
        moffat_x0_list.append(moffat_x0)
        moffat_y0_list.append(moffat_y0)
        chi2_moffat_list.append(chi2_moffat)
        sn_line_list.append(sn_line)
        chi2_line_list.append(chi2_line)
        linewidth_list.append(linewidth)

    dets = np.array(det_list)
    sn_im = np.array(sn_im_list)
    sn_max = np.array(sn_max_list)
    moffat_x0 = np.array(moffat_x0_list)
    moffat_y0 = np.array(moffat_y0_list)
    chi2_moffat = np.array(chi2_moffat_list)
    sn_line = np.array(sn_line_list)
    chi2_line = np.array(chi2_line_list)
    linewidth = np.array(linewidth_list)

    output = Table( [dets, sn_im, sn_max, moffat_x0, moffat_y0, chi2_moffat, sn_line, chi2_line, linewidth], 
                names=['detectid', 'sn_im', 'sn_max', 'moffat_x0', 'moffat_y0', 'chi2_moffat', 'sn_line', 'chi2_line', 'linewidth'])

    output.write(filename, overwrite=True)

    t1 = time.time()

    fileh.close()
    D_hdr4.close()
    D_hdr3.close()
    
    args.log.info('Completed indices {} to {} in {:3.2f} min'.format(args.start, args.end, (t1-t0)/60))

if __name__ == "__main__":
    main()
