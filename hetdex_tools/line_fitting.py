# -*- coding: utf-8 -*-
"""

Some line fitting functions using the specutils package
Functions with output spectra from get_spec.py

This is very preliminary and should not be used without
discussing it with Erin.

author = Erin Mentuch Cooper
"""
import os
import os.path as op
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, join, Column
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty

from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import estimate_line_parameters, fit_lines
from specutils.fitting import fit_generic_continuum
from specutils.manipulation import extract_region
from specutils.analysis import line_flux, snr


def calc_chi2(sub_spectrum, g_fit, reduced=True, n_free=2):
    '''
    Calculate the chi2 of a Spectrum1D region given
    a model spectrum.

    Parameters
    ----------
    sub_spectrum
        Spectral Region of Spectrum1D object
    g_fit
        Fitted model generated from fit_lines function from specutils.fitting
    reduced
        Boolean flag to calculate a reduced chi2 value. Default is True
    n_free
        number of parameters we are fitting, default is 2 assuming gaussian
        model
    '''

    fit = np.array(g_fit(sub_spectrum.spectral_axis))
    y = np.array(sub_spectrum.flux)
    yerr = np.array(sub_spectrum.uncertainty.array)
    N = np.size(fit)

    if reduced:
        chi2 = 1.0/(N-n_free) * np.sum(((fit - y)/yerr)**2)
    else:
        chi2 = np.sum(((fit - y)/yerr)**2)

    return chi2


def line_fit(spec, spec_err, wave_obj, dwave=10.*u.AA,
             dwave_cont=100.*u.AA, sigmamax=14.*u.AA):
    '''
    Function to fit a 1D gaussian to a HETDEX spectrum from get_spec.py

    Parameters
    ----------
    spec
        1D spectrum from a row in the table provided by get_spec.py.
        Will assume unit of 10**-17*u.Unit('erg cm-2 s-1 AA-1') if no units
        are provided.
    spec_err
        1D spectral uncertainty from table provided by get_spec.py.
        Will assume unit of 10**-17*u.Unit('erg cm-2 s-1 AA-1') if no units
        are provided.
    wave_obj
        wavelength you want to fit, an astropy quantity
    dwave
        spectral region above and below wave_obj to fit a line, an astropy quantity.
        Default is 10.*u.AA
    dwave_cont
        spectral region to fit continuum. Default is +/- 100.*u.AA
    sigmamax
        Maximum linewidth (this is sigma/stdev of the gaussian fit) to allow
        for a fit. Assumes unit of u.AA if not given

    Returns
    -------

    '''

    try:
        spectrum = Spectrum1D(flux=spec,
                              spectral_axis=(2.0*np.arange(1036)+3470.)*u.AA,
                              uncertainty=StdDevUncertainty(spec_err),
                              velocity_convention=None)
    except ValueError:
        spectrum = Spectrum1D(flux=spec * 10**-17 * u.Unit('erg cm-2 s-1 AA-1'),
                              spectral_axis=(2.0 * np.arange(1036) + 3470.) * u.AA,
                              uncertainty=StdDevUncertainty(spec_err * 10**-17 * u.Unit('erg cm-2 s-1 AA-1')),
                              velocity_convention=None)

    # measure continuum over 2*dwave_cont wide window first:
    cont_region = SpectralRegion((wave_obj-dwave_cont), (wave_obj+dwave_cont))
    cont_spectrum = extract_region(spectrum, cont_region)
    cont = np.median(cont_spectrum.flux)

    if np.isnan(cont):
        #set continuum if its NaN
        print('Continuum fit is NaN. Setting to 0.0')
        cont = 0.0*cont_spectrum.unit
        
    # now get region to fit the continuum subtracted line
    
    sub_region = SpectralRegion((wave_obj-dwave), (wave_obj+dwave))
    sub_spectrum = extract_region(spectrum, sub_region)

    try:
        line_param = estimate_line_parameters(sub_spectrum-cont, models.Gaussian1D())
    except:
        return None

    if np.isnan(line_param.amplitude.value):
        print('Line fit yields NaN result. Exiting.')
        return None
        
    try:
        sigma = np.minimum(line_param.stddev, sigmamax)
    except ValueError:
        sigma = np.minimum(line_param.stddev, sigmamax*u.AA)

    if np.isnan(sigma):
        sigma=sigmamax
        
    g_init = models.Gaussian1D(amplitude=line_param.amplitude,
                               mean=line_param.mean, stddev=sigma)
    
#    lineregion = SpectralRegion((wave_obj-2*sigma), (wave_obj+2*sigma))
#    cont = fit_generic_continuum(sub_spectrum, exclude_regions=lineregion,
#                                 model=models.Linear1D(slope=0))

    #r1 = SpectralRegion((wave_obj-dwave), (wave_obj-2*sigma))
    #r2 = SpectralRegion((wave_obj+2*sigma), (wave_obj+dwave))
    #fitcontregion = r1 + r2

    #fit_cont_spectrum = extract_region(sub_spectrum, fitcontregion)
    #cont = np.mean(np.hstack([fit_cont_spectrum[0].flux, fit_cont_spectrum[1].flux]))
    
    #contspec = cont(sub_spectrum.spectral_axis)

    g_fit = fit_lines(sub_spectrum - cont, g_init)

    x = np.arange(wave_obj.value-dwave.value,
                  wave_obj.value+dwave.value, 0.5)*u.AA
    y_fit = g_fit(x)

    line_flux_model = np.sum(y_fit*0.5*u.AA)

    chi2 = calc_chi2(sub_spectrum - cont, g_fit)

    sn = np.sum(np.array(sub_spectrum.flux )) / np.sqrt(np.sum(
        sub_spectrum.uncertainty.array**2))
    
    line_flux_data = line_flux(sub_spectrum-cont).to(u.erg * u.cm**-2 * u.s**-1)
    
    line_flux_data_err = np.sqrt(np.sum(sub_spectrum.uncertainty.array**2))

    #fitted_region = SpectralRegion((line_param.mean - 2*sigma),
    #                               (line_param.mean + 2*sigma))

    #fitted_spectrum = extract_region(spectrum, fitted_region)

    #line_param = estimate_line_parameters(fitted_spectrum, models.Gaussian1D())

    #sn = np.sum(np.array(fitted_spectrum.flux)) / np.sqrt(np.sum(
    #    fitted_spectrum.uncertainty.array**2))

    #line_flux_data = line_flux(fitted_spectrum).to(u.erg * u.cm**-2 * u.s**-1)

    #line_flux_data_err = np.sqrt(np.sum(fitted_spectrum.uncertainty.array**2))

    return line_param, sn, chi2, sigma, line_flux_data, line_flux_model, line_flux_data_err, g_fit, cont



def make_line_catalog(input_table, sources, shotidmatch=False):
    '''
    Function to perform line fitting and generate a catalog for
    a list of sources. Requires the input and output tables
    used to extract HETDEX spectra using hetdex_tools/get_spec.py

    Line fitting figures are created in a directory called line_fits
    in the working directory.
    
    Parameters
    ---------
    input_table
        astropy table with columns 'ID','ra','dec','wave'
        and optional 'shotid' to use with shotidmatch option
    sources
        output table from get_spectra or hetdex_get_spec
    shotidmatch
        boolean flag to match shotid in input_table to sources

    Output
    ------
    output_table
        astropy table with line fit parameters

    '''
    
    detectid = []
    shotid = []
    chi2_fit = []
    sn_fit = []
    wave_fit = []
    amp_fit = []
    sigma_fit = []
    lfdata_fit = []
    lfmodel_fit = []
    lfdata_err_fit = []
    g_fit_array = []
    cont_fit = []
    
    for row in sources:
        spec = row['spec']
        spec_err = row['spec_err']

        if shotidmatch == True:
            sel_obj = np.where( (input_table['ID'] == row['ID']) *
                                (input_table['shotid'] == row['shotid']))[0]
        else:
            sel_obj = np.where( (input_table['ID'] == row['ID']))[0]

        try:

            wave_obj = input_table['wave'][sel_obj]

            result = line_fit(spec*sources['spec'].unit,
                              spec_err*sources['spec_err'].unit,
                              wave_obj*u.AA)

            detectid.append(row['ID'])
            shotid.append(row['shotid'])

            if result is not None:
                line_param, sn, chi2, sigma, line_flux_data, line_flux_model, line_flux_data_err, g_fit, cont = result
            
                chi2_fit.append(chi2)
                sn_fit.append(sn)
                wave_fit.append(line_param.mean.value)
                amp_fit.append(line_param.amplitude.value)
                sigma_fit.append(line_param.stddev.value)
                lfdata_fit.append(line_flux_data.value)
                lfmodel_fit.append(line_flux_model.value)
                lfdata_err_fit.append(line_flux_data_err)
                cont_fit.append(cont.value)
            else:
                chi2_fit.append(np.nan)
                sn_fit.append(np.nan)
                wave_fit.append(np.nan)
                amp_fit.append(np.nan)
                sigma_fit.append(np.nan)
                lfdata_fit.append(np.nan)
                lfmodel_fit.append(np.nan)
                lfdata_err_fit.append(np.nan)
                cont_fit.append(np.nan)
        except:

            detectid.append(row['ID'])
            shotid.append(row['shotid'])
            chi2_fit.append(np.nan)
            sn_fit.append(np.nan)
            wave_fit.append(np.nan)
            amp_fit.append(np.nan)
            sigma_fit.append(np.nan)
            lfdata_fit.append(np.nan)
            lfmodel_fit.append(np.nan)
            lfdata_err_fit.append(np.nan)
            cont_fit.append(np.nan)

        output=Table()
        output.add_column(Column(detectid), name='ID')
        output.add_column(Column(shotid), name='shotid')
        output.add_column(Column(chi2_fit), name='chi2_fit')
        output.add_column(Column(wave_fit), name='wave_fit')
        output.add_column(Column(sn_fit), name='sn_fit')
        output.add_column(Column(amp_fit), name='amp_fit')
        output.add_column(Column(sigma_fit), name='sigma_fit')
        output.add_column(Column(lfdata_fit), name='line_flux_data')
        output.add_column(Column(lfmodel_fit), name='line_flux_model')
        output.add_column(Column(lfdata_err_fit), name='line_flux_data_err')
        output.add_column(Column(cont_fit), name='cont_fit')
        
        output_tab = join(input_table, output, keys=['ID'])

    return output_tab

def plot_line(objid, sources, wave_obj=None, shotid=None, save=False):
    """
    Function to plot up objid at a wavelength in an extracted source
    table from get_spectra
    
    Parameters
    ----------
    objid
       object name to match in the 'ID' column. Make sure dtype
       maches that in the sources table
    sources
       astropy table with spectra. Produced by get_spectra
    wave_obj
       wavelength you want to fit around. 
    shotid
       shotid to get spectrum from. Will return shotid list
       with spectra if not provided
    save
       boolean flag to save line fit to a png
    Returns
    -------
    a matplotlib figure
    """

    if shotid is None:
        sel_obj = sources['ID'] == objid
        shots = np.array( np.unique(sources['shotid'][sel_obj]))
        print('Source ' + str(objid) + ' is found in shotids: ', shots)
        return None
        
    sel_obj = (sources['ID'] == objid) * (sources['shotid'] == shotid)

    if np.sum(sel_obj) != 1:
        print('No unique match found')
        return None
    
    if True:
        plt.figure(figsize=(10,8))
        spec = np.array(sources['spec'][sel_obj]).flatten()
        spec_err = np.array(sources['spec_err'][sel_obj]).flatten()

        line_param, sn, chi2, sigma, line_flux_data, line_flux_model, line_flux_data_err, g_fit, cont=line_fit(
            spec*sources['spec'].unit,
            spec_err*sources['spec_err'].unit,
            wave_obj*u.AA)

        wave = (2.0 * np.arange(1036) + 3470.)
        sel_w = (wave >= wave_obj - 50) * ( wave <= wave_obj + 50)
        plt.errorbar(wave[sel_w], spec[sel_w], yerr=spec_err[sel_w],fmt='o', label='extract')
        x = np.arange(wave_obj-50, wave_obj+50, 0.5)*u.AA
                
        plt.plot(x, g_fit(x).value + cont.value*np.ones(np.size(x.value)), 'r', label='model')
        #plt.plot(x, cont(x),'b-', label='cont')
        plt.axhline(y = cont.value,label='cont', color='green', linestyle='dashed')
        
        plt.xlabel('Spectral Axis ({})'.format(u.AA))
        plt.ylabel('Flux Axis({})'.format(sources['spec'].unit))
        plt.title('ID = {:8s} SN = {:4.2f}  Chi2 = {:4.2f}  sigma = {:4.2f}'.format(str(objid), sn, chi2, line_param.stddev.value))
        plt.legend()

        if save:
            if not op.exist('line_fits'):
                os.makedirs('line_fits')
            plt.savefig('line_fit_ID' + str(row['ID']) + '_' + str(shotid) + '.png')


        return line_param, sn, chi2, sigma, line_flux_data, line_flux_model, line_flux_data_err, g_fit, cont

    else:
        print('No spectrum found for ' + str(objid) + ' in shotid = ' + str(shotid))
        return None
