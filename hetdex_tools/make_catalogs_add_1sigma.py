#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os.path as op
import time

from astropy.coordinates import SkyCoord
from astropy.table import Table, join

import extinction

from hetdex_api.extract import Extract
from hetdex_api.config import HDRconfig
from hetdex_api.flux_limits.shot_sensitivity import ShotSensitivity

from multiprocessing import Pool


version = '3.0.0'

config = HDRconfig()
# Enter the catalog version
version = '3.0.0'

config = HDRconfig()
catfile = op.join(config.hdr_dir['hdr3'], 'catalogs', 'source_catalog_' + version + '.fits')
source_table = Table.read(catfile)

print('Source catalog was found at {}'.format(catfile))


source_table = Table.read(catfile)

print('Source catalog was found at {}'.format(catfile))

shotlist = np.unique(source_table['shotid'])


def get_flim_sig_erin(detectid=None, coord=None, wave=None, datevobs=None, shotid=None):
    """
    Script to grab flim_1sigma and sn_sig on demand from calibrated fiber extractions
    
    Parameters
    ----------
    detectid int
       detectid for source 
    coord
        astropy SkyCoord object
    wave
        central wavelength
    shotid 
        observation ID
    
    Returns
    -------
    flim_1sigma
        the 1 sigma sensitivity calculated over 7 pixels of the PSF-weighted extracted spectra
    """
    
    detectid_obj = detectid
    det_info = source_table[source_table['detectid'] == detectid][0]
    shotid = det_info['shotid']
    wave = det_info['wave']
    coord = SkyCoord(ra=det_info['ra'], dec=det_info['dec'], unit='deg')        
    fwhm = det_info['fwhm']
    
    if datevobs is None:
        datevobs = '{}v{}'.format( str(shotid)[0:8], str(shotid)[8:])
        
    if shotid is None:
        shotid_obj = int(datevobs[0:8] + datevobs[9:])
    else:
        shotid_obj = shotid
    
    try:
        E = Extract()
        E.load_shot(datevobs, fibers=False)
    
        info_result = E.get_fiberinfo_for_coord(coord, radius=3.5, fiber_lower_limit=2 )

        ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result

        #print(len(ifux))

        moffat = E.moffat_psf(fwhm, 10.5, 0.25)
   
        I = None
        fac = None

        weights, I, fac = E.build_weights(xc, yc, ifux, ifuy, moffat, 
                                            I=I, fac=fac, return_I_fac = True)

        norm = np.sum(weights, axis=0)

        weights = weights/norm   

        result = E.get_spectrum(data, error, mask, weights, 
                                remove_low_weights=False)


        spec, spec_err = [res for res in result]

        w_index = np.where(E.wave >= wave)[0][0]

        nfib = np.shape(weights)[0]

        flim_1sigma = 2 * np.sqrt( np.nansum( (spec_err[w_index-3:w_index+4])**2))

        sn_sig = 2 * np.sum( spec[w_index-3:w_index+4]) / flim_1sigma

        npix = np.sum(np.isfinite(spec[w_index-3:w_index+4]))

        apcor = np.sum( norm[w_index-3:w_index+4])/len( norm[w_index-3:w_index+4])  

        E.close()

        # XXX divide by apcor .... 
        return flim_1sigma/apcor, apcor
    except:
        return 999


def get_noise_1sigma(shot):
    sncut = 1

    datevobs = str(shot)[0:8] + 'v' + str(shot)[8:12]
    
    sel_shot = (source_table['shotid'] == shot) * (source_table['det_type'] == 'line')

    detectid = source_table['detectid']
    
    if True:
        shot_sens_v1 = ShotSensitivity(datevobs,
                                       flim_model="v4",
                                       log_level="INFO")

        ra = source_table['ra'][sel_shot]
        dec = source_table['dec'][sel_shot]
        wave = source_table['wave'][sel_shot]
        
        f_1sigma, apcor = shot_sens_v1.get_f50(ra,
                                               dec,
                                               wave,
                                               sncut,
                                               direct_sigmas=True)
        shot_sens_v1.close()
        np.savetxt('f1sigma/{}.txt'.format(shot), np.array(detectid[sel_shot]), f_1sigma*10**17, apcor)
        return np.array(detectid[sel_shot]), f_1sigma*10**17, apcor
    else:#except:
        print('Failed to get sensitivity value for {}'.format(datevobs))


t0 = time.time()
p = Pool()
res = p.map(get_noise_1sigma, shotlist)
p.close()
t1 = time.time()

print('Sensitivity values complete in {}m'.format((t1-t0)/60))

np.shape(res)

dets = []
noise_1sigma = []
apcor = []
for r in res:
    if r is None:
        continue
    else:
        dets.extend(r[0])
        noise_1sigma.extend(r[1])
        apcor.extend(r[2])

new = Table([dets, noise_1sigma, apcor], names=['detectid',
                                                'flux_noise_1sigma_obs',
                                                'apcor_api'])
new.write('3.0.0_1sigma_api.fits')

bad = np.where(new['flux_noise_1sigma_obs'] > 990)[0]

for bad_i in bad:
    f1sigma, apcor = get_flim_sig_erin(new['detectid'][bad_i])
    if f1sigma == 0.0:
        f1sigma = 999
        print('Bad: {}'.format(new['detectid'][bad_i]))
    new['flux_noise_1sigma_obs'][bad_i] = f1sigma
    new['apcor_api'][bad_i] = apcor

join_new = join(source_table, new, join_type='left')

print(len(join_new), len(source_table))

# get intrinsic flux limit by dereddening

import extinction

Rv = 3.1
ext = []

for index in np.arange( np.size(join_new['detectid'])):
    src_wave = np.array([np.double(join_new['wave'][index])])
    ext_i = extinction.fitzpatrick99(src_wave, join_new['Av'][index], Rv)[0]
    ext.append(ext_i)

deredden = 10**(0.4*np.array(ext))

join_new['flux_noise_1sigma'] = deredden * join_new['flux_noise_1sigma_obs']
#join_new['flux_noise_1sigma'] = join_new['flux_noise_1sigma'].filled(999)
#join_new['flux_noise_1sigma_obs'] = join_new['flux_noise_1sigma_obs'].filled(999)
join_new['apcor_api'] = join_new['apcor_api']#.filled(0.0)

print(len(source_table), len(join_new))

join_new.write('source_catalog_3.0.0_1sigma.fits', overwrite=True)

print('Number of bad flim values is {}'.format(len(bad)))
