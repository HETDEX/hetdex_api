import sys
import numpy as np
import os.path as op
import time

from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.signal import correlate
from astropy.stats import sigma_clipped_stats

from hetdex_api.config import HDRconfig
from hetdex_api.shot import get_fibers_table

config = HDRconfig()

spec_median = np.load(config.cal_issue)

def get_mask_for_shotid(shotid, savefig=True, ifuslot=None):
    
    fibtab = get_fibers_table(shotid)['fiber_id','multiframe', 'expnum','ifuslot','calfib']
    
    # 5200 Feature
    start_index = 700 # too many peaks in the blue
    end_index = 900
    fiber_specs = np.array( fibtab['calfib'][:, start_index:end_index])
      
    mask = np.zeros_like( np.array( fibtab['calfib']) )
    
    mean, median, stddev = sigma_clipped_stats( np.sum(fiber_specs, axis=1))
    print(mean + 3*stddev)
    for i, fiber_spec in enumerate( fiber_specs):
        if np.sum(fiber_spec) > mean + 3*stddev:
            continue
            
        correlation = correlate( fiber_spec, spec_median, mode='same')

        correlation_peak_locs = np.where( correlation > 3*np.std(correlation))[0]

        if len(correlation_peak_locs) > 0:
            for peak_loc in correlation_peak_locs:
                
                mask[i, start_index + peak_loc-5: start_index + peak_loc+5] = 1
    
    # 5450 feature
    start_index = 900 # too many peaks in the blue
    end_index = 1032
    fiber_specs = np.array( fibtab['calfib'][:, start_index:end_index])
    
    mean, median, stddev = sigma_clipped_stats( np.sum(fiber_specs, axis=1))
    #print(mean + 3*stddev)
    for i, fiber_spec in enumerate( fiber_specs):
        if np.sum(fiber_spec) > mean + 3*stddev:
            continue
            
        correlation = correlate( fiber_spec, spec_median, mode='same')
    
        correlation_peak_locs = np.where( correlation > 3*np.std(correlation))[0]
    
        if len(correlation_peak_locs) > 0:
    
            for peak_loc in correlation_peak_locs:
                
                mask[i, start_index + peak_loc-5: start_index + peak_loc+5] = 1
                
    fibtab['masktmp'] = mask
    
    mask2 = np.ones_like(mask)
    
    mf_to_mask_5200 =[]
    expnum_to_mask_5200 =[]
    mf_to_mask_5450 =[]
    expnum_to_mask_5450 =[]
    mfs = np.unique( fibtab['multiframe'])
    
    for mf in mfs:
    
        for expn in [1,2,3]:
            sel = (fibtab['multiframe']==mf) * (fibtab['expnum']==expn)

            if np.sum( mask[sel][:, 860:880]) > 100:
                mask2[sel, 860:880] = 0
                mf_to_mask_5200.append(mf)
                expnum_to_mask_5200.append(expn)
                
            if np.sum( mask[sel][:, 990:1020]) > 130:
                mask2[sel, 990:1020] = 0
                mf_to_mask_5450.append(mf)
                expnum_to_mask_5450.append(expn)
    
    mask_amps_5200 = Table([ shotid*np.ones(len(mf_to_mask_5200), dtype=int), 
                       mf_to_mask_5200, 
                       expnum_to_mask_5200], 
                      names=['shotid', 'multiframe', 'expnum'])
    mask_amps_5200.write('amps_to_mask/flag5200_{}.txt'.format(shotid), 
                    format='ascii.no_header', overwrite=True)
    mask_amps_5450 = Table([ shotid*np.ones(len(mf_to_mask_5450), dtype=int), 
                        mf_to_mask_5450, 
                       expnum_to_mask_5450], 
                      names=['shotid', 'multiframe', 'expnum'])
    
    mask_amps_5450.write('amps_to_mask/flag5450_{}.txt'.format(shotid), 
                    format='ascii.no_header', overwrite=True)
    
    if savefig:
        fibtab['mask'] = mask2

        if ifuslot is not None:
            plt.figure(figsize=(8,8))
            
            sel_ifu = fibtab['ifuslot'] == ifuslot
            plt.subplot(121)
            plt.title('calfib')
            plt.imshow( fibtab['calfib'][sel_ifu], vmin=-0.1, vmax=0.2, origin='lower')
            
            plt.subplot(122)
            plt.title('cal mask')
            plt.imshow(fibtab['mask'][sel_ifu], origin='lower')
                        
            plt.savefig('figures/mask_{}_ifuslot={}.png'.format(shotid, ifuslot), dpi=150, bbox_inches='tight')

        plt.figure(figsize=(5,100))
        plt.subplot(121)
        plt.title('calfib')
        plt.imshow( fibtab['calfib'], vmin=-0.1, vmax=0.2, origin='lower')
        plt.subplot(122)
        plt.title('cal mask')
        plt.imshow(mask2, origin='lower')
        plt.savefig('figures/mask_{}.png'.format(shotid), dpi=150, bbox_inches='tight')
    return

shotid = int(sys.argv[1])
badshots = np.loadtxt(config.badshot, dtype=int)

if shotid in badshots:
    print('Shotid={} in badshot list. Exiting'.format(shotid))
    sys.exit()

print('Working on {}'.format(shotid))

if op.exists('figures/mask_{}.png'.format(shotid)):
    print('File already made for ', shotid)
    sys.exit()

t0 = time.time()

get_mask_for_shotid(shotid, savefig=True)

t1 = time.time()

print('Done {} in {:3.2f} m'.format(shotid, (t1-t0)/60))
