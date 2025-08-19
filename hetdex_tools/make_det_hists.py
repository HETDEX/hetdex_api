"""

Make histograms/tables of detected and undetected
sources from the source simulations

Author: D Farrow (Hull; 2025)

"""
import sys
from numpy import (linspace, digitize, 
                   log10, logspace, sum, array, 
                   int8, nan)
from scipy.special import erf
from numba import jit
from astropy.table import Table

from hetdex_api.flux_limits import flim_models

@jit
def make_frac_dets(sns, bwaves, blws, bflux, bsig, cmpls,
                   nwaves, nlws, nflux, nsig):
    """
    Count number of detected/undetected objects 
    in cells of flux, LW, etc.

    Slow in base python, if can get Numba working
    will be quicker

    """
    
    i = 0
    rows = []
    rows_f50 = []
    
    for sncut in [4.8, 5.0, 5.5, 6.0]:
        sel_sn = sns > sncut
        # model completeness as a test
        tcmpl = cmpls[i]
        i = i + 1
        for iw in range(1, nwaves):
            sel_w = (bwaves == iw)
            for ilw in range(1, nlws):
                sel_lw = sel_w & (blws == ilw)
                for isig in range(1, nsig):
                    sel_s = sel_lw & (bsig == isig)
        
                    f50 = 0
                    frac_last = 0
                    
                    for ib in range(1, nflux):

                        # Save total detected/undetected
                        sel = sel_s & (bflux == ib)
                        ndet = sum((sel & sel_sn))
                        ntot = sum(sel)
              
                        # total expected to be detected by model
                        model_ndet = sum(tcmpl[sel])

                        rows.append([sncut, iw, ilw, ib, isig, 
                                     ndet, ntot, model_ndet])

                        # Compute f50 for this bin
                        if ntot > 0 and ib > 1:
                            frac = ndet / ntot
                            if frac > 0.5 and f50 < 1e-19:
                                f50 = fbcens[ib - 1] + (0.5 - frac)*(fbcens[ib - 1] - fbcens[ib - 2]) / (frac - frac_last) 
                                rows_f50.append([sncut, iw, ilw, isig, f50])
                            else:
                                frac_last = frac

                    if f50 < 1e-19:
                        rows_f50.append([sncut, iw, ilw, isig, nan])
                                
    return rows, rows_f50


# Define the bins
wbins = linspace(3500, 5600, 8)
lwbins = linspace(2, 6, 5)
bins = linspace(1e-17, 5e-16, 15)
sigma_bins = logspace(log10(1.8e-17), 
                      log10(5e-17), 8)
sbcens = 0.5*(sigma_bins[1:] + sigma_bins[:-1])
fbcens = 0.5*(bins[1:] + bins[:-1])

# Read in the data
table = Table.read(sys.argv[1])
table["fluxin"] = table["fluxin"]*1e-17

# Compute the f50 value for the model_ver for testing
cmpls = []
for sncut in [4.8, 5.0, 5.5, 6.0]:
    model_ver = "v5_simple"
    model, sinterp, _ = flim_models.return_flux_limit_model(model_ver)
    f50s = model(table["sigma_3_3.5"], table["wavein"], sncut, 
                 linewidth=table["LWin"])
    frac_det = sinterp(table["fluxin"], f50s, 
                       table["wavein"], sncut)

    # Symbolic regression derived fitting function
    cmpls.append(frac_det)
    table[f'detected_{sncut}_{model_ver}'] = frac_det

# Need numpy array for Numba
cmpls = array(cmpls) 
print("Added model completeness...")


# Identify which of the bins each source is in
names = ["wavein", "LWin", "fluxin", "sigma"]
for n, tbins in zip(names, [wbins, lwbins, bins, sigma_bins]):
    
    if n == "sigma":
        table[n + "_bin"] = int8(digitize(table["sigma_3_3.5"], 
                                     bins=sigma_bins))
    else:
        table[n + "_bin"] = int8(digitize(table[n], bins=tbins))

if len(sys.argv) == 5:
    table.write(sys.argv[4])
   
# Make the 2D histograms import cast to arrays so Numba works
rows, rows_f50 = make_frac_dets(array(table["SN"], dtype="float64"), 
                                array(table["wavein_bin"], dtype="int32"),
                                array(table["LWin_bin"], dtype="int32"), 
                                array(table["fluxin_bin"], dtype="int32"), 
                                array(table["sigma_bin"], dtype="int32"),  
                                cmpls, len(wbins), 
                                len(lwbins), len(bins), len(sigma_bins))

# Save to a file
table_out = Table(array(rows), 
              names = ["SNcut", "wave_bin", "LW_bin", 
                       "flux_bin", "sigma_bin", "Ndet", 
                       "Ntot", f"model_ndet_{model_ver}"],
              dtype = ["float64", "int32", "int32",
                       "int32", "int32", "int32", "int32", 
                       "float64"]
             )
table_out.write(sys.argv[2])

table_out_f50 = Table(array(rows_f50), 
              names = ["SNcut", "wave_bin", "LW_bin", 
                       "sigma_bin", "f50"],
              dtype = ["float64", "int32", "int32",
                       "int32", "float64"]
             )
table_out_f50.write(sys.argv[3])
