# Sensitivity Cubes

## Overview 

The flux limits in HETDEX are stored in ra, dec, wavelength datacubes ("sensitivity cubes"). The 
value stored is the inverse of the PSF weighted 1-sigma noise, which can be used to
compute flux limits. This is how we can store the flux limit as a function of position and wavelength.


## The old, deprecated flux limit models

In the old, deprecated models (hdr1, hdr2pt1), to characterise the fraction of sources detected as a function 
of flux, a parameterisation from [Fleming et al 1995](http://adsabs.harvard.edu/abs/1995AJ....109.1044F) is used. 
One input to this parameterisation are the flux limits stored in the datacubes. The second and final parameter 
in this function controls how quickly the detection rate falls off with flux. This was set to a fixed 
value, calibrated from simulations. This wasn't used in the end as we developed a more sophisticated method.

## The new flux limit models

The new models store the wavelength dependent shape of the completeness versus flux for different S/N cuts 
in a series of files. The existing models and a description of them is given below. The "Curve Version" column gives
different versions of the completeness versus flux files (the `sn.use` files). These `sn.use` files also contain 50% 
completeness values, which some models use to derive the scaling between noise in the sensitivity cubes and 
the flux at 50% completeness. This scaling is needed to rescale the curves from the `sn.use` files as a function
of position and wavelength. Most users should just use the latest version that's usable for the data release 
they are using, the "Notes" column is mainly for the reference of the developers.

| Model Name | Curve Version | Data releases       | Noise to 50% flux scaling                      |   Notes                         | 
| ---------- | ------------- | --------------------| ---------------------------------------------- | -------------------             | 
| hdr2pt1pt1 |       v1      |      2.1.1-onwards  | `cal_on_karl_results.py` run on curves v0      | Near identical results to v1 for S/N < 6 or so  |
| hdr2pt1pt3 |       v1      |      2.1.1-onwards  | `cal_on_karl_results.py` run on curves v1      |                                 |
|    v1      |       v1      |      2.1.1-onwards  |  KG derived scaling (`snlist` file)            | Doesn't interpolate `sn.use` if completeness in lower bin is zero |
|    v1.1    |       v1      |      2.1.1-onwards  |  KG derived scaling (`snlist` file)            | v1 with interpolated sens. cubes |

These are the different flux versus completeness curve versions (the `sn.use` files)

| Curve Version | Notes                                           |
| ------------- | ----------------------------------------------  |
|      v0       | check.f bug didn't apply cuts to sims properly  |
|      v1       | check.f bug fixed                               |

## Command line tools

This package contains multiple command line tools. Two of the most important
deal with adding and extracting sensitivity cubes. To add a FITS file sensitivity cube(s)
to a HDF5 file, one can use the `add_sensitivity_cube_to_hdf5` command, like so

```
add_sensitivity_cube_to_hdf5 20181203v013_multi_*.fits 20181203v013_sensitivity_cubes.h5
```

This command extracts the shot and IFU information from the filename, using the 
regex capture groups specified with the optional `--regex` flag. The default
`--regex` flag should work with the standard HETDEX file-naming convention. To extract
the sensivitity as a FITs file use

`extract_sensitivity_cube 20181203v013_sensitivity_cubes.h5 087 outputfile.fits`

The second argument specifies the IFU slot you want to extract from the
cube. If there are multiple shots in a cube, you can specify which shot
you want to extract from using the `--datevshot` flag, e.g.

```
extract_sensitivity_cube --datevshot 20181203v013 20181203v013_sensitivity_cubes.h5 087 outputfile2.fits
```

You can also extract the biweight flux limits of all IFUs and shots, collapsed in ra and dec,
in an HDF5 with the following command

```
biweight_fluxlims_hdf5 20181203v013_sensitivity_cubes.h5
```

using the ``--wl`` flag you can specify what wavelength you want. It is also possible to output these limits
to a file, use the ``-h`` option on this command for details.

### Getting the flux limit in Python - the python API

To access the flux limit an application programming interface (API) is provided for
use in Python. An iPython notebook is provided in this repository which demonstrates
how to do this [here](../../notebooks/04-Getting_Flux_Limits_from_the_HDF5_Files.ipynb)


