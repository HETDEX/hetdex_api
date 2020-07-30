# Sensitivity Cubes

## Overview 

The flux limits in HETDEX are stored in ra, dec, wavelength datacubes ("sensitivity cubes"). The 
value stored is the inverse of the PSF weighted 1-sigma noise, which can be used to
compute flux limits. This is how we can store the flux limit as a function of position and wavelength.

To characterise the fraction of sources detected as a function of flux, a parameterisation
from [Fleming et al 1995](http://adsabs.harvard.edu/abs/1995AJ....109.1044F) is used. One input
to this paramterisation are the flux limits stored in the datacubes. The second and
final parameter in this function controls how quickly the detection rate falls off with 
flux. This is currently set to a fixed value, calibrated from simulations.


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


