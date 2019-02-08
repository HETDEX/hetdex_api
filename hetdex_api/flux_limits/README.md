# Sensitivity Cubes

## Overview 

The flux limits in HETDEX are stored in ra, dec, wavelength datacubes. The flux value stored
is the 6-sigma point source detection limit, which roughly corresponds to the flux at
which 50% of the sources are detected. This is how we can store the flux limit as a 
function of position and wavelength.

To characterise the fraction of sources detected as a function of flux, a parameterisation
from [Fleming et al 1995](http://adsabs.harvard.edu/abs/1995AJ....109.1044F) is used. One input
to this paramterisation are the flux limits stored in the datacubes. The second and
final parameter in this function controls how quickly the detection rate falls off with 
flux. This is currently set to a fixed value, calibrated from simulations.

## Using the software in this package

This package defines a ``SensitivityCube`` class, which handles the cubes. 

