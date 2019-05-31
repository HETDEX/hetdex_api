"""

Test the Sensitivity Cube class

AUTHOR: Daniel Farrow (MPE)


"""
import pytest
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube

@pytest.mark.parametrize("aper_corr", [0.25, 0.33, None])
def test_aper_corr(datadir, aper_corr):
    """
    Test the handling of the aperture
    correction
    """
    filename = datadir.join("test_sensitivity_cube.fits").strpath
    wavelengths = [3500.0, 5500.0]
    alphas = [-3.5, -3.5]

    scube = SensitivityCube.from_file(filename, wavelengths, alphas, 
                                      aper_corr=aper_corr)

    scube_corr1 = SensitivityCube.from_file(filename, wavelengths, alphas, 
                                            aper_corr=1.0)

    # If None, should be taken from the header
    if not aper_corr:
        aper_corr = scube.header["APCOR"]

    # The aper_corr should be multiplied though in the cubes
    ratio = scube.f50vals/scube_corr1.f50vals

    assert ratio == pytest.approx(aper_corr)


def test_flux_calib_corr(datadir, plot=True):
    """
    Test the handling of the flux calibration
    """

    #filename = datadir.join("test_sensitivity_cube.fits").strpath

    filename = datadir.join("test.fits").strpath
    fcalib_corr = datadir.join("polyvals.txt").strpath
    wavelengths = [3500.0, 5500.0]
    alphas = [-3.5, -3.5]

    scube = SensitivityCube.from_file(filename, wavelengths, alphas)
    scube_corr1 = SensitivityCube.from_file(filename, wavelengths, alphas)
    scube_corr1.apply_flux_recalibration(fcalib_corr) 

    # The aper_corr should be multiplied though in the cubes
    ratio = scube.f50vals/scube_corr1.f50vals

    if plot:
        import matplotlib.pyplot as plt
        ra, dec, wls = scube.wcs.wcs_pix2world(0, 0, range(scube.f50vals.shape[0]), 0)
        plt.plot(wls, (scube.f50vals[:,10,10] - scube_corr1.f50vals[:,10,10])/scube.f50vals[:,10,10])
        plt.show()

    # Check something happened
    assert ratio != pytest.approx(1.0)
