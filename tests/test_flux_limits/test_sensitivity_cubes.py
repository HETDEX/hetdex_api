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
