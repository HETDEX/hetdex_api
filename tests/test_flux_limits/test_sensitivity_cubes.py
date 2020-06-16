"""

Test the Sensitivity Cube class

AUTHOR: Daniel Farrow (MPE)


"""
import pytest
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube

@pytest.fixture(scope="module")
def sensitivity_cube(datadir):
    """ A sensitivity cube """
    filename = datadir.join("test_sensitivity_cube.fits").strpath
    wavelengths = [3500.0, 5500.0]
    alphas = [-3.5, -3.5]

    return SensitivityCube.from_file(filename, wavelengths, alphas, nsigma=1.0)

@pytest.fixture(scope="module")
def modified_sensitivity_cube(datadir):
    """ A sensitivity cube with a conversion polynomial set"""
    filename = datadir.join("test_sensitivity_cube.fits").strpath
    wavelengths = [3500.0, 5500.0]
    alphas = [-3.5, -3.5]
    p = [2.0, 0.0]

    return SensitivityCube.from_file(filename, wavelengths, alphas, 
                                     conversion_poly=p, nsigma=1.0)


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
    ratio = scube.sigmas/scube_corr1.sigmas

    assert ratio == pytest.approx(aper_corr)


def test_snr_conversion(modified_sensitivity_cube, sensitivity_cube):
    """
    Test that the conversion between the two SNR values
    is being applied
    """
    f501 = sensitivity_cube.get_f50(2e-16,161.4201, 50.8822, 4.0)
    f502 = modified_sensitivity_cube.get_f50(2e-16,161.4201, 50.8822, 4.0)

    assert f501 == pytest.approx(f502/2.0)


def test_completeness_func_hdr1(sensitivity_cube):
    """
    Test that a value is returned
    """
    c = sensitivity_cube.return_completeness_hdr1(2e-16,161.4201, 50.8822, 3478)
    assert c == pytest.approx(0.20411406725738124)

@pytest.mark.parametrize("sncut, expected", [(6.0, 0.20411406725738124),
                                             (5.0, 0.47937941516439253)])
def test_completeness_func(sensitivity_cube, sncut, expected):
    """
    Test that a value is returned
    """
    c = sensitivity_cube.return_completeness(2e-16, 161.4201, 50.8822, 3478, sncut)
    assert c == pytest.approx(expected)


def test_flux_calib_corr(datadir, plot=False):
    """
    Test the handling of the flux calibration
    """

    filename = datadir.join("test_sensitivity_cube.fits").strpath
    #filename = datadir.join("test.fits").strpath
    fcalib_corr = datadir.join("polyvals.txt").strpath
    wavelengths = [3500.0, 5500.0]
    alphas = [-3.5, -3.5]

    scube = SensitivityCube.from_file(filename, wavelengths, alphas)
    scube_corr1 = SensitivityCube.from_file(filename, wavelengths, alphas)
    scube_corr1.apply_flux_recalibration(1.0, flux_calib_correction_file=fcalib_corr) 

    # The aper_corr should be multiplied though in the cubes
    ratio = scube.sigmas/scube_corr1.sigmas

    if plot:
        import matplotlib.pyplot as plt
        ra, dec, wls = scube.wcs.wcs_pix2world(0, 0, range(scube.sigmas.shape[0]), 0)
        plt.plot(wls, (scube.f50vals[:,10,10] - scube_corr1.sigmas[:,10,10])/scube.sigmas[:,10,10])
        plt.show()

    # Check something happened
    assert ratio != pytest.approx(1.0)
