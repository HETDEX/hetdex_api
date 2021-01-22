"""

Test the Sensitivity Cube class

AUTHOR: Daniel Farrow (MPE)


"""
import pytest
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube


@pytest.mark.parametrize("aper_corr", [0.25, 0.33, 0.7])
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

    # The aper_corr should be multiplied though in the cubes
    ratio = scube.sigmas.filled()/scube_corr1.sigmas.filled()

    assert ratio == pytest.approx(aper_corr)


@pytest.mark.parametrize("model, flux, sncut, expected", [
                                                    ("hdr1", 2e-16, 6.0, 0.20411406725738124),
                                                    ("hdr1", 2e-16, 5.0, 0.47937941516439253),
                                                    ("hdr2pt1", 4e-16, 6.0, 0.2115161903693431),
                                                    ("hdr2pt1", 4e-16, 5.0, 0.4961144900645637),
                                                    ("hdr2pt1pt1", 4e-16, 6.0, 0.52132582),
                                                    ("hdr2pt1pt1", 4e-16, 5.0, 0.64842545)
                                                    ])
def test_completeness_func(datadir, flux, model, sncut, expected):
    """
    Test that a value is returned
    """
    filename = datadir.join("test_sensitivity_cube.fits").strpath
    wavelengths = [3500.0, 5500.0]
    alphas = [-3.5, -3.5]

    if model == "hdr1":
        scube = SensitivityCube.from_file(filename, wavelengths, alphas, nsigma=1.0,
                                          flim_model=model, aper_corr=None)
    else:
        scube = SensitivityCube.from_file(filename, wavelengths, alphas, nsigma=1.0,
                                          flim_model=model)
 
    c = scube.return_completeness(flux, 161.4201, 50.8822, 3478, sncut)
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
    ratio = scube.sigmas.filled()/scube_corr1.sigmas.filled()

    if plot:
        import matplotlib.pyplot as plt
        ra, dec, wls = scube.wcs.wcs_pix2world(0, 0, range(scube.sigmas.shape[0]), 0)
        plt.plot(wls, (scube.f50vals[:,10,10] - scube_corr1.sigmas[:,10,10])/scube.sigmas[:,10,10])
        plt.show()

    # Check something happened
    assert ratio != pytest.approx(1.0)

