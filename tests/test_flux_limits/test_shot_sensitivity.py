"""

Test the shot_sensitivity module

AUTHOR: Daniel Farrow

"""
import pytest
from hetdex_api.flux_limits.shot_sensitivity import ShotSensitivity


@pytest.fixture(scope="module")
def shot_sensitivity():
    s = ShotSensitivity("20190316v019", release="hdr2.1", 
                        flim_model="v1")
    return s

def test_extract_cube(shot_sensitivity):
    """ Test extracting a cube returns something """
    c = shot_sensitivity.extract_ifu_sensitivity_cube("ifuslot_094", nx=3, ny=3,
                                                      ifusize=62.0,
                                                      generate_sigma_array=True)

    assert any(c.sigmas.ravel() > 0)


# Reference values from Karl Gebhardt and Laurel Weiss
# test simulations
@pytest.mark.parametrize("ra, dec, wave, exptd",
                         [
                          (178.700104, 50.90477, 3785.28, 6.212),
                          (178.620972, 50.787903, 5243.32, 3.297),
                          (178.617676, 50.788761, 4509.88, 4.466),
                          (178.539658, 50.743679, 4163.34, 9.597),
                          (
                           [178.880859, 178.88652], [50.8312, 50.827698],
                           [4602.69, 3863.54], [2.319, 3.718]
                          )
                         ])
def test_get_f50(ra, dec, wave, shot_sensitivity, exptd):
    f50 = shot_sensitivity.get_f50(ra, dec, wave, 
                                   direct_sigmas=True)

    assert pytest.approx(exptd, rel=0.4) == f50*1e17


@pytest.mark.parametrize("ra, dec, wave, flux, exptd",
                          [
                          (178.700104, 50.90477, 3785.28, 7e-17, 0.7589),
                          (178.620972, 50.787903, 5243.32, 3e-17, 0.6118),
                          (178.617676, 50.788761, 4509.88, 2e-17, 0.0293),
                          (178.539658, 50.743679, 4163.34, 1e-16, 0.6652),
                          (
                           [178.880859, 178.88652], [50.8312, 50.827698],
                           [4602.69, 3863.54], [2e-17, 7e-17], [0.5351, 0.9812]
                          )
                         ])
def test_completeness(shot_sensitivity, ra, dec, wave, flux, exptd):

    c = shot_sensitivity.return_completeness(flux, ra, dec, wave, 5.0)

    assert pytest.approx(c, rel=0.005) == exptd

