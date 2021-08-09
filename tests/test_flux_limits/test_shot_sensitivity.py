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
def test_get_f50(ra, dec, wave, exptd, shot_sensitivity):
    f50 = shot_sensitivity.get_f50(ra, dec, wave, 5.0, 
                                   direct_sigmas=True)

    assert pytest.approx(exptd, rel=0.2) == f50*1e17


@pytest.mark.parametrize("ra, dec, wave, flux, exptd, flim_model",
                          [
                          (178.700104, 50.90477, 3785.28, 1e-16,  0.003964, "v1"),
                          (178.620972, 50.787903, 5243.32, 7e-17, 0.009348, "v1"),
                          (178.617676, 50.788761, 4509.88, 9e-17, 0.006173, "v1"),
                          (178.539658, 50.743679, 4163.34, 2e-16, 0.011419, "v1"),
                          (
                           [178.880859, 178.88652], [50.8312, 50.827698],
                           [4602.69, 3863.54], [8e-17, 9e-17], [0.15225, 0.03251], 
                           "v1"
                          )
                         ])
def test_completeness(ra, dec, wave, flux, exptd, flim_model):

    s = ShotSensitivity("20190316v019", release="hdr2.1", 
                        flim_model=flim_model)
 
    c = s.return_completeness(flux, ra, dec, wave, 5.0)

    assert pytest.approx(c, rel=0.005) == exptd

