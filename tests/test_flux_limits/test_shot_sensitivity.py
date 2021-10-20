"""

Test the shot_sensitivity module

AUTHOR: Daniel Farrow

"""
import pytest
from hetdex_api.flux_limits.shot_sensitivity import ShotSensitivity
from hetdex_api.extinction import get_2pt1_extinction_fix

# Fix for hdr2.1 and earlier
ext = get_2pt1_extinction_fix()

@pytest.fixture(scope="module")
def shot_sensitivity():
    s = ShotSensitivity("20190316v019", release="hdr2.1", 
                        flim_model="v2")
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

    exptd = exptd/ext(wave) 
    f50, norm = shot_sensitivity.get_f50(ra, dec, wave, 5.0, 
                                         direct_sigmas=True)

    assert pytest.approx(exptd, rel=0.2) == f50*1e17

@pytest.mark.parametrize("ra, dec, wave, linewidth, exptd",
                         [
                          (178.700104, 50.90477, 3785.28, 2.2, 30.31),
                          (178.700104, 50.90477, 3785.28, None, 30.31),
                          (178.620972, 50.787903, 5243.32, None, 17.84), 
                          (178.620972, 50.787903, 5243.32, 8.5, 44.99),
                          (178.620972, 50.787903, 5243.32, 4.5, 26.72),
                          (
                              [178.620972, 178.620972],
                              [50.787903, 50.787903],
                              [5243.32, 5243.32],
                              [2.5, 8.5],
                              [18.24, 44.99]

                          )
                          ]
                         )
def test_get_f50_lw(ra, dec, wave, linewidth, exptd, shot_sensitivity):
    f50, norm = shot_sensitivity.get_f50(ra, dec, wave, 5.0, 
                                         linewidth = linewidth)

    exptd = exptd/ext(wave) 
    assert pytest.approx(exptd, rel=0.02) == f50*1e17


@pytest.mark.parametrize("ra, dec, wave, flux, lw, exptd, flim_model",
                          [
                           (178.700104, 50.90477, 3785.28, 1e-16,  None, 0.00392, "v1"),
                           (178.620972, 50.787903, 5243.32, 7e-17, None, 0.009204, "v1"),
                           (178.617676, 50.788761, 4509.88, 9e-17, None, 0.006111, "v1"),
                           (178.539658, 50.743679, 4163.34, 2e-16, None, 0.011308, "v1"),
                           (
                            [178.880859, 178.88652], [50.8312, 50.827698],
                            [4602.69, 3863.54], [8e-17, 9e-17], None,
                            [0.1499, 0.03193], 
                            "v1"
                           ),
                           (178.700104, 50.90477, 3785.28, 1e-16,  None, 0.003580, "v2"),
                           (178.620972, 50.787903, 5243.32, 7e-17, None, 0.007610, "v2"),
                           (178.620972, 50.787903, 5243.32, 7e-17, 2.2,  0.007610, "v2"),
                           (178.620972, 50.787903, 5243.32, 7e-17, 6.6,  0.000603, "v2"),
                           (178.617676, 50.788761, 4509.88, 9e-17, None, 0.005033, "v2"),
                           (178.539658, 50.743679, 4163.34, 2e-16, None, 0.009651, "v2"),
                           (
                            [178.880859, 178.88652], [50.8312, 50.827698],
                            [4602.69, 3863.54], [8e-17, 9e-17], [1.2, 3.4],
                            [0.11211, 0.0086177], 
                            "v2"
                           )
                          ] 
                         )
def test_completeness(ra, dec, wave, flux, lw, exptd, flim_model):

    s = ShotSensitivity("20190316v019", release="hdr2.1", 
                        flim_model=flim_model, sclean_bad = False)
 
    c = s.return_completeness(flux/ext(wave), ra, dec, wave, 5.0,
                              linewidth=lw)

    assert pytest.approx(c, rel=0.005) == exptd


def test_chip_gap_mask():

    s = ShotSensitivity("20200129v029", release="hdr2.1", 
                        flim_model="v2")
 
    c, norm = s.get_f50(201.1692954, 50.7903255, 4512.0, 5.5)

    print(c)
    assert c > 998
 
