"""

Test the shot_sensitivity module

AUTHOR: Daniel Farrow

"""
import pytest
from numpy import ones, arange, ones, abs, isfinite, nanmax
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
                           ),
                           (178.700104, 50.90477, 3785.28, 1e-16,  None, 0.00392, "v4"),
                           (178.620972, 50.787903, 5243.32, 7e-17, None, 0.009204, "v4"),
                           (178.617676, 50.788761, 4509.88, 9e-17, None, 0.006111, "v4"),
                           (178.539658, 50.743679, 4163.34, 2e-16, None, 0.011308, "v4"),
                           (
                            [178.880859, 178.88652], [50.8312, 50.827698],
                            [4602.69, 3863.54], [8e-17, 9e-17], None,
                            [0.1499, 0.03193], 
                            "v4"
                           ),
                          ] 
                         )
def test_completeness(ra, dec, wave, flux, lw, exptd, flim_model):

    s = ShotSensitivity("20190316v019", release="hdr2.1", 
                        flim_model=flim_model, sclean_bad = False)
 
    c = s.return_completeness(flux/ext(wave), ra, dec, wave, 5.0,
                              linewidth=lw)

    assert pytest.approx(c, rel=0.005) == exptd



@pytest.mark.parametrize("ra, dec, sclean",
                          [
                            (201.5498985, 51.82207577, True),
                            (201.5498985, 51.82207577, False),

                          ]                           
                        )
def test_wave_none_mode(ra, dec, sclean):
    """ Compare passing wave=None to passing waves """

    s = ShotSensitivity("20190209v024", 
                        sclean_bad = sclean, 
                        log_level="INFO")

    all_sigmas, apcor1, mask  = s.get_f50(ra, dec, None, 5.0,
                                          direct_sigmas = True)


    print("All sigmas done!")

    waves = s.extractor.get_wave()
     
    # If wave is an array ra/dec has to be as well
    from_single_wave, apcor2 = s.get_f50(ra*ones(len(waves)), dec*ones(len(waves)), 
                                         waves, 5.0, 
                                         direct_sigmas=True)

    
    diff = (all_sigmas[0] - from_single_wave)/all_sigmas[0]
    diff_apcor = (apcor1[0] - apcor2)/apcor1[0]

    print(all_sigmas[0])
    print(from_single_wave)
    print(min(from_single_wave), max(from_single_wave))

    # test if good to within 0.05%
    assert nanmax(abs(diff_apcor)) < 5e-4
    assert nanmax(abs(diff)) < 5e-4


@pytest.mark.parametrize("ras, decs",
                         [
                          (
                           [229.95960512542143, 229.74236696267425],
                           [53.8360578364009, 54.0562276323938]
                          )
                         ]
                        )
def test_wave_none_array(ras, decs):
    """ 
    As test_wave_none_mode but multiple 
    ra/dec at once

    """
    s = ShotSensitivity("20200423v019")
    
    # All at once
    f50_arr, apcor_arr, mask = s.get_f50(ras, decs, None, 1.0,
                                         direct_sigmas = True)

    waves = s.extractor.get_wave()

    for i, (r, d) in enumerate(zip(ras, decs), 0):

        f50, apcor = s.get_f50(r*ones(len(waves)), d*ones(len(waves)), 
                               waves, 1.0,
                               direct_sigmas = True)

        diff = (f50 - f50_arr[i, :])/f50_arr[i, :]
        diff_apcor = (apcor - apcor_arr[i, :])/apcor_arr[i, :]
        
        print(apcor)
        print(apcor_arr[i, :])

        # test if good to within 0.05%
        assert nanmax(abs(diff_apcor)) < 5e-4
        assert nanmax(abs(diff)) < 5e-4


def test_chip_gap_mask():

    s = ShotSensitivity("20200129v029", release="hdr2.1", 
                        flim_model="v2")
 
    c, norm = s.get_f50(201.1692954, 50.7903255, 4512.0, 5.5)

    print(c)
    assert c > 998
 
