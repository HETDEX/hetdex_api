"""

Tests for the extract module

"""
import pytest
from astropy.coordinates import SkyCoord
from hetdex_api.extract import Extract


@pytest.fixture(scope="module")
def extract_class():
    """
    Set up an Extract class and load
    a particular shotid
    """
    e = Extract()
    e.load_shot("20190201012", survey="hdr2.1")
    return e


@pytest.mark.parametrize("fwhm, ra, dec, exptd_sum, exptd_1, exptd_2", 
                         [(1.1, 163.353, 51.6549, 950.0812383206099, 
                           0.037616255398805365, 0.09509532139159248),
                          (2.5, 163.400, 51.8978, 910.8985995435256,
                           2.7117054674636605, 7.310208323356028),
                          (3.0, 163.367, 51.647, 852.6199603368874,
                           0.8423417542640232, 1.6824348505231255)])
 
def test_build_weights(fwhm, ra, dec, exptd_sum, exptd_1, exptd_2, 
                       extract_class):
    """
    Test the build weights method against a previous
    run of it we think worked correctly

    """

    psf = extract_class.moffat_psf(fwhm, 10.5, 0.25)
    c = SkyCoord(ra=ra, dec=dec, unit="deg")
    info_result = extract_class.get_fiberinfo_for_coord(
                           c,
                           radius=3.5,
                           ffsky=False,
                           return_fiber_info=True,
                           )
    ifux, ifuy, xc, yc, ra, dec, data, error, mask, fiberid, \
               multiframe = info_result

    
    weights = extract_class.build_weights(xc, yc, ifux, ifuy, psf)
 
    assert pytest.approx(weights[5, 167]*100, rel=1e-5) == exptd_1
    assert pytest.approx(weights[2, 801]*100, rel=1e-5) == exptd_2
    assert pytest.approx(sum(weights.ravel()), rel=1e-5) == exptd_sum


