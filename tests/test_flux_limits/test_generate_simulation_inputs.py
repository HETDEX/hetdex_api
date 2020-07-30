"""
Test the generate simulation inputs module

AUTHOR: Daniel Farrow (MPE)

"""
import pytest
from hetdex_api.survey import Survey
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube
from hetdex_api.flux_limits.generate_simulation_inputs import (create_sensitivity_cube_from_astrom,
                                                               generate_sencube_hdf)



@pytest.fixture
def mock_hdf_cube(datadir):
    """ Generate a mock HDF cube for a given shot """
    datevobs = "20190201v019"
    survey_obj = Survey("hdr2")
    shot = survey_obj[survey_obj.datevobs == datevobs]
    sencube_hdf = generate_sencube_hdf(datevobs, shot.ra[0], 
                                       shot.dec[0], shot.pa[0], datadir, 
                                       31, 31, 1036, 62.0)

    return sencube_hdf



@pytest.mark.parametrize("x, y, xexpt, yexpt", [(25.0, 10.0, 209.60643351, 51.51866332),
                                                (10.0, 25.0, 209.61443495, 51.52934529),
                                                (-10.0, -25.0, 209.5663673, 51.53027137)])
def test_sensitivity_cube_astrometry(mock_hdf_cube, x, y, xexpt, yexpt):
    """ Test the astrometry in the header is same as in Karl's cube """

    scube = mock_hdf_cube.extract_ifu_sensitivity_cube("ifuslot_052")

    ix, iy, iz = scube.wcs.wcs_pix2world(x, y, 0, 0)
    print(scube.header)
    print(x, y, ix, iy, iz, xexpt, yexpt)

    # reference cubes not quite on IFU center
    # test only good to a few arceconds
    assert ix == pytest.approx(xexpt, abs=1e-3)
    assert iy == pytest.approx(yexpt, abs=1e-3)

