import pytest
from hetdex_api.config import HDRconfig
from hetdex_api.flux_limits.flim_models import read_karl_file, SimulationInterpolator


@pytest.mark.parametrize("release", ["hdr2.1"])
@pytest.mark.parametrize("sncut", [4.8, 5.0])
def test_simulation_interpolator(release, sncut):
    """ 
    Compare the interpolator to the file directly 
    and ensure they are the same
    """

    conf = HDRconfig(survey=release)
    fdir = conf.flim_sim_completeness
    sinterp = SimulationInterpolator(fdir,
                                     wl_collapse = False,
                                     cmax = None)
 
    ffile = fdir + "/sn{:2.1f}.use".format(sncut)
    waves, f50s, compl_curves, fluxes = read_karl_file(ffile)

    for wl, f50, compl_curve in zip(waves, f50s, compl_curves):
        model = sinterp(fluxes, f50, wl, sncut)
        print(wl, f50, (model == pytest.approx(compl_curve, rel=1e-2, abs=1e-2)))
        assert model == pytest.approx(compl_curve, rel=1e-2, abs=1e-2)

