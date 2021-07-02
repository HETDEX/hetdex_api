import pytest
from os.path import join
from hetdex_api.config import HDRconfig
from hetdex_api.flux_limits.flim_models import read_karl_file, SimulationInterpolator, SingleSNSimulationInterpolator


@pytest.mark.parametrize("dont_interp_to_zero", [True, False])
@pytest.mark.parametrize("release", ["hdr2.1"])
@pytest.mark.parametrize("sncut", [4.8, 5.0, 6.0])
def test_simulation_interpolator(release, sncut, dont_interp_to_zero):
    """ 
    Compare the interpolator to the file directly 
    and ensure they are the same
    """

    conf = HDRconfig(survey=release)
    fdir = join(conf.flim_sim_completeness, "curves_v1")
    sinterp = SimulationInterpolator(fdir, dont_interp_to_zero,
                                     wl_collapse = False,
                                     cmax = None)
 
    ffile = fdir + "/sn{:2.1f}.use".format(sncut)
    waves, f50s, compl_curves, fluxes = read_karl_file(ffile)

    for wl, f50, compl_curve in zip(waves, f50s, compl_curves):
        model = sinterp(fluxes, f50, wl, sncut)
        print(wl, f50, (model == pytest.approx(compl_curve, rel=1e-2, abs=1e-2)))
        assert model == pytest.approx(compl_curve, rel=1e-2, abs=1e-2)

@pytest.mark.parametrize("dont_interp_to_zero", [True, False])
def test_interpolation_toward_zero(datadir, dont_interp_to_zero):

    wave = 3957.0
    flux = 3.999
    f50 = 9.0

    fdir = join(datadir, "sntest.use")
    sinterp = SingleSNSimulationInterpolator(fdir, dont_interp_to_zero,
                                             wl_collapse = False,
                                             cmax = None)

    print(sinterp(flux, f50, wave), dont_interp_to_zero)
    if dont_interp_to_zero:
        assert sinterp(flux, f50, wave) < 1e-20
    else:
        assert sinterp(flux, f50, wave) > 0.0



 
