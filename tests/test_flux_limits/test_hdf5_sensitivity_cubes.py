"""

Tests for converting between HDF5
and FITs sensitivity cubes

AUTHOR: Daniel Farrow

"""

from shutil import copy
from os.path import isfile
import pytest
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import (SensitivityCubeHDF5Container, 
                                                           add_sensitivity_cube_to_hdf5,
                                                           extract_sensitivity_cube,
                                                           return_sensitivity_hdf_path) 
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube

@pytest.fixture(scope="function")
def hdf5_container(tmpdir):
    """ HDF5 we can write test data to """
    filename = tmpdir.join("test.h5").strpath
    hdcon = SensitivityCubeHDF5Container(filename, mode="w")

    # Clever trick to close the file when we're done with it 
    yield hdcon
    hdcon.close()
    
@pytest.fixture(scope="module")
def sensitivity_cube(datadir):
    """ A sensitivity cube to add or read"""
    filename = datadir.join("test_sensitivity_cube.fits").strpath
    wavelengths = [3500.0, 5500.0]
    alphas = [-3.5, -3.5]

    return SensitivityCube.from_file(filename, wavelengths, alphas)

@pytest.mark.parametrize("use_with", [True, False])
def test_hdf5_create_and_write(tmpdir, use_with):
    """
    Test writing the HDF5 file to a 
    file
    """

    filename = tmpdir.join("test.h5").strpath

    # Test with statement
    if use_with:
        with SensitivityCubeHDF5Container(filename, mode="w"):
            pass
    else:
        # Test explicitly closing
        hdcon = SensitivityCubeHDF5Container(filename, mode="w")
        hdcon.close()

    # Can we open it again?
    hdcon2 = SensitivityCubeHDF5Container(filename, mode="r")
    hdcon2.close()
    

@pytest.mark.parametrize('flush', [True, False])
def test_hdf5_add_extract(hdf5_container, sensitivity_cube, flush):
    """
    Test the adding and extracting of 
    sensitivity cubes into HDF5 
    """

    # Try adding it
    hdf5_container.add_sensitivity_cube("shot_19690720v11", "ifuslot_000", 
                                        sensitivity_cube, flush=flush)

    # Try getting it back again
    scube = hdf5_container.extract_ifu_sensitivity_cube("ifuslot_000")

    # Check it's the same (multiply by big number of approx thinks they're the same)
    vals_new = scube.sigmas.flatten()*1e16
    vals_old = sensitivity_cube.sigmas.flatten()*1e16

    assert vals_new == pytest.approx(vals_old)


def test_add_to_hdf5_cmd(tmpdir, datadir):
    """
    Test the command line tool to add
    to HDF5 files
    """
    filename_original = datadir.join("test_sensitivity_cube.fits").strpath
  
    # Make some files to input
    scube_fn1 = tmpdir.join("20181203v013_multi_324_062_055.fits").strpath
    scube_fn2 = tmpdir.join("20181203v013_multi_013_103_019.fits").strpath
    copy(filename_original, scube_fn1)
    copy(filename_original, scube_fn2)

    output = tmpdir.join("test_output.h5").strpath

    # Run with command line arguments passed
    args = ["--regex", ".*(2[0-9]{7}v[0-9]{3})_multi_[0-9]{3}_([0-9]{3})",
            scube_fn1, scube_fn2, output] 
    add_sensitivity_cube_to_hdf5(args=args)
 
    assert isfile(output)

@pytest.mark.parametrize("datevshot", [None, "20181203v013"])
def test_extract_sensitivity_cube(tmpdir, datadir, datevshot):
    """
    Test the command to extract sensitivity cubes
    from HDF5 containers
    """
    
    h5fn = datadir.join("test_hdf.h5").strpath
    outfn = tmpdir.join("test.fits").strpath 
    
    if datevshot:
        args = [h5fn, "--datevshot", datevshot, "063", outfn]
    else:
        args = [h5fn, "063", outfn] 

    extract_sensitivity_cube(args=args)
     
    assert isfile(outfn)


def test_flim_model(datadir):
    """
    Test that the flux limit model is 
    being passed to the sensitivity cubes
    """
    # ifuslot_063
    filename = datadir.join("test_hdf.h5").strpath
    hdcon1 = SensitivityCubeHDF5Container(filename, flim_model="hdr1")
    hdcon2 = SensitivityCubeHDF5Container(filename, flim_model="hdr2pt1")

    scube1 = hdcon1.extract_ifu_sensitivity_cube("ifuslot_063")
    scube2 = hdcon2.extract_ifu_sensitivity_cube("ifuslot_063")

    s1 = scube1.get_f50(161.4201, 50.8822, 3470.0, 5.5)
    s2 = scube2.get_f50(161.4201, 50.8822, 3470.0, 5.5)

    print(s1)
    # if different models passed should be different
    assert abs(s1 - s2) > 1e-19


@pytest.mark.parametrize("datevshot", ["20181203v013", "20191023v024"])
def test_return_sensitivity_hdf_path(datevshot):
    """ Test the correct path is returned """
    fn = return_sensitivity_hdf_path(datevshot)
    assert isfile(fn)    






