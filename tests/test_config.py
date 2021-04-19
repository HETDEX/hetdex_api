"""

"""
import pytest
from os.path import exists
from hetdex_api.config import HDRconfig


#@pytest.mark.parametrize("release", ["hdr1", "hdr2", "hdr2.1"])

@pytest.mark.parametrize("release", ["hdr2.1"])
@pytest.mark.parametrize("attribute", ["red_dir", "detect_dir",
                                       "elix_dir", "surveyh5",
                                       "detecth5", "badamp", 
                                       "baddetect"])
def test_path_exists(attribute, release):
    """
    Test that some of the 
    important paths in config exist. 
    """
    conf=HDRconfig(survey=release)
    assert exists(getattr(conf, attribute))
