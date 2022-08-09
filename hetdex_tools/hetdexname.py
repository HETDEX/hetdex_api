from astropy.coordinates import SkyCoord
import astropy.units as u

def get_source_name(ra, dec):
    """
    Convert ra, dec coordinates in degrees to a IAU-style object name.

    Example
    =======
    source_name = get_source_name(123.4512, 65.2341)
    
    """
    coord = SkyCoord(ra * u.deg, dec * u.deg)

    return "HETDEX_J{0}{1}".format(
        coord.ra.to_string(unit=u.hourangle, sep="", precision=2, pad=True),
        coord.dec.to_string(sep="", precision=1, alwayssign=True, pad=True),
    )
