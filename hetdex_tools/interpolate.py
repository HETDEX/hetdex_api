import os.path as op
import numpy as np
import tables
import tables as tb

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.table import Table

from hetdex_api.config import HDRconfig
from hetdex_api.extract import Extract

import datetime

try:  # using HDRconfig
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
    config = HDRconfig(LATEST_HDR_NAME)
    CONFIG_HDR2 = HDRconfig("hdr2.1")
    CONFIG_HDR3 = HDRconfig("hdr3")
    CONFIG_HDR4 = HDRconfig("hdr4")
    CONFIG_HDR5 = HDRconfig("hdr5")

except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr5"

OPEN_DET_FILE = None
DET_HANDLE = None
DET_FILE = None

current_hdr = LATEST_HDR_NAME
surveyh5 = None


def make_narrowband_image(
    detectid=None,
    coords=None,
    shotid=None,
    pixscale=None,
    imsize=None,
    wave_range=None,
    convolve_image=False,
    ffsky=False,
    ffsky_rescor=False,
    subcont=False,
    dcont=50.0,
    include_error=False,
    include_bitmask=False,
    survey=LATEST_HDR_NAME,
    extract_class=None,
    fiber_flux_offset=None,
    interp_kind="linear",
    apply_mask=False,
    mask_options=None,
    fill_value=0.0,
    include_grid=False,
    shot_h5=None,
):
    """
    Function to make narrowband image from either a detectid or from a
    coordinate/shotid combination.

    Paramaters
    ----------
    detectid: int
        detectid from the continuum or lines catalog. Default is
        None. Provide a coords/shotid combo if this isn't given
    coords: SkyCoords object
        coordinates to define the centre of the data cube
    pixscale: astropy angle quantity
         plate scale
    imsize: astropy angle quantity
        image size
    wave_range: list or None
        start and stop value for the wavelength range in Angstrom.
        If not given, the detectid linewidth is used
    convolve_image: bool
        option to convolve image with shotid seeing
    ffsky: bool
        option to use full frame calibrated fibers. Default is
        True.
    ffsky_rescor: bool
        option to use full frame calibrated fibers with residual
        correction. Default False
    subcont: bool
        option to subtract continuum. Default is False. This
        will measure the continuum 50AA below and above the
        input wave_range
    dcont
        width in angstrom to measure the continuum. Default is to
        measure 50 AA wide regions on either side of the line
    include_error bool
        option to include error array
    extract   Extract class object
        option to include a preloaded Extract class object.
        Default is to load extract class according to detection info
    fiber_flux_offset: 1036 array
        array of values in units of 10**-17 ergs/s/cm2/AA to add
        to each fiber spectrum used in the extraction. Defaults
        to None
    interp_kind: str
        Kind of interpolation to pixelated grid from fiber intensity.
        Options are 'linear', 'cubic', 'nearest'. Default is linear.
    apply_mask: bool
        Apply HETDEX fiber mask model. This will mask all fibers contributing
        to the spectral extraction before summation. Masked in place as NaN
    mask_options
        string or array of strings as options to select to mask. Default None
        will select all flags. Set this to 'BITMASK' to return the full bitmask array.
        Options are 'MAIN', 'FTF', 'CHI2FIB', 'BADPIX', 'BADAMP', 'LARGEGAL', 'METEOR',
        'BADSHOT', 'THROUGHPUT', 'BADFIB', 'SAT'
        If BITMASK appears as any element in the list, it overrides all others
        and returns the full bitmask array.
    include_bitmask
        option to include additional bitmask array. This will be added in the 3rd
        extension. Mask_options will be forced to bitmask. It will be applied in place
        if apply_mask==True, otherwise only the pipeline mask will be applied in place.
        If this is set, an error array is automatically applied in 2nd extension.
    fill_value: float, optional
        Value used to fill in for requested points outside of coverage or in a mask
        region. If not provided, then the default is 0.0.
    include_grid: bool
        Option to include xgrid, ygrid. This is used in lya_pyimfit.py. It is an array
        containing distance from center in arcsec matched to datagrid
    shot_h5: str
            optionally pass a specific <shot>.h5 fqfn

    Returns
    -------
    hdu: PrimaryHDU object
        the 2D summed data array and associated 2d header
        Units are '10^-17 erg cm-2 s-1'
        If include_error=True will include addiional hdu

    Examples
    --------

    For a specific detectid:
    >>> hdu = make_narrowband_image(detectid=2101046271)

    For a SkyCoords object. You must provide shotid and
    wavelength range

    >>> coords = SkyCoord(188.79312, 50.855747, unit='deg')
    >>> wave_obj = 4235.84 #in Angstrom
    >>> hdu = make_narrowband_image(coords=coords,
                                    shotid=20190524021,
                                    wave_range=[wave_obj-10, wave_obj+10])
    """
    global config, current_hdr, surveyh5
    global CONFIG_HDR2, CONFIG_HDR3, CONFIG_HDR4, CONFIG_HDR5, OPEN_DET_FILE, DET_HANDLE
    global DET_FILE

    if include_bitmask and not include_error:
        include_error = True
        print('Including bitmask and error arrays. Forcing include_error=True')

    if shot_h5 is not None:  # if shot_h5 is specified, use it instead of the survey
        if surveyh5 is not None:
            try:
                surveyh5.close()
            except:
                pass
    elif survey != current_hdr:
        config = HDRconfig(survey)
        current_hdr = survey

        try:
            if surveyh5 is not None:
                try:
                    surveyh5.close()
                except:
                    pass

            surveyh5 = tb.open_file(config.surveyh5, "r")
        except:
            pass
    else:
        # check to see if survey file is closed and open if so
        try:
            surveyh5.root
        except AttributeError:
            surveyh5 = tb.open_file(config.surveyh5, "r")

    if pixscale is None:
        pixscale = 0.25 * u.arcsec

    if imsize is None:
        imsize = 30.0 * u.arcsec

    if detectid is not None:
        # do not allow for detectid + any shotid/coords/wave_range
        if shotid is not None or wave_range is not None or coords is not None:
            print(
                "Can not override detectid paramters. Either enter a detectid or a shotid/coords/wave_range combination. Exiting"
            )
            return None

        if (detectid >= 2100000000) * (detectid < 2190000000):
            DET_FILE = CONFIG_HDR2.detecth5
        elif (detectid >= 2100000000) * (detectid < 3000000000):
            DET_FILE = CONFIG_HDR2.contsourceh5
        elif (detectid >= 3000000000) * (detectid < 3090000000):
            DET_FILE = CONFIG_HDR3.detecth5
        elif (detectid >= 3090000000) * (detectid < 3100000000):
            DET_FILE = CONFIG_HDR3.contsourceh5
        elif (detectid >= 4000000000) * (detectid < 4090000000):
            DET_FILE = CONFIG_HDR4.detecth5
        elif (detectid >= 4090000000) * (detectid < 4100000000):
            DET_FILE = CONFIG_HDR4.contsourceh5
        elif (detectid >= 5000000000) * (detectid < 5090000000):
            DET_FILE = CONFIG_HDR5.detecth5
        elif (detectid >= 5090000000) * (detectid < 5100000000):
            DET_FILE = CONFIG_HDR5.contsourceh5

        if OPEN_DET_FILE is None:
            DET_HANDLE = tb.open_file(DET_FILE, "r")
            OPEN_DET_FILE = DET_FILE
        else:
            if DET_FILE == OPEN_DET_FILE:
                pass
            else:
                DET_HANDLE.close()
                OPEN_DET_FILE = DET_FILE
                try:
                    DET_HANDLE = tb.open_file(DET_FILE, "r")
                except:
                    print("Could not open {}".format(DET_FILE))

        detectid_obj = detectid
        det_info = DET_HANDLE.root.Detections.read_where("detectid == detectid_obj")[0]

        shotid_obj = det_info["shotid"]
        wave_obj = det_info["wave"]
        linewidth = det_info["linewidth"]
        wave_range = [wave_obj - 2.0 * linewidth, wave_obj + 2.0 * linewidth]
        coords = SkyCoord(det_info["ra"], det_info["dec"], unit="deg")

    elif coords is not None:
        if shotid is not None:
            shotid_obj = shotid
        else:
            print("Provide a shotid")
        if wave_range is None:
            print(
                "Provide a wavelength range to collapse. \
            Example wave_range=[4500,4540]"
            )
    else:
        print("Provide a detectid or both a coords and shotid")

    if shot_h5 is not None:
        h5 = tables.open_file(shot_h5)
        fwhm = h5.root.Shot.read(field="fwhm_virus")[0]
        pa = h5.root.Shot.read(field="pa")[0]
        h5.close()
    else:
        fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]
        pa = surveyh5.root.Survey.read_where("shotid == shotid_obj")["pa"][0]

    if extract_class is None:
        E = Extract()
        E.load_shot(shotid_obj, fibers=False, survey=survey)
    else:
        E = extract_class
    # get spatial dims:
    ndim  = int(np.round((imsize / pixscale).to_value(u.dimensionless_unscaled)))
    center = (ndim + 1) / 2.0   # works for even/odd
    
    rad = imsize.to(u.arcsec).value  # convert to arcsec value, not quantity

    if include_bitmask:
        mask_options = "bitmask"
        apply_mask = True

        if include_error is False:
            include_error = True

    info_result = E.get_fiberinfo_for_coord(
        coords,
        radius=rad,
        ffsky=ffsky,
        ffsky_rescor=ffsky_rescor,
        fiber_flux_offset=fiber_flux_offset,
        add_mask=apply_mask,
        mask_options=mask_options,
        shot_h5=shot_h5
    )

    ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result

    # get ifu center:
    ifux_cen, ifuy_cen = E.convert_radec_to_ifux_ifuy(
        ifux, ifuy, ra, dec, coords.ra.deg, coords.dec.deg
    )

    #clip spectra edges
    w0, w1 = map(float, wave_range)
    w0 = max(w0, 3500.0); w1 = min(w1, 5500.0)
    wave_range = [w0, w1]

    if include_error:
        zarray = E.make_narrowband_image(
            ifux_cen,
            ifuy_cen,
            ifux,
            ifuy,
            data,
            mask,
            error=error,
            seeing_fac=fwhm,
            scale=pixscale.to(u.arcsec).value,
            boxsize=imsize.to(u.arcsec).value,
            wrange=wave_range,
            convolve_image=convolve_image,
            interp_kind=interp_kind,
            fill_value=fill_value,
            bitmask=include_bitmask,
        )

        imslice = zarray[0]
        imerror = zarray[1]

        if include_bitmask:
            imbitmask = zarray[2]
            xgrid = zarray[3]
            ygrid = zarray[4]
        else:
            xgrid = zarray[2]
            ygrid = zarray[3]
    else:
        zarray = E.make_narrowband_image(
            ifux_cen,
            ifuy_cen,
            ifux,
            ifuy,
            data,
            mask,
            seeing_fac=fwhm,
            scale=pixscale.to(u.arcsec).value,
            boxsize=imsize.to(u.arcsec).value,
            wrange=wave_range,
            convolve_image=convolve_image,
            interp_kind=interp_kind,
            fill_value=fill_value,
        )

        imslice = zarray[0]
        xgrid = zarray[1]
        ygrid = zarray[2]

    if subcont:
        wave_range_blue = [wave_range[0] - dcont - 10, wave_range[0] - 10]

        if wave_range_blue[0] <= 3500:
            wave_range_blue[0] = 3500
        if wave_range_blue[1] <= 3520:
            zarray_blue = None
        else:
            zarray_blue = E.make_narrowband_image(
                ifux_cen,
                ifuy_cen,
                ifux,
                ifuy,
                data,
                mask,
                seeing_fac=fwhm,
                scale=pixscale.to(u.arcsec).value,
                boxsize=imsize.to(u.arcsec).value,
                wrange=wave_range_blue,
                convolve_image=convolve_image,
                interp_kind=interp_kind,
                fill_value=fill_value,
            )

        wave_range_red = [wave_range[1] + 10, wave_range[1] + dcont + 10]

        if wave_range_red[1] >= 5500:
            wave_range_red[1] = 5500
        if wave_range_red[0] >= 5500:
            zarray_red = None
        else:
            zarray_red = E.make_narrowband_image(
                ifux_cen,
                ifuy_cen,
                ifux,
                ifuy,
                data,
                mask,
                seeing_fac=fwhm,
                scale=pixscale.to(u.arcsec).value,
                boxsize=imsize.to(u.arcsec).value,
                wrange=wave_range_red,
                convolve_image=convolve_image,
                interp_kind=interp_kind,
                fill_value=fill_value,
            )

        dwave = wave_range[1] - wave_range[0]

        if zarray_blue is None:
            im_cont = zarray_red[0] / dcont
        elif zarray_red is None:
            im_cont = zarray_blue[0] / dcont
        else:
            im_cont = (zarray_blue[0] + zarray_red[0]) / (2 * dcont)

        imslice = zarray[0] - dwave * im_cont

    w = wcs.WCS(naxis=2)
    imsize = imsize.to(u.arcsec).value
    w.wcs.crval = [coords.ra.deg, coords.dec.deg]
    w.wcs.crpix = [center, center]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cdelt = [-pixscale.to(u.deg).value, pixscale.to(u.deg).value]
    
    # get rotation:
    sys_rot = 1.55
    rot = 360.0 - (90.0 + pa + sys_rot)
    rrot = np.deg2rad(rot)

    w.wcs.crota = [ 0, rot]
    w.wcs.cunit = ["deg", "deg"]
    
    header = w.to_header()

    # add chosen variable info
    header["FILETIME"] = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    header["SURVEY"] = str(survey)
    header["MASK"] = str(apply_mask)
    header["CONVOLVE"] = str(convolve_image)
    header["INTERP"] = interp_kind
    header["SUBCONT"] = str(subcont)
    header["FFSKY"] = str(ffsky)
    header["FFSKYRC"] = str(ffsky_rescor)
        
    # Copy Shot table info
    if shot_h5 is not None:
        h5 = tables.open_file(shot_h5)
        shot_info_table = Table(h5.root.Shot.read())
        h5.close()
    else:
        shot_info_table = Table(surveyh5.root.Survey.read_where("shotid == shotid_obj"))
    for col in shot_info_table.colnames:
        if col.startswith("x") or col.startswith("y"):
            continue
        if "nstars" in col:
            continue
        if "darktime" in col:
            continue
        if "response" in col:
            header["RESPONSE"] = shot_info_table[col][0]
        elif "virus" in col:
            headername = col.replace("_virus", "")
            try:
                header[headername.upper()] = shot_info_table[col][0]
            except:
                header[headername.upper() + "1"] = shot_info_table[col][0][0]
                header[headername.upper() + "2"] = shot_info_table[col][0][1]
                header[headername.upper() + "3"] = shot_info_table[col][0][2]
        else:
            headername = str
            try:
                header[col.upper()] = shot_info_table[col][0]
            except:
                header[col.upper() + "1"] = shot_info_table[col][0][0]
                header[col.upper() + "2"] = shot_info_table[col][0][1]
                header[col.upper() + "3"] = shot_info_table[col][0][2]

    if extract_class is None:
        E.close()

    # create an empty Primary HDU for 1st extension
    hdu_primary = fits.PrimaryHDU()
    hdu_data = fits.ImageHDU(imslice.astype(np.float32), header=header, name="DATA")

    if include_grid:
        hdu_x = fits.ImageHDU(xgrid, header=header, name='XGRID')
        hdu_y = fits.ImageHDU(ygrid, header=header, name='YGRID')
        
    if include_error:
        
        hdu_error = fits.ImageHDU(imerror.astype(np.float32), header=header, name="ERROR")
        if include_bitmask:
            hdu_bitmask = fits.ImageHDU(imbitmask.astype(np.int16), header=header, name="BITMASK")
            if include_grid:
                return fits.HDUList([hdu_primary, hdu_data, hdu_error, hdu_bitmask, hdu_x, hdu_y])
            else:
                return fits.HDUList([hdu_primary, hdu_data, hdu_error, hdu_bitmask])
        else:
            if include_grid:
                return fits.HDUList([hdu_primary, hdu_data, hdu_error, hdu_x, hdu_y])
            else:
                return fits.HDUList([hdu_primary, hdu_data, hdu_error])
    else:
        if include_grid:
            return fits.HDUList([hdu_primary, hdu_data, hdu_x, hdu_y])
        else:
            return fits.HDUList([hdu_primary, hdu_data])
    

def make_data_cube(
    detectid=None,
    coords=None,
    shotid=None,
    pixscale=None,
    imsize=None,
    wave_range=None,
    dwave=2.0,
    convolve_image=False,
    ffsky=False,
    ffsky_rescor=False,
    survey=LATEST_HDR_NAME,
    fiber_flux_offset=None,
    interp_kind="linear",
    apply_mask=False,
    include_error=False,
    include_bitmask=False,
    mask_options=None,
    fill_value=0.0,
    extract_class=None,
):
    """
    Function to make a datacube from either a detectid or from a
    coordinate/shotid combination.

    Paramaters
    ----------
    detectid: int
        detectid from the continuum or lines catalog. Default is
        None. Provide a coords/shotid combo if this isn't given
    coords: SkyCoords object
        coordinates to define the centre of the data cube
    pixscale: astropy angle quantity
        plate scale
    imsize: astropy angle quantity
        spatial length of cube (equal dims is only option)
    wave_range: list
        start and stop value for the wavelength range in Angstrom
    dwave
        step in wavelength range in Angstrom
    convolve_image: bool
         option to convolve image with shotid seeing
    ffsky: bool
        option to use full frame calibrated fibers. Default False
    ffsky_rescor: bool
        option to use full frame calibrated fibers. Default is False
    fiber_flux_offset: 1036 array
        array of values in units of 10**-17 ergs/s/cm2/AA to add
        to each fiber spectrum used in the extraction. Defaults
        to None
    interp_kind: str
        Kind of interpolation to pixelated grid from fiber intensity.
        Options are 'linear', 'cubic', 'nearest'. Default is linear.
    add_mask: bool                                                                                        Apply HETDEX fiber mask model. This will mask all fibers contributing
        to the spectral extraction before summation. Masked in place according
        to fill_value
    fill_value: float, optional
        Value used to fill in for requested points outside of coverage or in a mask
        region. If not provided, then the default is 0.0.
    apply_mask: bool
        Apply HETDEX fiber mask model. This will mask all fibers contributing
        to the spectral extraction before summation. Masked in place as NaN
    mask_options
        string or array of strings as options to select to mask. Default None
        will select all flags. Set this to 'BITMASK' to return the full bitmask array.
        Options are 'MAIN', 'FTF', 'CHI2FIB', 'BADPIX', 'BADAMP', 'LARGEGAL', 'METEOR',
        'BADSHOT', 'THROUGHPUT', 'BADFIB', 'SAT'
        If BITMASK appears as any element in the list, it overrides all others
        and returns the full bitmask array.
    include_bitmask
        option to include additional bitmask array. This will be added in the 3rd
        extension. Mask_options will be forced to bitmask. It will be applied in place
        if apply_mask==True, otherwise only the pipeline mask will be applied in place.
        If this is set, an error array is automatically applied in 2nd extension.

    Returns
    -------
    hdu: PrimaryHDU object
        the data cube 3D array and associated 3d header
        Units are '10^-17 erg cm-2 s-1 per spaxel'

    Examples
    --------

    Can either pass in a detectid:

    >>> detectid_obj=2101602788
    >>> hdu = make_data_cube( detectid=detectid_obj)
    >>> hdu.writeto( str(detectid_obj) + '.fits', overwrite=True)

    or can put in an SkyCoord object:

    >>> star_coords = SkyCoord(9.625181, -0.043587, unit='deg')
    >>> hdu = make_data_cube( coords=star_coords[0], shotid=20171016108, dwave=2.0)
    >>> hdu.writeto( 'star.fits', overwrite=True)

    """
    global config, current_hdr, surveyh5
    global CONFIG_HDR2, CONFIG_HDR3, OPEN_DET_FILE, DET_HANDLE
    global DET_FILE

    if pixscale is None:
        pixscale = 0.25 * u.arcsec

    if imsize is None:
        imsize = 30.0 * u.arcsec

    if wave_range is None:
        wave_range = [3470.0, 5540.0]
        
    if survey != current_hdr:
        config = HDRconfig(survey)
        current_hdr = survey

        try:
            if surveyh5 is not None:
                try:
                    surveyh5.close()
                except:
                    pass
            surveyh5 = tb.open_file(config.surveyh5, "r")
        except:
            pass
    else:
        # check to see if survey file is closed and open if so
        try:
            surveyh5.root
        except AttributeError:
            surveyh5 = tb.open_file(config.surveyh5, "r")

    if detectid is not None:
        if (detectid >= 2100000000) * (detectid < 2190000000):
            DET_FILE = CONFIG_HDR2.detecth5
        elif (detectid >= 2100000000) * (detectid < 2190000000):
            DET_FILE = CONFIG_HDR2.contsourceh5
        elif (detectid >= 3000000000) * (detectid < 3090000000):
            DET_FILE = CONFIG_HDR3.detecth5
        elif (detectid >= 3090000000) * (detectid < 3100000000):
            DET_FILE = CONFIG_HDR3.contsourceh5
        elif (detectid >= 4000000000) * (detectid < 4090000000):
            DET_FILE = CONFIG_HDR4.detecth5
        elif (detectid >= 4090000000) * (detectid < 4100000000):
            DET_FILE = CONFIG_HDR4.contsourceh5
        elif (detectid >= 5000000000) * (detectid < 5090000000):
            DET_FILE = CONFIG_HDR5.detecth5
        elif (detectid >= 5090000000) * (detectid < 5100000000):
            DET_FILE = CONFIG_HDR5.contsourceh5

        if OPEN_DET_FILE is None:
            OPEN_DET_FILE = DET_FILE
            DET_HANDLE = tb.open_file(DET_FILE, "r")

        else:
            if DET_FILE == OPEN_DET_FILE:
                pass
            else:
                DET_HANDLE.close()
                OPEN_DET_FILE = DET_FILE
                try:
                    DET_HANDLE = tb.open_file(DET_FILE, "r")
                except Exception:
                    print("Could not open {}".format(DET_FILE))

        detectid_obj = detectid
        det_info = DET_HANDLE.root.Detections.read_where("detectid == detectid_obj")[0]

        shotid = det_info["shotid"]
        coords = SkyCoord(det_info["ra"], det_info["dec"], unit="deg")

    if coords is None or shotid is None:
        print("Provide a detectid or both a coords and shotid")

    if extract_class is None:
        E = Extract()
        E.load_shot(shotid, fibers=False, survey=survey)
    else:
        E = extract_class

    # get spatial dims:
    #ndim = int(imsize / pixscale)
    #center = int(ndim / 2)

    # updated 2025-08-18 by EMC
    ndim = int(np.round((imsize / pixscale).to_value(u.dimensionless_unscaled)))
    center = (ndim + 1) / 2.0   # may be half-integer if ndim is even

    # get wave dims: (edited 2025-08-14 by EMC)
    #nwave = int((wave_range[1] - wave_range[0]) / dwave + 1)

    # Number of spectral bins, not edges
    nwave = int(np.round((wave_range[1] - wave_range[0]) / dwave)) + 1
    wave_centers = np.linspace(wave_range[0], wave_range[1], nwave)

    # --- WCS ---
    # Bin centers (added 2025-08-14 by EMC)
    
    w = wcs.WCS(naxis=3)
    w.wcs.crval = [coords.ra.deg, coords.dec.deg, float( wave_range[0])] 
    w.wcs.crpix = [center, center, 1.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
    w.wcs.cdelt = [-pixscale.to(u.deg).value, pixscale.to(u.deg).value, float( dwave)]
    #added 2025-08-14 by EMC
    w.wcs.cunit = ["deg", "deg", "Angstrom"]
    
    rad = imsize.to(u.arcsec).value

    if include_bitmask:
        mask_options = "bitmask"
        apply_mask = True

        if include_error is False:
            include_error = True

    info_result = E.get_fiberinfo_for_coord(
        coords,
        radius=rad,
        ffsky=ffsky,
        ffsky_rescor=ffsky_rescor,
        fiber_flux_offset=fiber_flux_offset,
        add_mask=apply_mask,
        mask_options=mask_options,
    )

    ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result

    # get ifu center:
    ifux_cen, ifuy_cen = E.convert_radec_to_ifux_ifuy(
        ifux, ifuy, ra, dec, coords.ra.deg, coords.dec.deg
    )

    # get FWHM and PA
    shotid_obj = shotid

    pa = surveyh5.root.Survey.read_where("shotid == shotid_obj")["pa"][0]

    if convolve_image:
        fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]
    else:
        fwhm = 1.8  # just a dummy variable as convolve_image=False

    # add in rotation
    sys_rot = 1.55
    rot = 360.0 - (90.0 + pa + sys_rot)

    w.wcs.crota = [0, rot, 0]

    # positive magnitudes for scales
    sx = pixscale.to(u.deg).value
    sy = pixscale.to(u.deg).value
    theta = np.deg2rad(rot)

    lam0 = float(wave_range[0])   # 3470                                                           
    dlam = float(dwave)           # 2.0           
    
    im_cube = np.zeros((nwave, ndim, ndim), dtype=np.float32)

    if include_error:
        im_errorcube = np.zeros((nwave, ndim, ndim))
        if include_bitmask:
            im_bitmaskcube = np.zeros((nwave, ndim, ndim), dtype=np.uint16)

    for i, wave_i in enumerate(wave_centers):
        
        half = 0.5 * dwave
        wrange=[wave_i - half, wave_i + half]

        try:
            zarray = E.make_narrowband_image(
                ifux_cen,
                ifuy_cen,
                ifux,
                ifuy,
                data,
                mask,
                scale=pixscale.to(u.arcsec).value,
                wrange=wrange,
                seeing_fac=fwhm,
                convolve_image=convolve_image,
                boxsize=imsize.to(u.arcsec).value,
                interp_kind=interp_kind,
                fill_value=fill_value,
                error=error,
                bitmask=include_bitmask,
            )
            im_cube[i] = zarray[0]

            if include_error:
                im_errorcube[i] = zarray[1]
                if include_bitmask:
                    im_bitmaskcube[i] = zarray[2]
        except ValueError:
            pass

    # CREATE HEADER
    header = w.to_header()

    # add chosen variable info
    header["FILETIME"] = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    header["SURVEY"] = str(survey)
    header["MASK"] = str(apply_mask)
    header["CONVOLVE"] = str(convolve_image)
    header["INTERP"] = interp_kind
    header["FFSKY"] = str(ffsky)
    header["FFSKYRC"] = str(ffsky_rescor)

    # Copy Shot table info
    shot_info_table = Table(surveyh5.root.Survey.read_where("shotid == shotid_obj"))
    for col in shot_info_table.colnames:
        if col.startswith("x") or col.startswith("y"):
            continue
        if "nstars" in col:
            continue
        if "darktime" in col:
            continue
        if "response" in col:
            header["RESPONSE"] = shot_info_table[col][0]
        elif "virus" in col:
            headername = col.replace("_virus", "")
            try:
                header[headername.upper()] = shot_info_table[col][0]
            except:
                header[headername.upper() + "1"] = shot_info_table[col][0][0]
                header[headername.upper() + "2"] = shot_info_table[col][0][1]
                header[headername.upper() + "3"] = shot_info_table[col][0][2]
        else:
            headername = str
            try:
                header[col.upper()] = shot_info_table[col][0]
            except:
                header[col.upper() + "1"] = shot_info_table[col][0][0]
                header[col.upper() + "2"] = shot_info_table[col][0][1]
                header[col.upper() + "3"] = shot_info_table[col][0][2]

    # enforce header cards for spectral axis
    header["CTYPE3"] = "WAVE"
    header["CUNIT3"] = "Angstrom"
    header["CRVAL3"] = lam0
    header["CDELT3"] = dlam
    header["CRPIX3"] = 1.0

    # clean out any CD/PC matrix elements for axis 3 that could override
    for key in list(header.keys()):
        if key.startswith("PC3_") or key.startswith("CD3_"):
            del header[key]

    # close extract_class if just accessing for one call
    if extract_class is None:
        try:
            E.close()
        except:
            pass

    # create an empty Primary HDU for 1st extension
    hdu_primary = fits.PrimaryHDU()
    hdu_data = fits.ImageHDU(im_cube.astype(np.float32), header=header, name="DATA")

    if include_error:
        hdu_error = fits.ImageHDU(im_errorcube.astype(np.float32), header=header, name="ERROR")

        if include_bitmask:
            hdu_bitmask = fits.ImageHDU(im_bitmaskcube.astype(np.uint16), header=header, name="BITMASK")
            return fits.HDUList([hdu_primary, hdu_data, hdu_error, hdu_bitmask])
        else:
            return fits.HDUList([hdu_primary, hdu_data, hdu_error])
    else:
        return fits.HDUList([hdu_primary, hdu_data])

