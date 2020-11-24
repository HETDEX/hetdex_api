import numpy as np
import tables as tb

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from astropy.io import fits

from hetdex_api.config import HDRconfig
from hetdex_api.extract import Extract

LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME

config = HDRconfig()
surveyh5 = tb.open_file(config.surveyh5, "r")
detecth5 = tb.open_file(config.detecth5, "r")


def make_narrowband_image(
    detectid=None,
    coords=None,
    shotid=None,
    pixscale=0.25 * u.arcsec,
    imsize=30.0 * u.arcsec,
    wave_range=None,
    convolve_image=True,
    ffsky=True,
    subcont=False,
    dcont=50.,
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
    subcont: bool
        option to subtract continuum. Default is False. This
        will measure the continuum 50AA below and above the
        input wave_range
    dcont
        width in angstrom to measure the continuum. Default is to
        measure 50 AA wide regions on either side of the line
         
    
    Returns
    -------
    hdu: PrimaryHDU object
        the 2D summed data array and associated 2d header
        Units are '10^-17 erg cm-2 s-1'
    
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
    global config, detecth5, surveyh5

    if detectid is not None:
        
        detectid_obj = detectid
        det_info = detecth5.root.Detections.read_where("detectid == detectid_obj")[0]
        shotid_obj = det_info["shotid"]
        wave_obj = det_info["wave"]
        linewidth = det_info["linewidth"]
        wave_range = [wave_obj - 2.0 * linewidth,
                      wave_obj + 2.0 * linewidth]
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

    fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]

    E = Extract()
    E.load_shot(shotid_obj)

    # get spatial dims:
    ndim = int(imsize / pixscale)
    center = int(ndim / 2)

    rad = imsize
    info_result = E.get_fiberinfo_for_coord(coords, radius=rad, ffsky=ffsky)
    ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result

    # get ifu center:
    ifux_cen, ifuy_cen = E.convert_radec_to_ifux_ifuy(
        ifux, ifuy, ra, dec, coords.ra.deg, coords.dec.deg
    )

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
    )

    imslice = zarray[0]

    if subcont:
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
            wrange=[wave_range[0]-dcont, wave_range[0]],
            convolve_image=convolve_image,
        )

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
            wrange=[wave_range[1], wave_range[1]+dcont],
            convolve_image=convolve_image,
        )
        
        dwave = wave_range[1]-wave_range[0]
        im_cont = (zarray_blue[0] + zarray_red[0])/(2*dcont)

        imslice = zarray[0] - dwave*im_cont


    w = wcs.WCS(naxis=2)
    imsize = imsize.to(u.arcsec).value
    w.wcs.crval = [coords.ra.deg, coords.dec.deg]
    w.wcs.crpix = [center, center]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cdelt = [-pixscale.to(u.deg).value, pixscale.to(u.deg).value]

    hdu = fits.PrimaryHDU(imslice, header=w.to_header())

    return hdu


def make_data_cube(
    detectid=None,
    coords=None,
    shotid=None,
    pixscale=0.25 * u.arcsec,
    imsize=30.0 * u.arcsec,
    wave_range=[3470, 5540],
    dwave=2.0,
    dcont=50.0,
    convolve_image=True,
    ffsky=True,
    subcont=False,
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
        option to use full frame calibrated fibers. Default is
        True.
    subcont: bool
        option to subtract continuum. Default is False. This
        will measure the continuum 50AA below and above the
        input wave_range
    dcont
        width in angstrom to measure the continuum. Default is to
        measure 50 AA wide regions on either side of the line  

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
    global config, detecth5, surveyh5

    if detectid is not None:
        detectid_obj = detectid
        det_info = detecth5.root.Detections.read_where("detectid == detectid_obj")[0]
        shotid = det_info["shotid"]
        coords = SkyCoord(det_info["ra"], det_info["dec"], unit="deg")

        if coords is None or shotid is None:
            print("Provide a detectid or both a coords and shotid")

    E = Extract()
    E.load_shot(shotid)

    # get spatial dims:
    ndim = int(imsize / pixscale)
    center = int(ndim / 2)

    # get wave dims:
    nwave = int((wave_range[1] - wave_range[0]) / dwave + 1)

    w = wcs.WCS(naxis=3)
    w.wcs.crval = [coords.ra.deg, coords.dec.deg, wave_range[0]]
    w.wcs.crpix = [center, center, 1]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
    w.wcs.cdelt = [-pixscale.to(u.deg).value, pixscale.to(u.deg).value, dwave]
    
    rad = imsize
    info_result = E.get_fiberinfo_for_coord(coords, radius=rad, ffsky=False)
    ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result

    # get ifu center:
    ifux_cen, ifuy_cen = E.convert_radec_to_ifux_ifuy(
        ifux, ifuy, ra, dec, coords.ra.deg, coords.dec.deg
    )

    if convolve_image:
        surveyh5 = tb.open_file(config.surveyh5, "r")
        shotid_obj = shotid
        fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]
        surveyh5.close()
    else:
        fwhm = 1.8  # just a dummy variable as convolve_image=False

    im_cube = np.zeros((nwave, ndim, ndim))

    wave_i = wave_range[0]
    i = 0

    while wave_i <= wave_range[1]:
        try:
            im_src = E.make_narrowband_image(
                ifux_cen,
                ifuy_cen,
                ifux,
                ifuy,
                data,
                mask,
                scale=pixscale.to(u.arcsec).value,
                wrange=[wave_i, wave_i + dwave],
                nchunks=1,
                seeing_fac=fwhm,
                convolve_image=convolve_image,
                boxsize=imsize.to(u.arcsec).value,
            )

            im_slice = im_src[0]

            if subcont:
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
                    nchunks=2,
                    wrange=[wave_i-dcont, wave_i],
                    convolve_image=convolve_image,
                )
                zarray_red = E.make_narrowband_image(
                    ifux_cen,
                    ifuy_cen,
                    ifux,
                    ifuy,
                    data,
                    mask,
                    seeing_fac=fwhm,
                    nchunks=2,
                    scale=pixscale.to(u.arcsec).value,
                    boxsize=imsize.to(u.arcsec).value,
                    wrange=[wave_i + dwave, wave_i + dwave + dcont],
                    convolve_image=convolve_image,
                )

                im_cont = (zarray_blue[0] + zarray_red[0])/(2*dcont)
                im_slice = im_src[0] - dwave*im_cont
    
            im_cube[i, :, :] = im_slice

        except Exception:
            im_cube[i, :, :] = np.zeros((ndim, ndim))
        wave_i += dwave
        i += 1

    hdu = fits.PrimaryHDU(im_cube, header=w.to_header())

    E.close()

    return hdu
