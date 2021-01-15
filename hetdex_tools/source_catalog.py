import numpy as np
import time
import glob
import argparse as ap

from astropy.table import Table, vstack, join, Column, unique
from astropy.coordinates import SkyCoord
import astropy.units as u

import matplotlib

matplotlib.use("agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from regions import LineSkyRegion, PixCoord, LinePixelRegion

from hetdex_api.config import HDRconfig
from hetdex_api.survey import Survey
from hetdex_api.detections import Detections
import hetdex_tools.fof_kdtree as fof

from hetdex_api.elixer_widget_cls import ElixerWidget

from elixer import catalogs

catlib = catalogs.CatalogLibrary()
config = HDRconfig()


def create_source_catalog(version="2.1.2", make_continuum=True, save=True, dsky=5.0):

    global config

    detects = Detections(curated_version=version)
    detects_line_table = detects.return_astropy_table()
    detects_line_table.add_column(Column("line", name="det_type", dtype=str))
    detects_cont = Detections(catalog_type="continuum")

    sel1 = detects_cont.remove_bad_amps()
    sel2 = detects_cont.remove_meteors()
    sel3 = detects_cont.remove_shots()

    detects_cont_table = detects_cont[sel1 * sel2 * sel3].return_astropy_table()
    detects_cont_table.add_column(Column("cont", name="det_type", dtype=str))

    if make_continuum:
        detects_cont_table.write("continuum_" + version + ".fits", overwrite=True)

    dets_all = Detections().refine()
    dets_all_table = dets_all.return_astropy_table()
    agn_tab = Table.read(config.agncat, format="ascii", include_names=["detectid"])

    # add in continuum sources to match to Chenxu's combined catalog
    detects_cont_table_orig = detects_cont[sel1 * sel2 * sel3].return_astropy_table()
    dets_all_table = vstack([dets_all_table, detects_cont_table_orig])

    detects_broad_table = join(
        agn_tab, dets_all_table, join_type="inner", keys=["detectid"]
    )
    detects_broad_table.add_column( Column('agn', name='det_type', dtype=str))
    dets_all.close()
    del dets_all_table

    detect_table = vstack([detects_line_table,
                           detects_broad_table,
                           detects_cont_table])

    del detects_cont_table_orig, detects_cont_table, detects_broad_table
    
    kdtree, r = fof.mktree(
        detect_table["ra"],
        detect_table["dec"],
        np.zeros_like(detect_table["ra"]),
        dsky=dsky,
    )
    t0 = time.time()
    print("starting fof ...")
    friend_lst = fof.frinds_of_friends(kdtree, r, Nmin=1)
    t1 = time.time()

    print("FOF analysis complete in {:3.2f} minutes \n".format((t1 - t0) / 60))

    # get fluxes to derive flux-weighted distribution of group

    detect_table["gmag"][np.isnan(detect_table["gmag"])] = 27
    gmag = detect_table["gmag"] * u.AB
    flux = gmag.to(u.Jy).value

    friend_table = fof.process_group_list(
        friend_lst,
        detect_table["detectid"],
        detect_table["ra"],
        detect_table["dec"],
        0.0 * detect_table["wave"],
        flux,
    )

    print("Generating combined table \n")
    memberlist = []
    friendlist = []
    for row in friend_table:
        friendid = row["id"]
        members = np.array(row["members"])
        friendlist.extend(friendid * np.ones_like(members))
        memberlist.extend(members)
    friend_table.remove_column("members")

    detfriend_tab = Table()
    detfriend_tab.add_column(Column(np.array(friendlist), name="id"))
    detfriend_tab.add_column(Column(memberlist, name="detectid"))

    detfriend_all = join(detfriend_tab, friend_table, keys="id")

    del detfriend_tab
    
    expand_table = join(detfriend_all, detect_table, keys="detectid")

    del detfriend_all, detect_table, friend_table

    if version == '2.1.2':#will want to improve this logic
        starting_id = int(2.12 * 10 ** 12)

    expand_table.add_column(
        Column(expand_table["id"] + starting_id, name="source_id"), index=0
    )

    expand_table.remove_column("id")

    gaia_stars = Table.read(config.gaiacat)

    gaia_coords = SkyCoord(ra=gaia_stars["ra"] * u.deg, dec=gaia_stars["dec"] * u.deg)
    src_coords = SkyCoord(
        ra=expand_table["ra"] * u.deg, dec=expand_table["dec"] * u.deg
    )

    idx, d2d, d3d = src_coords.match_to_catalog_sky(gaia_coords)

    sel = d2d < 5.0 * u.arcsec

    gaia_match_name = np.zeros_like(expand_table["source_id"], dtype=int)
    gaia_match_name[sel] = gaia_stars["source_id"][idx][sel]

    expand_table["gaia_match_id"] = gaia_match_name
    expand_table.rename_column("size", "n_members")
    expand_table.rename_column("icx", "ra_mean")
    expand_table.rename_column("icy", "dec_mean")

    expand_table.sort("source_id")

#    source_table = add_z_guess(expand_table)
#
#    del expand_table, gaia_stars, c

    if save:
        expand_table.write("source_catalog_" + version + ".fits", overwrite=True)

    return expand_table


def z_guess_3727(group, cont=False):
    sel_good_lines = None
    if np.any((group["sn"] > 15) * (group["wave"] > 3727)):
        sel_good_lines = (group["sn"] > 15) * (group["wave"] > 3727)
    elif np.any((group["sn"] > 10) * (group["wave"] > 3727)):
        sel_good_lines = (group["sn"] > 8) * (group["wave"] > 3727)
    elif np.any((group["sn"] > 8) * (group["wave"] > 3727)):
        sel_good_lines = (group["sn"] > 8) * (group["wave"] > 3727)
    elif np.any((group["sn"] > 6) * (group["wave"] > 3727)):
        sel_good_lines = (group["sn"] > 6) * (group["wave"] > 3727)
    elif np.any((group["sn"] > 5.5) * (group["wave"] > 3727)):
        if cont:
            pass
        else:
            sel_good_lines = (group["sn"] > 5.5) * (group["wave"] > 3727)
    else:
        if cont:
            pass
        else:
            sel_good_lines = group["wave"] > 3727

    if sel_good_lines is not None:
        try:
            wave_guess = np.min(group["wave"][sel_good_lines])
            z_guess = wave_guess / 3727.0 - 1
        except:
            z_guess = -2.0
        else:
            z_guess = 0.0

    return z_guess


def guess_source_wavelength(source_id, source_table):

    sel_group = source_table["source_id"] == source_id
    group = source_table[sel_group]
    z_guess = -1.0

    # for now assign a plae_classification for -1 values
    if np.any(group["det_type"] == "agn"):
        # get proper z's from Chenxu's catalog
        agn_tab = Table.read(config.agncat, format="ascii")
        agn_dets = group["detectid"][group["det_type"] == "agn"]
        sel_det = chenxu_tab["detectid"] == agn_dets[0]
        z_guess = chenxu_tab["z"][sel_det][0]

    elif np.any(group["gaia_match_id"] > 0):
        z_guess = z_guess_3727(group, cont=True)

    elif np.any(group["det_type"] == "cont"):
        z_guess = z_guess_3727(group, cont=True)

    elif np.any(group["plae_classification"]) < 0.3:
        z_guess = z_guess_3727(group)

    elif np.nanmedian(group["plae_classification"] < 0.5):
        z_guess = z_guess_3727(group)

    elif np.nanmedian(group["plae_classification"] > 0.5):
        if np.any(group["sn"] > 15):
            sel_good_lines = group["sn"] > 15
            wave_guess = np.min(group["wave"][sel_good_lines])
        elif np.any(group["sn"] > 10):
            sel_good_lines = group["sn"] > 10
            wave_guess = np.min(group["wave"][sel_good_lines])
        elif np.any(group["sn"] > 8):
            sel_good_lines = group["sn"] > 8
            wave_guess = np.min(group["wave"][sel_good_lines])
        elif np.any(group["sn"] > 6.5):
            sel_good_lines = group["sn"] > 6.5
            wave_guess = np.min(group["wave"][sel_good_lines])
        elif np.any(group["sn"] > 5.5):
            sel_good_lines = group["sn"] > 5.5
            wave_guess = np.min(group["wave"][sel_good_lines])
        else:
            argmaxsn = np.argmax(group["sn"])
            wave_guess = group["wave"][argmaxsn]
            z_guess = wave_guess / 1216.0 - 1
    else:
        argmaxsn = np.argmax(group["sn"])
        wave_guess = group["wave"][argmaxsn]
        z_guess = wave_guess / 1216.0 - 1

    return z_guess


def add_z_guess(source_table):

    try:
        # remove z_guess column if it exists
        print("Removing existing z_guess column")
        source_table.remove_column("z_guess")
    except Exception:
        pass

    print("Assigning a best guess redshift")

    from multiprocessing import Pool

    src_list = np.unique(source_table["source_id"])

    t0 = time.time()
    p = Pool()
    src_z = p.map(guess_source_wavelength, src_list)
    t1 = time.time()
    p.close()

    print("Finished in {:3.2f} minutes".format((t1 - t0) / 60))

    z_guess = np.array(src_z)
    z_table = Table([src_list, z_guess], names=["source_id", "z_guess"])

    all_table = join(source_table, z_table, join_type="left")

    return all_table


def plot_source_group(source_id=None,
                      source_table=None,
                      k=3.5, vmin=3, vmax=99,
                      label=True,
                      save=False):
    """
    Plot a unique source group from the HETDEX
    unique source catalog
    
    Parameters
    ----------
    source_id: int
    
    """

    if source_table is None:
        print("Please provide source catalog (an astropy table)")
    else:
        sel = source_table["source_id"] == source_id
        group = source_table[sel]

    if source_id is None:
        print("Please provide source_id (an integer)")

    ellipse = False

    # get ellipse parameters if more than 1 source
    if np.size(group) > 1:
        ellipse = True

    if ellipse:
        a, b, pa, a2, b2, pa2 = (
            group["a"][0],
            group["b"][0],
            group["pa"][0],
            group["a2"][0],
            group["b2"][0],
            group["pa2"][0],
        )

        cosd = np.cos(np.deg2rad(group["dec_mean"][0]))
        imsize = (
            np.max(
                [
                    (np.max(group["ra"]) - np.min(group["ra"])) * cosd,
                    (np.max(group["dec"]) - np.min(group["dec"])),
                ]
            )
            * 1.7
        )

    else:
        imsize = 10.0 / 3600.0

    coords = SkyCoord(
        ra=group["ra_mean"][0] * u.deg, dec=group["dec_mean"][0] * u.deg
    )

    cutout = catlib.get_cutouts(
        position=coords,
        side=imsize,
        filter=["f606W", "r", "g"],
        first=True,
        allow_bad_image=False,
        allow_web=True,
    )[0]
    
    im = cutout["cutout"].data
    wcs = cutout["cutout"].wcs
    plt.figure(figsize=(8, 8))

    plt.subplot(projection=wcs)
    ax = plt.gca()
    ax.coords[0].set_major_formatter("d.dddd")
    # ax.coords[0].set_ticks(spacing=1. * u.arcsec)
    ax.coords[1].set_major_formatter("d.dddd")
    # ax.coords[0].set_ticks(spacing=1. * u.arcsec)

    impix = im.shape[0]
    pixscale = imsize / impix  # in degrees/pix
    m = np.percentile(im, (vmin, vmax))
    plt.imshow(im, vmin=m[0], vmax=m[1], origin="lower", cmap="gray_r")
    plt.text(
        0.95,
        0.05,
        cutout["instrument"] + cutout["filter"],
        transform=ax.transAxes,
        fontsize=20,
        color="red",
        horizontalalignment="right",
    )


    # plot the group members
    sel_line = group['det_type'] == 'line'
    if np.sum(sel_line) >= 1:
        plt.scatter(
            group['ra'][sel_line],
            group['dec'][sel_line],
            transform=ax.get_transform("world"),
            marker="o",
            color="orange",
            linewidth=4,
            s=group['sn'][sel_line],
            zorder=100,
            label='line emission'
        )
        
        
    sel_cont = group['det_type'] == 'cont'
    if np.sum(sel_cont)>= 1:
        plt.scatter(
            group['ra'][sel_cont],
            group['dec'][sel_cont],
            transform=ax.get_transform("world"),
            marker="o",
            color="green",
            linewidth=4,
            s=10,
            label='continuum'
        )
        
    sel_agn = group['det_type'] == 'agn'
    if np.sum(sel_agn)>= 1:
        plt.scatter(
            group['ra'][sel_agn],
            group['dec'][sel_agn],
            transform=ax.get_transform("world"),
            marker="o",
            color="red",
            linewidth=4,
            s=10,
            label='agn')
        
    # plot and elliptical kron-like aperture representing the group. Ellipse
    # in world coords does not work, so plot in pixelcoordinates...
    # East is to the right in these plots, so pa needs transform

    if ellipse:
        ellipse = Ellipse(
            xy=(impix // 2, impix // 2),
            width=k * a2 / pixscale,
            height=k * b2 / pixscale,
            angle=180 - pa2,
            edgecolor="r",
            fc="None",
            lw=1,
        )
        ax.add_patch(ellipse)

        ellipse = Ellipse(
            xy=(impix // 2, impix // 2),
            width=a / pixscale,
            height=b / pixscale,
            angle=180 - pa,
            edgecolor="b",
            fc="None",
            lw=1,
        )
        ax.add_patch(ellipse)

    # add 5 arcsec scale bar
    x1 = 0.05 * impix
    y1 = 0.05 * impix
    y2 = y1 + (5.0 / 3600) / pixscale
    start = PixCoord(x=x1, y=y1)
    end = PixCoord(x=x1, y=y2)
    reg = LinePixelRegion(start=start, end=end)
    plt.text(x1, 0.025 * impix, "5 arcsec", color="blue")
    patch = reg.as_artist(facecolor="none", edgecolor="blue", lw=4)
    ax.add_patch(patch)

    if label:
        # plot detecid labels
        for row in group:
            plt.text(
                row["ra"],
                row["dec"],
                str(row["detectid"]),
                transform=ax.get_transform("world"),
                fontsize=9,
                color="red",
            )
#    z_guess = guess_source_wavelength(source_id, source_table)
    z_guess = group['z_guess'][0]

    plt.title(
        "source_id:%d n:%d ra:%6.3f dec:%6.3f z:%4.3f"
        % (
            source_id,
            group["n_members"][0],
            group["ra_mean"][0],
            group["dec_mean"][0],
            z_guess,
        )
    )

    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.legend()
    
    if save:
        plt.savefig("figures/source-%03d.png" % source_id, format="png")
        plt.close()
    else:
        plt.show()


def get_parser():
    """ function that returns a parser from argparse """

    parser = ap.ArgumentParser(
        description="""Extracts 1D spectrum at specified RA/DEC""", add_help=True
    )

    parser.add_argument(
        "-v",
        "--version",
        help="""source catalog version you want to create""",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-dsky", "--dsky", help="""Spatial linking length in arcsec""", default=5.0,
    )

    return parser


def main(argv=None):
    """ Main Function """

    parser = get_parser()
    args = parser.parse_args(argv)
#    args.log = setup_logging()# set up later

    create_source_catalog(version=args.version, dsky=args.dsky)

if __name__ == "__main__":
    main()
