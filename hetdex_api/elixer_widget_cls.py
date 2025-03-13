from __future__ import print_function

"""
ELiXer Widget to Classify Detections

Based on Dr. Erin Mentuch Cooper's original elixer_widgets.py


"""

ALLOW_MANUAL_CATALOG_UPDATE = False
ALLOW_CLASSIFICATION_BUTTONS = True

import numpy as np
import matplotlib.pyplot as plt
import os
import os.path as op
import traceback

import astropy.units as u
from astropy.io import ascii
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import hetdex_api.sqlite_utils as sql

import ipywidgets as widgets
from ipywidgets import interact, Layout  # Style #, interactive

from IPython.display import display, Image, Javascript, HTML
from PIL import Image as PILImage
import io

# from IPython.display import clear_output
import tables

from hetdex_tools.get_spec import get_spectra
from hetdex_api.detections import Detections

if ALLOW_MANUAL_CATALOG_UPDATE:
    try:
        from h5tools import source_catalog_manual_updates
    except:
        ALLOW_MANUAL_CATALOG_UPDATE = False
        print("Cannot import necessary packages for source catalog manual updating options.")
        print("This is disabled until the issue is resolved.")
        print("For Jupyter-hub users ...")
        print("you may need to run:   !pip install filelock")
        print("then restart your kernel. This will hold until you restart the server.")

from elixer import spectrum_utilities as SU

try:
    from hetdex_api.config import HDRconfig
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
    CONFIG_HDR2 = HDRconfig('hdr2.1')
    CONFIG_HDR3 = HDRconfig('hdr3')
    CONFIG_HDR4 = HDRconfig('hdr4')
    CONFIG_HDR5 = HDRconfig('hdr5')
    
except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr5"
    CONFIG_HDR2 = None
    CONFIG_HDR3 = None

HDR_NAME_DICT = {10: "hdr1", 20: "hdr2", 21: "hdr2.1", 30: "hdr3"}


try:  # using HDRconfig
    HETDEX_API_CONFIG = HDRconfig(survey=LATEST_HDR_NAME)
    HDR_BASEPATH = HETDEX_API_CONFIG.hdr_dir[LATEST_HDR_NAME]
    HETDEX_DETECT_HDF5_FN = HETDEX_API_CONFIG.detecth5
    HETDEX_DETECT_HDF5_HANDLE = None
    CONT_H5_FN = HETDEX_API_CONFIG.contsourceh5
    CONT_H5_HANDLE = None
    #HETDEX_ELIXER_HDF5_FN = HETDEX_API_CONFIG.elixerh5
    #ELIXER_H5 = None
    #elix_dir = None
    
except Exception as e:
    HETDEX_API_CONFIG = None
    # needed only if detection observered wavelength is not supplied
    try:
        HETDEX_DETECT_HDF5_FN = None #"/work/03946/hetdex/hdr2.1/detect/detect_hdr2.1.h5"
    except:
        HETDEX_DETECT_HDF5_FN = None

    HETDEX_DETECT_HDF5_HANDLE = None

OPEN_DET_FILE = None

SOURCECAT_ROOTPATHS_DICT = {"scratch":"/scratch/projects/hetdex/",
                            "cluster":"/corral-repl/utexas/Hobby-Eberly-Telesco/",
                            "hub":"/home/jovyan/Hobby-Eberly-Telesco/"}
SOURCECAT_SUBPATH = "hdr5/catalogs/"
SOURCECAT_FILE = "source_catalog_5.0.0.h5"
HETDEX_ELIXER_SUBPATH = "hdr5/detect"
HETDEX_ELIXER_FILE = "elixer.h5"

#elix_dir = None

# set up classification dictionary and associated widget
# the widget takes an optional detection list as input either
# as an array of detectids or a text file that can be loaded in
# You may also initiate with no variables to just open any ELiXeR
# on demand


current_wavelength = -1.0

line_id_dict_lae = "1216 LyA"
line_id_dict_lae_wave = 1215.67
line_id_dict_sep = "--------------"
line_id_dict_other = "Other (xxxx)"

line_id_dict_default = "Unk (????)"
line_id_dict = {
    line_id_dict_default: -1.0,
    "1035 OVI": 1034.7625,  #vacuum
    line_id_dict_lae: 1215.67, #vacuum
    "1241 NV": 1240.7077, #vacuum
    "1263 SiII": 1263.40075, #vacuum
    "1398 SiIV": 1397.7617, #vacuum
    "1549 CIV": 1549.4115, #vacuum
    "1640 HeII": 1640.420, #vacuum
    "1909 CIII": 1909.734, #vacuum
    "2326 CII": 2324.095,
    "2799 MgII": 2798.6944,
    "3346 NeV": 3345.821,
    "3426 NeVI": 3425.881,
    "3727 OII": 3727.8,
    "3835 H-eta": 3835.391,
    "3868 NeIII":3868.760,
    "3889 H-zeta": 3889.064,
    "3934 CaII(K)": 3933.6614,
    "3967 NeIII": 3967.470,
    "3968 CaII(H)": 3968.4673,
    "3970 H-epsilon": 3970.079,
    "4101 H-delta": 4101.742,
    "4340 H-gamma": 4340.471,
    "4861 H-beta": 4861.333,
    "4959 OIII": 4958.911,
    "4981 NaI": 4981.3894,
    "5007 OIII": 5006.843,
    "5153 NaI": 5153.4024,
    line_id_dict_sep: -1.0,
    line_id_dict_other: -1.0,
}

classification_labels = ["","AGN","gal","LAB","LzG","meteor","PNe","sat","star","WD","*clear*"]
real_fake_default = "--"
real_fake_sep = "--------------"
real_fake_dict = {
    real_fake_default:0,
    "Confirmed Emission Line": 1,
    "Confirmed Absorption Line": 2,
    "Confirmed Continuum Source": 3,
    real_fake_sep:-999,
    # "True Positive":1,
    "Not Line, Continuum Feature":4, #e.g. return to continuum between two absorption features
    "False Positive (Noise)":10,
    "Artifact (generic)": 11,
    "Artifact (interference)": 12,
    "Artifact (sky line)": 13,
    "Artifact (bad sky sub)": 14,
    "Artifact (cosmic ray)": 15,
    "Artifact (charge trap)": 16,
    "Artifact (stuck pixel)": 17,
    "Artifact (hot pixel/column)": 18,
    "Artifact (diffraction)": 19,

}



class ElixerWidget:
    def __init__(
        self,
        detectfile=None,
        detectlist=None,
        savedfile=None,
        #outfile=None, #redundant with savedfile
        resume=False,
        img_dir=None,
        counterpart=False,
        detect_h5=None,
        elixer_h5=None,
        cutoutpath=None,
        show_cls_buttons=ALLOW_CLASSIFICATION_BUTTONS
    ):

        #global elix_dir, HETDEX_DETECT_HDF5_FN, HETDEX_ELIXER_HDF5_FN
        global HETDEX_DETECT_HDF5_FN, ALLOW_CLASSIFICATION_BUTTONS

        ALLOW_CLASSIFICATION_BUTTONS = show_cls_buttons

        self.elixer_conn_mgr = sql.ConnMgr()
        self.current_idx = 0
        self.show_counterpart_btns = counterpart
        self.detectid = None
        self.shotid = None
        self.ra = None
        self.dec = None
        self.constructor_status = ""
        self.flag = []
        self.vis_class = []
        self.z = []
        self.comment = []
        self.counterpart = []
        self.cutoutpath = cutoutpath

        self.HETDEX_ELIXER_HDF5_FN = None

        self.detections_interface_lines = None
        self.detections_interface_cont = None
        self.source_catalog_detid_revision = -1

        self.sourcecat_root_path = None
        self.source_catalog_h5 = None
        self.sourcecat_z_hetdex = None
        self.sourcecat_z_hetdex_src = None
        self.sourcecat_z_hetdex_conf = None
        self.sourcecat_class_labels = ""
        self.sourcecat_obswave = None
        self.sourcecat_source_id = None
        self.sourcecat_parent_detectid  = None



        # try:
        #     self.HETDEX_API_CONFIG.elixerh5
        # except:
        #     pass

        #setup for source catalog manual updates
        try:
            self.sct = source_catalog_manual_updates.SrcCatUpdateTable()
            self.sct.set_paths()
        except:
            self.sct = None


        try:
            if os.path.isfile(os.path.join(SOURCECAT_ROOTPATHS_DICT['scratch'], SOURCECAT_SUBPATH, SOURCECAT_FILE)):
                h5fn = os.path.join(SOURCECAT_ROOTPATHS_DICT['scratch'], SOURCECAT_SUBPATH, SOURCECAT_FILE)
                self.sourcecat_root_path = SOURCECAT_ROOTPATHS_DICT['scratch']
            elif os.path.isfile(os.path.join(SOURCECAT_ROOTPATHS_DICT['hub'], SOURCECAT_SUBPATH, SOURCECAT_FILE)):
                h5fn = os.path.join(SOURCECAT_ROOTPATHS_DICT['hub'], SOURCECAT_SUBPATH, SOURCECAT_FILE)
                self.sourcecat_root_path = SOURCECAT_ROOTPATHS_DICT['hub']
            elif os.path.isfile(os.path.join(SOURCECAT_ROOTPATHS_DICT['cluster'], SOURCECAT_SUBPATH, SOURCECAT_FILE)):
                h5fn = os.path.join(SOURCECAT_ROOTPATHS_DICT['cluster'], SOURCECAT_SUBPATH, SOURCECAT_FILE)
                self.sourcecat_root_path = SOURCECAT_ROOTPATHS_DICT['cluster']
            else:
                h5fn = None
                self.sourcecat_root_path = None

            if h5fn is not None:
                try:
                    self.source_catalog_h5 = tables.open_file(h5fn)
                except:
                    # note: status box note created yet
                    print(f"Cannot open: {h5fn}")
                    self.source_catalog_h5 = None
            else:  # note: status box note created yet
                print(f"Does not exist: {h5fn}")
                self.source_catalog_h5 = None
        except:
            self.source_catalog_h5 = None

        self.elix_dir = None
        self.ELIXER_H5 = None
        
        if img_dir is not None:
            if op.exists(img_dir):
                self.elix_dir = img_dir
                self.elixer_conn_mgr.extra_db_paths += self.elix_dir
                # also prepend to the database search directories
                # so will look there for alternate databases
                # for key in sql.DICT_DB_PATHS.keys():
                #     if self.elix_dir not in sql.DICT_DB_PATHS[key]:
                #         sql.DICT_DB_PATHS[key].insert(0, self.elix_dir)



        if detect_h5 is not None:
            if op.isfile(detect_h5):
                HETDEX_DETECT_HDF5_FN = detect_h5

        if elixer_h5 is not None:
            if op.isfile(elixer_h5):
                self.HETDEX_ELIXER_HDF5_FN = elixer_h5


        #change order, attempt to load the savedfile first. If it exists, ignore the detectfile and detectlist,
        #if it does not exist, create it and use detectfile or detectlist if they exist

        if savedfile:  # it is specified
            try:
                if not op.exists(savedfile):
                    #it does not exist, though, so we will assume it as the new file to save, but not load from it
                    self.constructor_status += f"savedfile ({savedfile}) does not exist.\n"
                else:
                    saved_data = ascii.read(savedfile)
                    try:
                        self.detectid = np.array(saved_data["detectid"], dtype=int)
                    except:
                        self.detectid = np.array(saved_data["detectid"], dtype=float).astype(int)
                    try:
                        u_detectid = np.unique(self.detectid)
                        if len(self.detectid) != len(u_detectid):
                            self.constructor_status += "!!! WARNING. Provided DetectIDs are NOT UNIQUE !!!\nWill FORCE to UNIQUE"
                            self.detectid = u_detectid
                    except Exception as e:
                        self.constructor_status += "Exception enforcing uniqueness for detectIDs\n"
                        self.constructor_status += str(e)

                    try:
                        self.vis_class = np.array(saved_data["vis_class"], dtype=int)
                    except:
                        self.vis_class = np.zeros(np.size(self.detectid), dtype=int)

                    # could have a flag
                    try:
                        self.flag = np.array(saved_data["flag"], dtype=int)
                    except:
                        self.flag = np.zeros(np.size(self.detectid), dtype=int)

                    # could have z
                    try:
                        self.z = np.array(saved_data["z"], dtype=float)
                    except:
                        self.z = np.full(np.size(self.detectid), -1.0)

                    # could have comment
                    try:
                        self.comment = np.array(
                            saved_data["comments"], dtype="|S80"
                        ).astype(str)
                    except:
                        self.comment = np.full(
                            np.size(self.detectid), "?", dtype="|S80"
                        ).astype(str)

                    # could have counterpart
                    try:
                        self.counterpart = np.array(saved_data["counterpart"], dtype=int)
                    except:
                        self.counterpart = np.full(np.size(self.detectid), -1, dtype=int)
            except:
                print(
                    "Could not open and read in savedfile. Are you sure its in astropy table format"
                )

        #if we have already specified the detection ID list from a saved file, ignore then detectfile and detectlist
        #otherwise attempt to lost the file first and if that files then the list

        #if they both exist, print a warning

        if self.detectid is None or len(self.detectid) == 0:
            if detectfile is not None and len(detectfile) > 0: #assume just a list of detectIDs
                try:
                    self.detectid = np.loadtxt(detectfile, dtype=np.int64, ndmin=1)

                    try:
                        u_detectid = np.unique(self.detectid)
                        if len(self.detectid) != len(u_detectid):
                            self.constructor_status += "!!! WARNING. Provided DetectIDs are NOT UNIQUE !!!\nWill FORCE to UNIQUE"
                            self.detectid = u_detectid
                    except Exception as e:
                        self.constructor_status += "Exception enforcing uniqueness for detectIDs\n"
                        self.constructor_status += str(e)

                    self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
                    self.flag = np.zeros(np.size(self.detectid), dtype=int)
                    self.z = np.full(np.size(self.detectid), -1.0)
                    self.comment = np.full(np.size(self.detectid), "?", dtype="|S80").astype(
                        str
                    )
                    self.counterpart = np.full(np.size(self.detectid), -1, dtype=int)
                    # hidden flag, distinguish vis_class 0 as unset vs reviewed & fake
                    # and possible future use as followup

                    if detectlist is not None and len(detectlist) > 0:
                        self.constructor_status += "detectfile and detectlist both specified. Using detectfile.\n"
                except Exception as e:
                    print (f"Failed to load detectfile {detectfile}. Possible bad format. Expect one detectid per line.\n")
                    print(e)
                    return

            elif detectlist is None:
                global HETDEX_DETECT_HDF5_HANDLE

                if HETDEX_DETECT_HDF5_HANDLE is None:
                    try:
                        print( f"Reading from ({HETDEX_DETECT_HDF5_FN}). ")
                        HETDEX_DETECT_HDF5_HANDLE = tables.open_file(
                            HETDEX_DETECT_HDF5_FN, "r"
                        )
                    except:
                        print("Could not open detections database")

                if HETDEX_DETECT_HDF5_HANDLE is not None:
                    print(f"Attempting to load all detections. This may take a while ...")
                    #these must already be unique so no need to waste time through np.unique
                    self.detectid = HETDEX_DETECT_HDF5_HANDLE.root.Detections.cols.detectid[:]
                    self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
                    self.flag = np.zeros(np.size(self.detectid), dtype=int)
                    self.z = np.full(np.size(self.detectid), -1.0)
                    self.comment = np.full(
                        np.size(self.detectid), "?", dtype="|S80"
                    ).astype(str)
                    self.counterpart = np.full(np.size(self.detectid), -1, dtype=int)

                else:
                    self.detectid = []
            else: #detectlist is specified
                self.detectid = np.array(detectlist).flatten()
                try:
                    u_detectid = np.unique(self.detectid)
                    if len(self.detectid) != len(u_detectid):
                        self.constructor_status += "!!! WARNING. Provided DetectIDs are NOT UNIQUE !!!\nWill FORCE to UNIQUE"
                        self.detectid = u_detectid
                except Exception as e:
                    self.constructor_status += "Exception enforcing uniqueness for detectIDs\n"
                    self.constructor_status += str(e)

                self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
                self.flag = np.zeros(np.size(self.detectid), dtype=int)
                self.z = np.full(np.size(self.detectid), -1.0)
                self.comment = np.full(np.size(self.detectid), "?", dtype="|S80").astype(
                    str
                )
                self.counterpart = np.full(np.size(self.detectid), -1, dtype=int)
        else: #we loaded from a savedfile
            if detectfile is not None and len(detectfile) > 0:
                self.constructor_status += f"detectfile ({detectfile}) ignored in favor of savedfile ({savedfile}).\n"

            if detectlist is not None and len(detectlist) > 0:
                self.constructor_status += f"detectlist ignored in favor of savedfile ({savedfile}).\n"
        # store outfile name if given
        # if outfile:
        #     print(
        #         "Careful with this option, it likely won't work properly. You are better off using the savedfile option"
        #     )
        #     self.outfilename = outfile
        #
        if savedfile:
            self.outfilename = savedfile
        else:
            self.outfilename = "elixer_cls.dat"

        #can be duplicates, but will this wreck the saved file??
        #duplicates mess up the on next or on previous logic which advances/retreats +/- 1 from the matched
        #detectid ... so if it hits a duplicate, you get stuck

        try:
            if self.constructor_status is not None and len(self.constructor_status) > 0:
                print(self.constructor_status)
        except:
            pass

        self.resume = resume
        self.setup_widget()

        interact(self.main_display, x=self.detectbox)

    # def interact(self):
    #     self.setup_widget()
    #     interact(self.main_display, x=self.detectbox)

    def main_display(self, x):

        #global ELIXER_H5, HETDEX_ELIXER_HDF5_FN


        try:
            detectid = np.int64(x)
        except:
            self.status_box.value = f"Invalid entry in DetectID: {x}"
            detectid = 0

        show_selection_buttons = ALLOW_CLASSIFICATION_BUTTONS

        if ALLOW_MANUAL_CATALOG_UPDATE:
            self.updatecat_box = widgets.Output()#layout={"border": "1px solid black"})
            display(self.updatecat_box)

        try:
            objnum = np.where(self.detectid == detectid)[0][0]
            print(
                "On ELiXer Report "
                + str(objnum + 1)
                + "/"
                + str(np.size(self.detectid))
            )
        except:
            print(
                f"Current object {self.detectbox.value} not in original list. Go to Next or Previous DetectID to return to input Detectlist"
            )
            show_selection_buttons = False
            objnum = None

        self.previousbutton.on_click(self.on_previous_click)
        self.nextbutton.on_click(self.on_next_click)
        self.elixerNeighborhood.on_click(self.on_elixer_neighborhood)

        self.updateCatalog.on_click(self.on_update_catalog)
        self.doUpdateCatalog.on_click(self.on_do_update_catalog)

        #print("main_display", objnum)
        self.reset_widget_values(objnum)

        if show_selection_buttons:

            if ALLOW_MANUAL_CATALOG_UPDATE:
                display(
                    widgets.HBox(
                        [
                            self.previousbutton,
                            self.nextbutton,
                            self.elixerNeighborhood,
                            # self.line_id_drop,
                            # self.wave_box,
                            # self.z_box,
                            self.z_hetdex_box,
                            self.comment_box,
                            self.updateCatalog,
                        ]
                    )
                )
            else:
                display(
                    widgets.HBox(
                        [
                            self.previousbutton,
                            self.nextbutton,
                            self.elixerNeighborhood,
                            # self.line_id_drop,
                            # self.wave_box,
                            # self.z_box,
                            self.z_hetdex_box,
                            self.comment_box,
                            #self.updateCatalog,
                        ]
                    )
                )

            if self.show_counterpart_btns:
                display(
                    widgets.HBox(
                        [
                            widgets.Label(value="Counterpart:"),
                            self.c_none_button,
                            self.c_aper_button,
                            self.c_blue_button,
                            self.c_red_button,
                            self.c_green_button,
                            self.c_other_button,
                        ]
                    )
                )

                self.c_none_button.on_click(self.c_none_button_click)
                self.c_aper_button.on_click(self.c_aper_button_click)
                self.c_blue_button.on_click(self.c_blue_button_click)
                self.c_red_button.on_click(self.c_red_button_click)
                self.c_green_button.on_click(self.c_green_button_click)
                self.c_other_button.on_click(self.c_other_button_click)

            display(
                widgets.HBox(
                    [
                        self.sm1_button,
                        self.lowz_button,
                        self.other_button,
                        self.s1_button,
                        self.s2_button,
                        self.s3_button,
                        self.s4_button,
                        self.s5_button,
                    ]
                )
            )

            self.sm1_button.on_click(self.sm1_button_click)
            # self.s0_button.on_click(self.s0_button_click)
            self.s1_button.on_click(self.s1_button_click)
            self.s2_button.on_click(self.s2_button_click)
            self.s3_button.on_click(self.s3_button_click)
            self.s4_button.on_click(self.s4_button_click)
            self.s5_button.on_click(self.s5_button_click)
            self.lowz_button.on_click(self.lowz_button_click)
            self.other_button.on_click(self.other_button_click)
        else:

            if ALLOW_MANUAL_CATALOG_UPDATE:
                display(
                    widgets.HBox(
                        [
                            self.previousbutton,
                            self.nextbutton,
                            self.elixerNeighborhood,
                            # self.line_id_drop,
                            # self.wave_box,
                            # self.z_box,
                            self.z_hetdex_box,
                            self.comment_box,
                            self.updateCatalog,
                        ]
                    )
                )
            else:
                display(
                    widgets.HBox(
                        [
                            self.previousbutton,
                            self.nextbutton,
                            self.elixerNeighborhood,
                            # self.line_id_drop,
                            # self.wave_box,
                            # self.z_box,
                            self.z_hetdex_box,
                            self.comment_box,
                            #self.updateCatalog,
                        ]
                    )
                )

        try:
            try:
                # display(Image(sql.fetch_elixer_report_image(sql.get_elixer_report_db_path(detectid),detectid)))
                # display(Image(sql.fetch_elixer_report_image(self.elixer_conn_mgr.get_connection(detectid),detectid)))

                 #added 2023-11-02 by EMC   
                if self.cutoutpath is not None:
                    display(Image(filename=op.join(self.cutoutpath, '{}.png'.format(detectid))))

                display(Image(self.elixer_conn_mgr.fetch_image(detectid)))

            except Exception as e:

                if self.elix_dir:
                    fname = op.join(self.elix_dir, "%d.png" % (detectid))
                    if op.exists(fname):
                        display(Image(fname))
                    else:  # try the archive location
                        self.status_box.value = (
                            "Report databases not found and specified img_dir does not exist\n"
                            + str(e)
                            + "\n"
                            + traceback.format_exc()
                            + "\n"
                        )
                        # print("Cannot load ELiXer Report image: ", fname)
                else:
                    # dbpaths = ""
                    # for key in sql.DICT_DB_PATHS.keys():
                    #     dbpaths += str(sql.DICT_DB_PATHS[key])

                    self.status_box.value = str(e) + "\n" + traceback.format_exc() #+ "\n" + dbpaths


                    # print("Cannot load ELiXer Report image: ", str(detectid))
        except Exception as e:
            self.status_box.value = (
                "Report databases not found and no other img__dir not specified\n"
                + str(e)
                + "\n"
                + traceback.format_exc()
            )
            # pass
            # print("Cannot load ELiXer Report image: ", str(detectid))

        if self.ELIXER_H5 is None:
            if self.HETDEX_ELIXER_HDF5_FN is not None and op.exists(self.HETDEX_ELIXER_HDF5_FN):
                try:
                    self.ELIXER_H5 = tables.open_file(self.HETDEX_ELIXER_HDF5_FN, "r")
                except Exception as e:
                    self.status_box.value = str(e) + "\n" + traceback.format_exc()
                    self.ELIXER_H5 = None
            else:
                pass
                # print('No counterparts found in ' + self.HETDEX_ELIXER_HDF5_FN)
                # return

        # only execute the below if we have ELIXER_H5 ... the return just above exits this func otherwise
        if self.ELIXER_H5:
            detectid_i = detectid

            try:
                self.CatalogMatch = self.ELIXER_H5.root.CatalogMatch.read_where(
                    "detectid == detectid_i"
                )

                if np.size(self.CatalogMatch) == 1:
                    display(
                        widgets.HBox(
                            [
                                widgets.Label(value="Extract Counterpart:  "),
                                self.e_blue_button,
                            ]
                        )
                    )
                elif np.size(self.CatalogMatch) == 2:
                    display(
                        widgets.HBox(
                            [
                                widgets.Label(value="Extract Counterpart:  "),
                                self.e_blue_button,
                                self.e_red_button,
                            ]
                        )
                    )
                elif np.size(self.CatalogMatch) == 3:
                    display(
                        widgets.HBox(
                            [
                                widgets.Label(value="Extract Counterpart:  "),
                                self.e_blue_button,
                                self.e_red_button,
                                self.e_green_button,
                            ]
                        )
                    )
                else:
                    pass
            except Exception as e:
                self.status_box.value = str(e) + "\n" + traceback.format_exc()
                self.CatalogMatch = []

            #        display(widgets.HBox([widgets.Label(value="Manual Entry:  "),
            #                              self.e_manual_ra,
            #                              self.e_manual_dec,
            #                              self.e_manual_button]))

            self.e_blue_button.on_click(self.e_blue_button_click)
            self.e_red_button.on_click(self.e_red_button_click)
            self.e_green_button.on_click(self.e_green_button_click)
            self.e_manual_button.on_click(self.e_manual_button_click)

        self.det_table_button.on_click(self.det_table_button_click)
        self.get_mini_button.on_click(self.get_mini_button_click)
        self.get_sdss_button.on_click(self.get_sdss_button_click)
        self.get_legacy_button.on_click(self.get_legacy_button_click)

        display(widgets.HBox( [self.det_table_button,
                               self.get_sdss_button,
                               self.get_legacy_button,
                               self.get_mini_button]))



        # always show the status box
        display(self.status_box)

        self.bottombox = widgets.Output(layout={"border": "1px solid black"})
        display(self.bottombox)




    def setup_widget(self):
        if self.resume:
            try:
                i_start = np.max(np.where(self.flag != 0)) + 1

                if i_start is None:
                    i_start = 0
                    detectstart = self.detectid[i_start]
                elif i_start < np.size(self.detectid):
                    detectstart = self.detectid[i_start]
                else:
                    i_start = 0
                    detectstart = self.detectid[i_start]  # np.min(self.detectid)
            except:
                i_start = 0
                detectstart = self.detectid[i_start]

        else:
            i_start = 0
            detectstart = self.detectid[i_start]  # np.min(self.detectid)

        self.current_idx = i_start

        # self.detectbox = widgets.BoundedIntText(
        #     value=detectstart,
        #     # min=1,
        #     min=1000000000,
        #     max=int(8.9e18),#100000000000, #99999999999,
        #     step=1,
        #     description="DetectID:",
        #     disabled=False,
        # )

        self.detectbox = widgets.Text(
            value=str(detectstart),
            placeholder = str(detectstart),
            # min=1,
            description="DetectID:",
            disabled=False,
        )

        self.detid = np.int64(self.detectbox.value)
        self.neighbor_list = []
        
        self.previousbutton = widgets.Button(
            layout=Layout(width="5%")
        )  # description='Previous DetectID')
        self.nextbutton = widgets.Button(
            layout=Layout(width="5%")
        )  # description='Next DetectID')

        # see https://fontawesome.com/icons?d=gallery
        self.previousbutton.style.button_color = "darkgray"
        self.previousbutton.icon = "arrow-circle-left"
        self.nextbutton.style.button_color = "darkgray"
        self.nextbutton.icon = "arrow-circle-right"
        # self.nextbutton.layout = Layout()

        self.elixerNeighborhood = widgets.Button(
            description="Neighbors", layoout=Layout(width="10%")
        )  # , button_style='info')
        self.elixerNeighborhood.style.button_color = "darkgray"
        # self.detectwidget = widgets.HBox([self.detectbox, self.nextbutton])


        self.s1_button = widgets.Button(
            description="Noise",
            tooltip="Very sure (90%+ confident) this is just noise (or nothing)",
            button_style="success",
        )

        self.updateCatalog = widgets.Button(
            #description="Update Catalog", layoout=Layout(width="10%")
            description = "Update Catalog",
            tooltip = "Manually update catalog info for this detection",
            button_style = "success",
        )  # , button_style='info')
        self.updateCatalog.style.button_color = "black"
        # self.updateCatalog.add_class("left-spacing-class")
        # display(HTML(
        #     "<style>.left-spacing-class {margin-left: 10px;}</style>"
        # ))

        self.doUpdateCatalog = widgets.Button(
            #description="Update Catalog", layoout=Layout(width="10%")
            description = "Update",
            tooltip = "Manually update catalog info for this detection",
            button_style = "success",
        )
        self.doUpdateCatalog.style.button_color = "black"
        self.doUpdateCatalog.add_class("left-spacing-class")
        display(HTML(
            "<style>.left-spacing-class {margin-left: 50px;}</style>"
        ))

        # self.class_labels_drop = widgets.Dropdown(
        #     options=classification_labels,
        #     value="",
        #     description="Label",
        #     layout=Layout(width="17%"),
        #     disabled=False,
        # )

        self.class_labels_drop = widgets.Combobox(
            options=classification_labels,
            placeholder='labels',
            ensure_option=False,
            description="Label",
            layout=Layout(width="25%"),
            disabled=False,
        )

        #self.class_labels_drop.observe(self._handle_class_labels_selection, names="value")

        self.real_fake_drop = widgets.Dropdown(
            options=real_fake_dict.keys(),
            value=real_fake_default,
            description="Conf",
            layout=Layout(width="30%"),
            disabled=False,
        )

        #self.real_fake_drop.observe(self._handle_real_fake_selection, names="value")

        self.line_id_drop = widgets.Dropdown(
            options=line_id_dict.keys(),
            value=line_id_dict_default,
            description="Line",
            layout=Layout(width="20%"),
            disabled=False,
        )

        self.line_id_drop.observe(self._handle_line_id_selection, names="value")

        #restwave
        self.wave_box = widgets.FloatText(
            value=-1.0,
            step=0.00001,
            description= "wave rest", #r"$\lambda$ rest",
            layout=Layout(width="20%"),
            disabled=False,
            indent=False,
        )
        # self.wave_box = widgets.Text(
        #     value='-1.0',
        #     description=r"$\lambda$ rest",
        #     layout=Layout(width="20%")
        #     disabled=False),

        self.wave_box.observe(self._handle_wave_box_selection, names="value")
        self.z_box = widgets.FloatText(
            value=-1.0,
            step=0.00001,
            description="z",
            layout=Layout(width="20%"),
            disabled=False,
            indent=False,
        )
        #self.z_box.add_class("left-spacing-class")
        # display(HTML("<style>.left-spacing-class {margin-left: 10px;}</style>"))

        #show the z_hetdex value AND its source (z_hetdex, z_hetdex_soruce) as combined string
        self.z_hetdex_box = widgets.Text(
            value="",
            placeholder="-1 N/A",
            description="z_hetdex:",
            disabled=True, #this is a display ONLY not user editable
            layout=Layout(width='42%')
        )

        self.comment_box = widgets.Text(
            value="",
            placeholder="Enter any comments here",
            description="Comments:",
            disabled=False,
            layout=Layout(width='40%')
        )

        # self.clusterid_box = widgets.BoundedIntText(
        #     value=0,
        #     min=0,
        #     max=int(9e9),
        #     description="ClusterID:",
        #     layout=Layout(width="17%"),
        #     disabled=False,
        #     #layout=Layout(width='50%')
        # )

        self.clusterid_box = widgets.Text(
            value="",
            description="Parent ID:",
            layout=Layout(width="20%"),
            disabled=False,
            #layout=Layout(width='50%')
        )

        self.catalog_sourceid_box = widgets.Text(
            value="",
            description="Source ID:",
            layout=Layout(width="20%"),
            disabled=True,
            #layout=Layout(width='50%')
        )

        self.catalog_comment_box = widgets.Text(
            value="",
            placeholder="Enter any comments here",
            description="Comments:",
            disabled=False,
            layout=Layout(width='79%')
        )

        self.catalog_status_box = widgets.Textarea(
            value="",
            placeholder= "",
            description="Status:",
            disabled=False,
            layout=Layout(width="80%"),
        )

        self.status_box = widgets.Textarea(
            value="",
            placeholder= "OK",
            description="Status:",
            disabled=False,
            layout=Layout(width="80%"),
        )

        if self.constructor_status is not None and len(self.constructor_status) > 0:
            self.status_box.value = self.constructor_status
            self.constructor_status = None


        # buttons as classification selection
        # self.s0_button = widgets.Button(description=' No Imaging ', button_style='success')
        # self.s1_button = widgets.Button(description=' Spurious ', button_style='success')
        # self.s2_button = widgets.Button(description=' Not LAE ', button_style='success')
        # self.s3_button = widgets.Button(description=' Unknown ', button_style='success')
        # self.s4_button = widgets.Button(description=' Maybe LAE ', button_style='success')
        # self.s5_button = widgets.Button(description=' Definite LAE ', button_style='success')

        self.sm1_button = widgets.Button(
            description="Artifact",
            button_style="danger",
            tooltip="Bad pixels, wavy lines, streaks, or other problems",
        )
        # layout=Layout(width='10%'))

        self.other_button = widgets.Button(
            description="Star/Meteor",
            button_style="warning",
            tooltip="Star, meteor or other non-galaxy object",
        )
        self.lowz_button = widgets.Button(
            description="Nearby Galaxy",
            button_style="warning",
            tooltip="Nearby galaxy (large and/or has apparent structure",
        )

        # self.s0_button = widgets.Button(description='  Not LAE (0) ', button_style='success')
        self.s1_button = widgets.Button(
            description="Noise",
            tooltip="Very sure (90%+ confident) this is just noise (or nothing)",
            button_style="success",
        )
        self.s1_button.style.button_color = "blue"
        self.s2_button = widgets.Button(
            description="Likely Noise",
            tooltip="Somewhat sure (70%-90% confident) this is just noise (or nothing)",
            button_style="success",
        )
        self.s2_button.style.button_color = "DodgerBlue"
        self.s3_button = widgets.Button(
            description="Maybe Distant",
            tooltip="Not sure, but more likely a distant galaxy than noise",
            button_style="success",
        )
        self.s3_button.style.button_color = "CadetBlue"
        self.s4_button = widgets.Button(
            description="Likely Distant",
            tooltip="Somewhat sure (70%-90% confident) this is a distant galaxy",
            button_style="success",
        )
        self.s4_button.style.button_color = "MediumSeaGreen"
        self.s5_button = widgets.Button(
            description="Distant Galaxy",
            tooltip="Very sure (90%+ confident) this is a distant galaxy",
            button_style="success",
        )
        self.s5_button.style.button_color = "green"

        self.c_none_button = widgets.Button(
            description="Not visible"
        )  # , button_style='success')
        self.c_none_button.style.button_color = "lightgray"
        self.c_aper_button = widgets.Button(
            description="Aperture"
        )  # , button_style='success')
        self.c_aper_button.style.button_color = "gold"
        self.c_blue_button = widgets.Button(
            description="Blue"
        )  # , button_style='success')
        self.c_blue_button.style.button_color = "blue"
        self.c_red_button = widgets.Button(
            description="Red"
        )  # , button_style='success')
        self.c_red_button.style.button_color = "red"
        self.c_green_button = widgets.Button(
            description="Green"
        )  # , button_style='success')
        self.c_green_button.style.button_color = "green"
        self.c_other_button = widgets.Button(
            description="Other"
        )  # , button_style='success')
        self.c_other_button.style.button_color = "orange"

        self.e_blue_button = widgets.Button(
            description="Blue"
        )  # , button_style='success')
        self.e_blue_button.style.button_color = "blue"
        self.e_red_button = widgets.Button(
            description="Red"
        )  # , button_style='success')
        self.e_red_button.style.button_color = "red"
        self.e_green_button = widgets.Button(
            description="Green"
        )  # , button_style='success')
        self.e_green_button.style.button_color = "green"
        self.e_manual_ra = widgets.FloatText(
            value=0.0, description="RA (deg):", layout=Layout(width="20%")
        )
        self.e_manual_dec = widgets.FloatText(
            value=0.0, description="DEC (deg):", layout=Layout(width="20%")
        )
        self.e_manual_button = widgets.Button(description="Go")

        self.det_table_button = widgets.Button(
            description="Get Detection Table Info", layout=Layout(width="30%")
        )

        self.get_mini_button = widgets.Button(
            description="Get the zooniverse mini Image", layout=Layout(width="30%"))

        self.get_sdss_button = widgets.Button(
            description="Open SDSS Nav Tool",
            layout=Layout(width="30%"))

        self.get_legacy_button = widgets.Button(
            description="Open legacysurvey.com",
            layout=Layout(width="30%"))
        
        # self.submitbutton = widgets.Button(description="Submit Classification", button_style='success')
        # self.savebutton = widgets.Button(description="Save Progress", button_style='success')

    def get_observed_wavelength(self):
        global HETDEX_DETECT_HDF5_HANDLE, HETDEX_DETECT_HDF5_FN, current_wavelength

        if self.sourcecat_obswave is not None:
            return self.sourcecat_obswave

        current_wavelength = -1.0

        if (HETDEX_DETECT_HDF5_HANDLE is None) and (HETDEX_DETECT_HDF5_FN is not None):
            try:
                HETDEX_DETECT_HDF5_HANDLE = tables.open_file(HETDEX_DETECT_HDF5_FN, "r")
            except Exception as e:
                self.status_box.value = str(e) + "\n" + traceback.format_exc()
                # pass

        if HETDEX_DETECT_HDF5_HANDLE:
            try:
                #dtb = HETDEX_DETECT_HDF5_HANDLE.root.Detections
                try:
                    q_detectid = np.int64(self.detectbox.value)
                except:
                    q_detectid = 0
                    self.status_box.value = f"Invalid entry in DetectID box"

                rows =  HETDEX_DETECT_HDF5_HANDLE.root.Detections.read_where("detectid==q_detectid", field="wave")
                if (rows is not None) and (rows.size == 1):
                    current_wavelength = rows[0]
            except Exception as e:
                self.status_box.value = str(e) + "\n" + traceback.format_exc()

        return current_wavelength

    def _handle_wave_box_selection(self, event):

        global current_wavelength, line_id_dict_lae_wave

        _old_wave = event["old"]
        _new_wave = event["new"]

        if _old_wave == _new_wave:
            return

        if _new_wave < 0:  # i.e. a -1.00
            # do not change the z and lambda_rest values
            self.z_box.value = -1.0
            self.wave_box.value = -1.0
        else:
            if current_wavelength < 0:  # need to find it
                self.get_observed_wavelength()

            #self.z_box.value = np.round(current_wavelength / _new_wave - 1.0, 5)
            #self.z_box.value = np.round(current_wavelength / _new_wave - 1.0, 5)
            z = current_wavelength / _new_wave - 1.0
            #def z_correction(z,w_obs,vcor=None,shotid=None):#,*args):
            self.z_box.value = SU.z_correction(z,current_wavelength,shotid=self.shotid)

            self.wave_box.value = np.round(_new_wave, 5)

            # do NOT reset the line id label ... can cause a loop
            # self.get_line_match(self.z_box.value, current_wavelength )


    def _handle_class_labels_selection(self, event):
        pass

    def _handle_real_fake_selection(self, event):
        pass

    def _handle_line_id_selection(self, event):
        """
        event': 'value',
        'old': 'Unk (????)',
        'new': '4960 OIII',
        'owner': Dropdown(description='Line ID', index=13, options=('Unk (????)', 'Other (xxxx)', '1216 LyA', '1241 NV', '1260 SiII', '1549 CIV', '1640 HeII', '1909 CII', '2799 MgII', '3727 OII', '4102 H-delta', '4342 H-gamma', '4863 H-beta', '4960 OIII', '5008 OIII'), value='4960 OIII'),
        'type': 'change'
        :param event:
        :return:
        """
        global current_wavelength, line_id_dict_lae_wave

        _old = event["old"]
        _old_wave = line_id_dict[_old]
        _new = event["new"]
        _new_wave = line_id_dict[_new]

        if _old_wave == _new_wave:
            return

        if _new_wave < 0:  # i.e. a -1.00
            # do not change the z and lambda_rest values
            if _new == line_id_dict_default:
                self.z_box.value = -1.0
                self.wave_box.value = -1.0
        else:
            if current_wavelength < 0:  # need to find it
                self.get_observed_wavelength()

            #self.z_box.value = np.round(current_wavelength / _new_wave - 1.0, 5)
            z = current_wavelength / _new_wave - 1.0
            #def z_correction(z,w_obs,vcor=None,shotid=None):#,*args):
            self.z_box.value = SU.z_correction(z,current_wavelength,shotid=self.shotid)

            self.wave_box.value = np.round(_new_wave, 5)

            # if (self.z_box.value < 0.0) and (self.line_id_drop.value != line_id_dict_default):
            #     self.z_box.style = 'warning'

            # what button value to hit?
            # LAE lines would be LyA, CIV, HeII
            # set these to '5' or '4' or what?
            # This causes some syncrhonization problems ... probably should skip it
            # if abs(self.wave_box.value -line_id_dict_lae_wave) < 0.1:
            #     if self.s5_button.icon == '': #only trip the button if not already at 5
            #         self.s5_button_click(None)

        # causes dirty flag
        # self.current_idx = ix
        # self.detectbox.value = self.detectid[ix]

        # print("Drop down event",event)
        # print(event['old'],event['new'],line_id_dict[event['new']])

    def goto_previous_detect(self):

        try:
            if np.int64(self.detectbox.value) in self.detectid:
                ix = np.where(self.detectid == np.int64(self.detectbox.value))[0][0]
                old_id = self.detectid[ix]
                ix -= 1
                while self.detectid[ix] == old_id and ix < len(self.detectid):
                    ix -= 1
            else:
                ix = np.max(np.where(self.detectid <= np.int64(self.detectbox.value)))

            if ix < 0:
                ix = 0
                print("At the beginning of the DetectID List")
                return
        except:
            # invalid index ... the report displayed is not in the operating list
            # so use the last good index:
            ix = self.current_idx



        # causes dirty flag
        self.current_idx = ix
        self.detectbox.value = str(self.detectid[ix])
        self.detid = np.int64(self.detectbox.value)
        # print("Prev detect idx", ix)

        self.reset_widget_values(idx=ix)

        # print("Prev detect post reset idx", )
        
    def goto_next_detect(self):
        try:

            if np.int64(self.detectbox.value) in self.detectid:
                ix = np.where(self.detectid == np.int64(self.detectbox.value))[0][0]
                old_id = self.detectid[ix]
                ix += 1
                while ix < len(self.detectid) and self.detectid[ix] == old_id:
                    ix += 1
            else:
                ix = np.where(self.detectid > np.int64(self.detectbox.value))[0][0]

            if ix >= np.size(self.detectid):
                ix = np.argmax(self.detectid)
                print("At the end of the DetectID List")
                self.status_box.value += "At the end of the DetectID List"
                return
        except Exception as e:
            # invalid index ... the report displayed is not in the operating list
            # so use the last good index:
            #print(e)
            self.status_box.value += "EXCEPTION"
            self.status_box.value += str(e)
            ix = self.current_idx

        self.current_idx = ix
        self.detectbox.value = str(self.detectid[ix])
        self.detid = np.int64(self.detectbox.value)
        self.reset_widget_values(idx=ix)

        
    def set_counterpart(self, value=0):
        self.counterpart[self.current_idx] = value
        self.c_none_button.icon = ""
        self.c_aper_button.icon = ""
        self.c_blue_button.icon = ""
        self.c_red_button.icon = ""
        self.c_green_button.icon = ""
        self.c_other_button.icon = ""

        if value == 0:
            self.c_none_button.icon = "check"
        elif value == 1:
            self.c_aper_button.icon = "check"
        elif value == 2:
            self.c_blue_button.icon = "check"
        elif value == 3:
            self.c_red_button.icon = "check"
        elif value == 4:
            self.c_green_button.icon = "check"
        elif value == 99:
            self.c_other_button.icon = "check"

    def set_classification(self, value=0):
        self.current_idx = np.where(self.detectid == np.int64(self.detectbox.value))[0][
            0
        ]  # current position

        self.vis_class[self.current_idx] = value
        self.flag[self.current_idx] = 1

        # make sure the values match
        if self.line_id_drop.value == line_id_dict_default:
            self.wave_box.value = -1.0
            self.z_box.value = -1.0

        self.z[self.current_idx] = self.z_box.value
        self.comment[self.current_idx] = self.comment_box.value

        self.on_save_click(None)
        self.goto_next_detect()

    def get_line_match(self, z, obswave=None, restwave=None):
        global current_wavelength

        if obswave is None and restwave is None:
            self.line_id_drop.value = line_id_dict_default
            return line_id_dict_default, -1

        if z > -0.1:
            self.line_id_drop.value = line_id_dict_other
            if restwave is None and obswave is not None:
                restwave =obswave / (1.0 + z)
        elif restwave is not None and obswave is not None:
            z = obswave/restwave - 1.0
            z = SU.z_correction(z,obswave,shotid=self.shotid)
        else:
            self.line_id_drop.value = line_id_dict_default
            self.wave_box.value = -1.0 #rest_wave

        if z is None or z < 0 or obswave is None or obswave < 0 or restwave is None or restwave < 0:
            return self.line_id_drop.value, self.wave_box.value #rest_wave


        #w_rest = np.round(obswave / (1.0 + z), 5)
        # find the match in the line_id_dict
        best_k = None
        for k in line_id_dict.keys():
            if abs(line_id_dict[k] - restwave) < 6.0:
                if best_k is None:
                    best_k = k
                else:
                    if abs(line_id_dict[k] - restwave) < abs(line_id_dict[best_k] - restwave):
                        best_k = k

        if best_k is not None:
            self.line_id_drop.value = best_k
            self.wave_box.value = restwave

        return self.line_id_drop.value, self.wave_box.value

    def reset_widget_values(self, idx=0):

        global current_wavelength

        # print("Reset idx",idx,"Current w",current_wavelength)

        # DD 2020-06-20
        # given the potential use of low ELiXer assigned detectid values, and the unrealized use of this logic,
        # just comment out this following check for an index vs a detectID
        # literal value is int(1e9), the lowest possible HETDEX detectID
        # if self.detectbox.value < 1000000000: #assume an index
        #     self.detectbox.value = self.detectid[idx]
        #     return

        self.updateCatalog.disabled = False
        #if idx is not None:
        #    self.z_box.value = self.z[idx]  # -1.0

        #reset prior to update
        self.shotid = None
        self.sourcecat_z_hetdex = None
        self.sourcecat_z_hetdex_src = None
        self.sourcecat_z_hetdex_conf = None
        self.sourcecat_class_labels = ""
        self.sourcecat_obswave = None
        self.sourcecat_source_id = None
        self.sourcecat_parent_detectid  = 0
        self.source_catalog_detid_revision = -1

        self.wave_box.value = -1
        self.line_id_drop.value = line_id_dict_default
        self.z_box.value = -1
        self.catalog_status_box.value = ""
        self.catalog_comment_box.value = ""
        self.real_fake_drop.value = real_fake_default
        self.status_box.value = ""

        self.read_source_catalog_basic(detectid=np.int64(self.detectbox.value))

        if self.sourcecat_z_hetdex is not None:
            self.z_box.value = self.sourcecat_z_hetdex


        current_wavelength = self.get_observed_wavelength()
       # print("Line match before", current_wavelength, self.line_id_drop.value, self.wave_box.value, self.z_box.value)

        self.line_id_drop.value, self.wave_box.value = self.get_line_match(self.z_box.value, current_wavelength)
        # if self.sourcecat_z_hetdex is not None:
        #     #the get_line_match may apply z_corrections, but in this case we don't want that
        #     #since it was ALREADY applied
        #     self.z_box.value = self.sourcecat_z_hetdex

        self.wave_box.value = line_id_dict[self.line_id_drop.value]

        if (self.wave_box.value == -1) and (self.sourcecat_obswave is not None) and \
                ((self.sourcecat_z_hetdex is not None) and  (self.sourcecat_z_hetdex > -1)):
            self.wave_box.value = self.sourcecat_obswave/(self.sourcecat_z_hetdex + 1.0)

       # print("Line match after", current_wavelength,self.line_id_drop.value, self.wave_box.value,self.z_box.value)

        if idx is not None:
            self.comment_box.value = self.comment[idx]

        # print("Updated Reset idx", idx, "Current w", current_wavelength)

        # reset all to base
        self.lowz_button.icon = ""
        self.other_button.icon = ""
        self.sm1_button.icon = ""
        # self.s0_button.icon = ''
        self.s1_button.icon = ""
        self.s2_button.icon = ""
        self.s3_button.icon = ""
        self.s4_button.icon = ""
        self.s5_button.icon = ""

        self.c_none_button.icon = ""
        self.c_aper_button.icon = ""
        self.c_blue_button.icon = ""
        self.c_red_button.icon = ""
        self.c_green_button.icon = ""
        self.c_other_button.icon = ""

        # mark the label on the button
        if idx is not None and self.flag[idx] != 0:
            # if self.vis_class[idx] == 0:
            #    self.s0_button.icon = 'check'
            if self.vis_class[idx] == 1:
                self.s1_button.icon = "check"
            elif self.vis_class[idx] == 2:
                self.s2_button.icon = "check"
            elif self.vis_class[idx] == 3:
                self.s3_button.icon = "check"
            elif self.vis_class[idx] == 4:
                self.s4_button.icon = "check"
            elif self.vis_class[idx] == 5:
                self.s5_button.icon = "check"
            elif self.vis_class[idx] == 11:
                self.lowz_button.icon = "check"
            elif self.vis_class[idx] == 12:
                self.other_button.icon = "check"
            elif self.vis_class[idx] == -1:
                self.sm1_button.icon = "check"

            if self.counterpart[idx] == 0:
                self.c_none_button.icon = "check"
            elif self.counterpart[idx] == 1:
                self.c_aper_button.icon = "check"
            elif self.counterpart[idx] == 2:
                self.c_blue_button.icon = "check"
            elif self.counterpart[idx] == 3:
                self.c_red_button.icon = "check"
            elif self.counterpart[idx] == 4:
                self.c_green_button.icon = "check"
            elif self.counterpart[idx] == 99:
                self.c_other_button.icon = "check"

    def on_previous_click(self, b):
        self.goto_previous_detect()

    def on_next_click(self, b):
        self.goto_next_detect()

    def c_none_button_click(self, b):
        self.set_counterpart(0)

    def c_aper_button_click(self, b):
        self.set_counterpart(1)

    def c_blue_button_click(self, b):
        self.set_counterpart(2)

    def c_red_button_click(self, b):
        self.set_counterpart(3)

    def c_green_button_click(self, b):
        self.set_counterpart(4)

    def c_other_button_click(self, b):
        self.set_counterpart(99)

    def sm1_button_click(self, b):
        global line_id_dict_lae
        self.line_id_drop.value = line_id_dict_default
        self.wave_box.value = -1.0
        self.z_box.value = -1.0

        self.set_classification(-1)

    def lowz_button_click(self, b):
        global line_id_dict_lae
        self.set_classification(11)

    def other_button_click(self, b):
        global line_id_dict_lae
        self.set_classification(12)

    # def s0_button_click(self, b):
    #     global line_id_dict_lae
    #     if self.is_consistent_with_lae():
    #         #NOT okay, if you say is not LAE
    #         self.line_id_drop.value == line_id_dict_default
    #         self.wave_box.value = -1.0
    #         self.z_box.value = -1.0
    #     else:
    #         pass #okay, already NOT consistent with LAE
    #
    #     self.set_classification(0)

    def s1_button_click(self, b):
        global line_id_dict_lae
        if self.is_consistent_with_lae():
            # NOT okay, if you say is not LAE
            self.line_id_drop.value = line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0
        else:
            pass  # okay, already NOT consistent with LAE

        self.set_classification(1)

    def s2_button_click(self, b):
        global line_id_dict_lae

        if self.is_consistent_with_lae():
            # NOT okay, if you say is not LAE
            self.line_id_drop.value = line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0
        else:
            pass  # okay, already NOT consistent with LAE

        self.set_classification(2)

    def s3_button_click(self, b):
        # leave the drop box where ever it is
        self.set_classification(3)

    def s4_button_click(self, b):
        global line_id_dict_lae
        if self.is_consistent_with_lae():
            pass  # already okay, could be CIV, etc
        else:  # not consistent with LAE, so reset
            self.line_id_drop.value = line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0

        self.set_classification(4)

    def s5_button_click(self, b):
        global line_id_dict_lae

        if self.is_consistent_with_lae():
            pass  # already okay, could be CIV, etc
        else:  # not consistent with LAE, so reset
            self.line_id_drop.value = line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0

        self.set_classification(5)

    def is_consistent_with_lae(self):
        # is the emission line consistent with this being an LAE in HETDEX
        # basically, a cheat ... is  1.9 < z < 3.5
        if 1.8 < self.z_box.value < 3.6:
            return True
        else:
            return False

    def on_save_click(self, b):
        self.output = Table()
        self.output.add_column(Column(self.detectid, name="detectid", dtype=int))
        self.output.add_column(Column(self.vis_class, name="vis_class", dtype=int))
        self.output.add_column(Column(self.flag, name="flag", dtype=int))
        self.output.add_column(Column(self.z, name="z", dtype=float))
        self.output.add_column(Column(self.counterpart, name="counterpart", dtype=int))
        self.output.add_column(Column(self.comment, name="comments", dtype="|S80"))

        ascii.write(self.output, self.outfilename, overwrite=True)

    def on_elixer_neighborhood(self, b):
        detectid = np.int64(self.detectbox.value)

        #try:
            # display(Image(sql.fetch_elixer_report_image(sql.get_elixer_report_db_path(detectid,report_type="nei"), detectid)))
            # display(Image(sql.fetch_elixer_report_image(self.elixer_conn_mgr.get_connection(detectid,report_type="nei"), detectid)))
        nei_imag = self.elixer_conn_mgr.fetch_image(detectid, report_type="nei")

        if nei_imag is not None:

            try:
                pil_img = PILImage.open(io.BytesIO(nei_imag))
                if 'Neighbors' in pil_img.info.keys():
                    self.neighbor_list = pil_img.info['Neighbors']
                    if len(self.neighbor_list) > 0:
                        self.status_box.value = f"Neighbors: {self.neighbor_list}"
                else:
                    self.neighbor_list = []
                del pil_img
            except: #Exception as e:
                pass #self.status_box.value = str(e) + "\n" + traceback.format_exc()

            with self.bottombox:
                display(
                    #Image(self.elixer_conn_mgr.fetch_image(detectid, report_type="nei"))
                    Image(nei_imag)
                )
        else:#except Exception as e:
            # self.status_box.value = str(e) + "\n" + traceback.format_exc()

            if self.elix_dir:
                path = op.join(self.elix_dir, "%d_nei.png" % (detectid))

                if not op.isfile(path):
                    path = op.join(self.elix_dir, "%dnei.png" % (detectid))  # try w/o '_'
                    if not op.isfile(path):
                        print("%s not found" % path)
                else:
                    try:
                        pil_img = PILImage.open(path)
                        if 'Neighbors' in pil_img.info.keys():
                            self.neighbor_list = pil_img.info['Neighbors']
                            if len(self.neighbor_list) > 0:
                                self.status_box.value = f"Neighbors: {self.neighbor_list}"
                        else:
                            self.neighbor_list = []
                        del pil_img
                    except: # Exception as e:
                        pass #self.status_box.value = str(e) + "\n" + traceback.format_exc()

                    with self.bottombox:
                        display(Image(path))
            else:
                print("neighborhood not found")

    def on_update_catalog(self, b):

        #print("On Update Catalog --- click")

        detectid = np.int64(self.detectbox.value)
        #print("On Update Catalog --- click -- ID: ",detectid)
        is_continuum = self.detectbox.value[2] == '9'
        hdr = f"hdr{self.detectbox.value[0]}" #20, 21, 30, 40, 50 but really 2.1 is not used

        #print("On Update Catalog:",detectid,is_continuum)

        try:

#            self.read_source_catalog_basic(detectid)
            self.read_source_catalog_update(detectid)

            if self.shotid is None:
                if self.detections_interface_lines is None or self.detections_interface_lines.survey != hdr:
                    # print("On Update Catalog --- get DI")
                    self.detections_interface_lines = Detections(survey=hdr, catalog_type='lines', searchable=False)
                    self.detections_interface_cont = Detections(survey=hdr, catalog_type='continuum', searchable=False)
                    # just assume this works, if not the try will trigger and report

                # print("calling into DI ...")
                if is_continuum:
                    detinfo = self.detections_interface_cont.get_detection_info(detectid)
                else:
                    detinfo = self.detections_interface_lines.get_detection_info(detectid)

                # print("return from DI ...")
                self.shotid = detinfo['shotid'][0]


            self.line_id_drop.value, self.wave_box.value = self.get_line_match(self.z_box.value, self.sourcecat_obswave)
            #print(f"testing: ",self.line_id_drop.value, self.wave_box.value,self.z_box.value, self.sourcecat_obswave)

            self.updateCatalog.disabled=True

            if is_continuum:
                self.wave_box.disabled = True
                self.line_id_drop.disabled = True
            else:
                self.wave_box.disabled = False
                self.line_id_drop.disabled = False

            if b is not None:
                with self.updatecat_box:
                    display(widgets.HBox(
                        [self.real_fake_drop,
                         self.class_labels_drop,
                         self.clusterid_box,
                         self.catalog_sourceid_box,
                         ]))

                    display(widgets.HBox(
                        [self.line_id_drop,
                         self.wave_box,
                         self.z_box,
                         ]))

                    display(widgets.HBox(
                        [ self.catalog_comment_box,
                         ]))

                    display(widgets.HBox(
                        [ self.catalog_status_box,
                          self.doUpdateCatalog,
                         ]))

                #self.focus()
        except Exception as e:
            print("Exception in on_update_catalog()",e)
            print(str(e) + "\n" + traceback.format_exc())


    def on_do_update_catalog(self,b):
        #todo: perform the update to the manual catalog and audit file
        if self.doUpdateCatalog.description == "Reload":
            self.doUpdateCatalog.description = "Update"
            self.on_update_catalog(None)
        else:
            self.record_source_catalog_update()

    def get_mini_button_click(self, b):
        detectid = np.int64(self.detectbox.value)
        try:
            with self.bottombox:
                display(
                    Image(self.elixer_conn_mgr.fetch_image(
                        detectid, report_type="mini"))
                )
        except Exception as e:
            self.status_box.value = str(e) + "\n" + traceback.format_exc()
            print("mini not found")

    def get_hdr_name_for_detectid(self):
        hdr_name = LATEST_HDR_NAME
        try:
            hdr_prefix = int(np.int64(self.detectbox.value) / 1e8)
            hdr_name = HDR_NAME_DICT[hdr_prefix]
        except:
            print("Could not identify or map detectid to hdr version")

        return hdr_name

    def plot_spec(self, matchnum):

        global current_wavelength

        match = self.CatalogMatch["match_num"] == matchnum

        if matchnum > 0:
            col_name = ["blue", "red", "green"]

            object_label = col_name[matchnum - 1] + " Counterpart"

            try:
                coords = SkyCoord(
                    ra=self.CatalogMatch["cat_ra"][match] * u.deg,
                    dec=self.CatalogMatch["cat_dec"][match] * u.deg,
                    frame="icrs",
                )
            except:
                coords = SkyCoord(
                    ra=self.CatalogMatch["ra"][match] * u.deg,
                    dec=self.CatalogMatch["dec"][match] * u.deg,
                    frame="icrs",
                )
        else:
            object_label = (
                "RA="
                + str(self.e_manual_ra.value).zfill(3)
                + " DEC="
                + str(self.e_manual_dec.value).zfill(3)
            )
            coords = SkyCoord(
                ra=self.e_manual_ra.value * u.deg,
                dec=self.e_manual_dec.value * u.deg,
                frame="icrs",
            )

        spec_table = get_spectra(
            coords, ID=np.int64(self.detectbox.value), survey=self.get_hdr_name_for_detectid()
        )

        # if current_wavelength < 0:
        # current_wavelength = self.get_observed_wavelength()

        with self.bottombox:
            for row in spec_table:
                plt.figure(figsize=(15, 2))
                plt.plot(row["wavelength"], row["spec"])
                plt.title(object_label + "        SHOTID = " + str(row["shotid"]))
                plt.xlabel("wavelength (A)")
                plt.ylabel("spec")
            plt.axvline(x=self.get_observed_wavelength(), color="r", linestyle="--")

    def e_blue_button_click(self, b):
        self.plot_spec(1)

    def e_red_button_click(self, b):
        self.plot_spec(2)

    def e_green_button_click(self, b):
        self.plot_spec(3)

    def e_manual_button_click(self, b):
        self.plot_spec(0)

    def get_spec(self, b):
        pass

    def get_det_info(self):
        global CONFIG_HDR2, CONFIG_HDR3, OPEN_DET_FILE, DET_HANDLE

        detid = np.int64(self.detectbox.value)
        
        if (detid >= 2100000000) & (detid < 2190000000):
            self.det_file = CONFIG_HDR2.detecth5
        elif (detid >= 2100000000) & (detid < 2190000000):
            self.det_file = CONFIG_HDR2.contsourceh5
        elif (detid >= 3000000000) & (detid < 3090000000):
            self.det_file = CONFIG_HDR3.detecth5
        elif (detid >= 3090000000) & (detid < 3100000000):
            self.det_file = CONFIG_HDR3.contsourceh5
        elif (detid >= 4000000000) & (detid < 4090000000):
            self.det_file = CONFIG_HDR4.detecth5
        elif (detid >= 4090000000) & (detid < 4100000000):
            self.det_file = CONFIG_HDR4.contsourceh5
        elif (detid >= 5000000000) & (detid < 5090000000):
            self.det_file = CONFIG_HDR5.detecth5
        elif (detid >= 5090000000) & (detid < 5100000000):
            self.det_file = CONFIG_HDR5.contsourceh5

        if OPEN_DET_FILE is None:
            OPEN_DET_FILE = self.det_file
            DET_HANDLE = tables.open_file(self.det_file, 'r')
        else:
            if self.det_file == OPEN_DET_FILE:
                pass
            else:
                DET_HANDLE.close()
                OPEN_DET_FILE = self.det_file
                try:
                    DET_HANDLE = tables.open_file(self.det_file, 'r')
                except Exception:
                    with self.status_box:
                        print("Could not open {}".format(self.det_file))

        if DET_HANDLE is not None:
            detid = np.int64(self.detectbox.value)
            try:
                self.det_row = Table(
                    DET_HANDLE.root.Detections.read_where(
                        "detectid == detid"
                    )
                )
            except Exception as e:
                self.status_box.value = str(e) + "\n" + traceback.format_exc()
                            
    def det_table_button_click(self, b):

        self.get_det_info()
        
        try:
            with self.bottombox:
                display(self.det_row.show_in_notebook())
        except Exception as e:
            self.status_box.value = str(e) + "\n" + traceback.format_exc()        

    def get_sdss_button_click(self, b):
        self.get_det_info()
                    
        url = 'http://skyserver.sdss.org/dr16/en/tools/chart/navi.aspx?ra={:6.4f}&dec={:6.4f}&scale=0.05&opt=G'.format(self.det_row['ra'][0], self.det_row['dec'][0])
        with self.bottombox:
            display(Javascript('window.open("{url}");'.format(url=url)))

    def get_legacy_button_click(self, b):

        self.get_det_info()

        url = 'https://www.legacysurvey.org/viewer?ra={:6.4f}&dec={:6.4f}&layer=ls-dr9&zoom=16'.format(self.det_row['ra'][0], self.det_row['dec'][0])
        with self.bottombox:
            display(Javascript(f'window.open("{url}");'.format(url=url))) 




    def read_source_catalog_basic(self,detectid=None):
        """

        Returns
        -------

        """
        try:

            check_elixer_h5 = False

            if self.source_catalog_h5 is not None:

                #print("read_source_catalog_basic",detectid,self.detectbox.value)

                # if detectid is not None:
                #     q_detectid=np.int64(detectid)
                # else:
                q_detectid=np.int64(self.detectbox.value)

                #print("read_source_catalog_basic q_detectid", q_detectid)

                rows = self.source_catalog_h5.root.SourceCatalog.read_where("detectid==q_detectid")
                ct = len(rows)
                if len(rows)==1:
                    #exactly what we want
                    #print("read rows:", ct)

                    self.sourcecat_z_hetdex = rows['z_hetdex'][0]
                    self.sourcecat_z_hetdex_src = rows['z_hetdex_src'][0].decode()
                    self.sourcecat_z_hetdex_conf = rows['z_hetdex_conf'][0]
                    self.shotid = rows['shotid'][0]
                    self.ra = rows['ra'][0]
                    self.dec = rows['dec'][0]
                    self.sourcecat_obswave = rows['wave'][0]
                    #print(f"here: read_source_catalog_basic", rows['wave'])
                    self.sourcecat_source_id = rows['source_id'][0]
                    self.sourcecat_class_labels = rows['classification_labels'][0]
                    selected = rows['selected_det'][0] #boolean
                    #or can use
                    #selected_flag = rows['flag_seldet'][0]



                    q_source_id = self.sourcecat_source_id
                    if not selected:
                        # print("For now parent checking is TURNED OFF!!!")
                        #
                        #
                        # self.z_hetdex_box.value = f"{self.sourcecat_z_hetdex_src}: {self.sourcecat_z_hetdex:0.4f} " \
                        #                           f"(other source parent)"

                        try:
                            src_rows = self.source_catalog_h5.root.SourceCatalog.read_where(
                                "(source_id==q_source_id) & (selected_det==True)")
                            #must be exactly one
                            if len(src_rows) == 1: #all good
                                self.sourcecat_parent_detectid = src_rows['detectid'][0]

                                self.z_hetdex_box.value = f"{self.sourcecat_z_hetdex_src}: {self.sourcecat_z_hetdex:0.4f} " \
                                                          f"cpid: {self.sourcecat_parent_detectid}"

                            else:
                                self.z_hetdex_box.value = f"{self.sourcecat_z_hetdex_src}: {self.sourcecat_z_hetdex:0.4f} " \
                                                          f"(!catalog error!)"
                        except:
                            print(traceback.format_exc())
                    else:
                        self.sourcecat_parent_detectid = q_detectid

                        self.z_hetdex_box.value = f"{self.sourcecat_z_hetdex_src}: {self.sourcecat_z_hetdex:0.4f} " \
                                                  f"(q:{self.sourcecat_z_hetdex_conf:0.2f})"

                else:
                    #print("fail read rows:", ct)
                    self.status_box.value = f"Source Catalog failure. {ct} entries for {q_detectid}"
                    self.sourcecat_z_hetdex = None
                    self.sourcecat_z_hetdex_src = None
                    self.sourcecat_z_hetdex_conf = None
                    self.sourcecat_class_labels = ""
                    self.sourcecat_parent_detectid = 0
                    self.z_hetdex_box.value = "No Source Catalog Entry"
                    self.shotid = None
                    self.ra = None
                    self.dec = None
                    check_elixer_h5 = True

            elif self.ELIXER_H5 is not None:
                check_elixer_h5 = True

            else:
                check_elixer_h5 = True
                self.status_box.value = f"Source Catalog failure. File not found."
                self.sourcecat_z_hetdex = None
                self.sourcecat_z_hetdex_src = None
                self.sourcecat_z_hetdex_conf = None
                self.z_hetdex_box.value = "Source Catalog Unavailable"
                self.sourcecat_class_labels = ""
                self.sourcecat_parent_detectid = 0
                self.shotid=None
                self.ra = None
                self.dec = None

            if check_elixer_h5:

                if self.ELIXER_H5 is None and self.sourcecat_root_path is not None:
                    #print("read_source_catalog_update()    --- 2")
                    try:
                        elixer_h5fn = os.path.join(self.sourcecat_root_path,
                                                           HETDEX_ELIXER_SUBPATH,HETDEX_ELIXER_FILE)
                        self.ELIXER_H5 = tables.open_file(elixer_h5fn)
                    except:
                        #print("read_source_catalog_update()    --- 2b")
                        print(traceback.format_exc())
                        pass

                if self.ELIXER_H5 is not None:
                    q_detectid = detectid  # self.detid
                    rows = self.ELIXER_H5.root.Detections.read_where("detectid==q_detectid")
                    ct = len(rows)
                    if len(rows) == 1:  # there is an elixer entry
                        self.catalog_status_box.value += f"\nELiXer catalog success. {ct} entries for {detectid}"

                        self.sourcecat_z_hetdex = rows['z_best'][0]
                        self.sourcecat_z_hetdex_src = "elixer.h5"
                        self.sourcecat_z_hetdex_conf = rows['z_best_pz'][0]
                        self.sourcecat_obswave = rows['wavelength_obs'][0]
                        self.sourcecat_source_id = 0
                        self.sourcecat_class_labels = rows['classification_labels'][0]
                        self.shotid = rows['shotid'][0]
                        self.ra = rows['ra'][0]
                        self.dec = rows['dec'][0]
                        self.sourcecat_parent_detectid = rows['cluster_parent'][0]
                        # self.sourcecat_sourceid = rows['sourceid'] #not in elixer

                       # print("elixer z value", self.sourcecat_z_hetdex)
                    else:
                        self.status_box.value += f"\nSource Catalog fail. ELiXer catalog fail. No/invalid entries for {detectid}"
                else:
                    self.status_box.value += f"\nSource Catalog fail. ELiXer catalog fail. No entries for {detectid}"


        except Exception as e:
            #print("Exception!!!", e)
            self.status_box.value = str(e) + "\n" + traceback.format_exc()


    def read_source_catalog_update(self, detectid=None, refresh=True):
        """
        get data to display in the source catalog update dialog

        use the official source catalog with any updates from the update table

        Returns
        -------

        """

        try:

            #first get the main source catalog info
            #this is the new h5 file

            #print("read_source_catalog_update()")

            #print("read_source_catalog_update: detectid: ", detectid)

            if detectid is None:
                detectid = np.int64(self.detectbox.value)

            #print("read_source_catalog_update: detectid: ", detectid)

            self.catalog_status_box.value = ""

            #read_source_catalog_basic() should have already happened
            if self.sourcecat_z_hetdex is None: #no entry was found, try elixer??
                #print("read_source_catalog_update()    --- 1")
                self.catalog_status_box.value = f"Source Catalog fail. No entries for {detectid}"
                #try the elixer catalog, if there is one

                if self.ELIXER_H5 is None and self.sourcecat_root_path is not None:
                    #print("read_source_catalog_update()    --- 2")
                    try:
                        elixer_h5fn = os.path.join(self.sourcecat_root_path,
                                                           HETDEX_ELIXER_SUBPATH,HETDEX_ELIXER_FILE)
                        self.ELIXER_H5 = tables.open_file(elixer_h5fn)
                    except:
                        #print("read_source_catalog_update()    --- 2b")
                        print(traceback.format_exc())
                        pass


                if self.ELIXER_H5 is not None:
                    #print("read_source_catalog_update()    --- 3")
                    q_detectid=detectid #self.detid
                    rows = self.ELIXER_H5.root.Detections.read_where("detectid==q_detectid")
                    ct = len(rows)
                    if len(rows) == 1: #there is an elixer entry
                        print("read_source_catalog_update()    --- 4")
                        self.catalog_status_box.value += f"\nELiXer catalog success. {ct} entries for {q_detectid}"
                        self.sourcecat_z_hetdex = rows['z_best'][0]
                        self.sourcecat_z_hetdex_src = "elixer.h5"
                        self.sourcecat_z_hetdex_conf = rows['z_best_pz'][0]
                        self.sourcecat_obswave = rows['wavelength_obs'][0]
                        self.sourcecat_source_id = 0
                        self.sourcecat_class_labels = rows['classification_labels'][0]
                        self.shotid = rows['shotid'][0]
                        self.ra = rows['ra'][0]
                        self.dec = rows['dec'][0]
                        self.sourcecat_parent_detectid = rows['cluster_parent'][0]
                        #self.sourcecat_sourceid = rows['sourceid'] #not in elixer

                        #print("elixer z value", self.sourcecat_z_hetdex)
                    else:
                        #print("elixer entry not found.")
                        self.catalog_status_box.value += f"\nELiXer catalog fail. No entries for {q_detectid}"
            else:
                #print("testing. sourcecat_z has value ", self.sourcecat_z_hetdex)
                pass

            if self.sourcecat_parent_detectid is not None:
                self.clusterid_box.value = str(self.sourcecat_parent_detectid)
            else:
                self.clusterid_box.value = "0"


            if self.sourcecat_source_id is not None:
                self.catalog_sourceid_box.value = str(self.sourcecat_source_id)
            else:
                self.catalog_sourceid_box.value = "N/A"

            #remember, can't just keep this around, since we only lock when we want to update
            #and the table can become stale if someone else edits


            sct_row,*_ = self.sct.get_row(detectid,refresh=refresh) #don't need to refresh here (yet), don't care about the index

            if self.sct.status != 0:
                self.status_box.value = self.sct.status_msg
                self.catalog_status_box.value = self.sct.status_msg
            elif sct_row is None: #there is no such row, so no existing update
                self.catalog_status_box.value = "No prior updates"
            else: #the sct_row is good, load the data
                self.catalog_status_box.value = f"Last update: ({sct_row['revision']}) {sct_row['revision_user']} " \
                                                f"at {sct_row['revision_date']} UTC"


                #print(f"Testing: classification_labels: ", sct_row['classification_labels'], type(sct_row['classification_labels']))

                if sct_row['classification_labels'] is not None and not np.ma.is_masked(sct_row['classification_labels']) and len(str(sct_row['classification_labels'])) > 0:
                    try:
                        self.sourcecat_class_labels = sct_row['classification_labels'].decode()
                    except:
                        self.sourcecat_class_labels = sct_row['classification_labels']

                if sct_row['z'] is not None and sct_row['z'] > -1:
                    self.sourcecat_z_hetdex = sct_row['z']
                if sct_row['parent_detectid'] is not None and sct_row['parent_detectid'] > -1:
                    self.sourcecat_parent_detectid = sct_row['parent_detectid']
                if sct_row['conf_status'] is not None and sct_row['conf_status'] != 0:
                    #this is a short list, but this is otherwise a dumb brute force
                    #assume (and this should be true) that the keys and values are each unique
                    for key,val in zip(real_fake_dict.keys(),real_fake_dict.values()):
                        if val == sct_row['conf_status']:
                            self.real_fake_drop.value = key
                            break
                if sct_row['comments'] is not None and not np.ma.is_masked(sct_row['comments']) and len(sct_row['comments']) > 0:
                    self.catalog_comment_box.value = sct_row['comments']
                self.source_catalog_detid_revision = sct_row['revision']


            line_id, rest_wave = self.get_line_match(self.sourcecat_z_hetdex,obswave=self.sourcecat_obswave)

            #print(f"Testing:",rest_wave, line_id, self.sourcecat_z_hetdex)

            self.wave_box.value = rest_wave
            self.line_id_drop.value = line_id
            if self.sourcecat_class_labels is not None:
                if str(self.sourcecat_class_labels) == '--':
                    self.class_labels_drop.value = None
                else:
                    try:
                        self.class_labels_drop.value = self.sourcecat_class_labels.decode()
                    except:
                        self.class_labels_drop.value = self.sourcecat_class_labels

            if self.sourcecat_z_hetdex is not None:
                self.z_box.value = self.sourcecat_z_hetdex



        except Exception as e:
            self.catalog_status_box.value = str(e) + "\n" + traceback.format_exc()

    def record_source_catalog_update(self):
        """

        Returns
        -------

        """

        try:

            #get a dictionary, update the fields then apply the update
            row = self.sct.get_row_dict()


            row['detectid'] = np.int64(self.detectbox.value)
            row['ra'] = self.ra
            row['dec'] = self.dec
            row['shotid'] = self.shotid
            row['obs_wave'] = self.sourcecat_obswave
            row['revision'] = self.source_catalog_detid_revision
            #note: revision_date, revision_user set by the source_catalog updater

            #for real/fake prohibit recording of accidental selection of the separator (value = -999)
            #since all proper values are non-zero and positive, if an invalid (negative) value is there, set to 0 (the default)
            row['conf_status'] = max(real_fake_dict[real_fake_default],real_fake_dict[self.real_fake_drop.value])
            row['z'] = self.z_box.value
            row['parent_detectid'] = self.clusterid_box.value
            if self.class_labels_drop.value is None or len(self.class_labels_drop.value) == 0:
                row['classification_labels'] = ""
            else:
                row['classification_labels'] = self.class_labels_drop.value
            row['comments'] = self.catalog_comment_box.value


            #print(f"Testing: update row: ", row)
            self.sct.update_table(row)

            if self.sct.status == 0:

                #reload to update the revision
                self.on_update_catalog(None)

                self.catalog_status_box.value = f"Success."
                try:
                    sct_row, *_ = self.sct.get_row(row['detectid'], refresh=True)  # don't need to refresh here (yet), don't care about the index

                    self.catalog_status_box.value += f"\nLast update: ({sct_row['revision']}) {sct_row['revision_user']} " \
                                                    f"at {sct_row['revision_date']} UTC"
                except:
                    pass
            else:
                self.catalog_status_box.value = self.sct.status_msg
                if "stale" in self.sct.status_msg:
                    self.doUpdateCatalog.description ="Reload"


        except Exception as e:
            self.catalog_status_box.value = "Update Failed.\n"
            self.catalog_status_box.value += str(e) + "\n" + traceback.format_exc()





