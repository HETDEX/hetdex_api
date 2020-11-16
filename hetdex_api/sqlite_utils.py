"""
Utility functions for SQLite3 ELiXer Report imaging databases

Each HETDEX Data release will contain multiple serverless ELiXer report imaging databsaes with approx 100,000 images
  per database file with each file containing one type of report image.

The database files should be named as: elixer_reports_<detectid prefix>_<option report type>.db
    For example: "elixer_reports_20003.db" contains HDR2 (from the leading '2') main report images that all begin with 20003
                        ... so, images 2000300000.png to 2000399999.png

                "elixer_reports_10006_nei.db" contain the Neighborhood report images for HDR1 that fit the  10006*nei.png pattern

                "elixer_reports_20000_mini.db" contain the mini report images (e.g. for smart phone app) for HDR2, with
                        detectIDs 2000000000 to 2000099999
"""


import sqlite3
from sqlite3 import Error
import numpy as np

import os.path as op
FILENAME_PREFIX = "elixer_reports_" #keep the trailing underscore
REPORT_TYPES = ["report","nei","mini"]

try:
    from hetdex_api.config import HDRconfig
except:
    # print("Warning! Cannot find or import HDRconfig from hetdex_api!!")
    pass

#key is the HDR version number, value is list of directories that contain ELiXer imaging databses
#Base paths
#  xx0 = standard hetdex
#  xx6 = broad emission lines (still in with the xx0 detections as of hdr2.1)
#  xx9 = continuum sources
DICT_DB_PATHS = {10: ["/work/03946/hetdex/hdr1/detect/image_db"
                     ],
                 20: ["/scratch/03261/polonius/hdr2/detect/image_db",
                      "/work/03261/polonius/hdr2/detect/image_db",
                      "/work/03946/hetdex/hdr2/detect/image_db"
                     ],
                 21: ["/scratch/03946/hetdex/hdr2.1/detect/image_db",
                      "/scratch/03261/polonius/hdr2.1/detect/image_db",
                      "/work/03946/hetdex/hdr2.1/detect/image_db",
                      "/work/03261/polonius/hdr2.1/detect/image_db"
                      ],
                 }
#
# add paths from hetdex_api to search (place in first position)
#
for v in DICT_DB_PATHS.keys():
    try:
        release_number = v/10.0
        if v % 10 == 0:
            release_string = "hdr{:d}".format(int(release_number))
        else:
            release_string = "hdr{:2.1f}".format(release_number)

        DICT_DB_PATHS[v].insert(0,op.join(HDRconfig(survey=release_string).elix_dir))
    except:# Exception as e:
        #print(e)
        continue


def get_elixer_report_db_path(detectid,report_type="report"):
    """
    Return the top (first found) path to database file based on the detectid (assumes the HDR version is part of the
    prefix, i.e. HDR1 files are 1000*, HDR2 are 2000*, and so on)
    :param detectid:
    :param report_type: choose one of "report" (normal ELiXer report image) [default]
                               "nei" (ELiXer neighborhood image)
                               "mini" (ELiXer mini-report image for phone app)
    :return: None or database filename
    """

    detect_prefix = None
    db_path = None
    try:
        detect_prefix = int(np.int64(detectid) / 1e5)
        hdr_prefix = int(np.int64(detectid)/1e8)

        #keep the leading underscore
        if report_type == "report":
            ext = ""
        elif report_type == "nei":
            ext = "_nei"
        elif report_type == "mini":
            ext = "_mini"
        else: #assume same as report
            ext = ""

        if detect_prefix is not None:
            if hdr_prefix in DICT_DB_PATHS.keys():
                paths = DICT_DB_PATHS[hdr_prefix]
                for p in paths:
                    if op.exists(p):
                        fqfn = op.join(p, FILENAME_PREFIX + str(detect_prefix) + ext + ".db")
                        if op.isfile(fqfn):
                            db_path = fqfn
                        break
            else:
                #print("Invalid HDR version")
                return None
    except Error as e:
        print(e)

    return db_path



def get_db_connection(fn,readonly=True):
    """
    return a SQLite3 databse connection object for the provide databse filename

    assumes file exists (will trap exception if not and return None)
    :param fn:
    :return: None or connection object
    """

    conn = None
    try:
        if fn is not None:
            if readonly:
                conn = sqlite3.connect("file:" +fn + "?mode=ro",uri=True)
            else:
                conn = sqlite3.connect(fn)
    except Error as e:
        print(e)

    return conn


def fetch_elixer_report_image(conn,detectid):
    """
    Return a single image (image type (png, jpg) and report type (report, neighborhood, mini) depend on the database connection

    :param conn: a sqlite3.Connection object or a path to a database
    :param detectid: HETDEX detectid (int64 or string)
    :return: None or single image
    """

    try:
        keep_conn_open = True
        if type(conn) != sqlite3.Connection:
            #could be a file
            if op.isfile(conn):
                conn = get_db_connection(conn,readonly=True)

                if type(conn) != sqlite3.Connection:
                    print("Invalid databse connection.")
                    return None
                keep_conn_open = False
            else:
                print("Invalid database connection.")
                return None

        cursor = conn.cursor()
        sql_read_blob = """SELECT image from report where detectid = ?"""

        cursor.execute(sql_read_blob, (str(detectid),))
        image = cursor.fetchall()
        #get back a list of tuples (each list entry is one row, the tuple is the row, so
        #we want the 1st entry (detectid is unique, so should be one or none, and the second column which is the image)
        cursor.close()
        if not keep_conn_open:
            conn.close()

        if image is not None:
            if len(image) == 1:
                return image[0][0]
            elif len(image) == 0:
                print("No matching detectid found")
                return None
            else:
                print("Unexpected number of images returned")
                return None
        else:
            print("None returned from image fetch.")
            return None
    except Error as e:
        print(e)

    return None


def build_elixer_report_image_db(db_name,img_dir,img_regex):
    """
    Not for general use. Should normally be called once per image type per grouping of (100,000) per
    data release.

    If db already exists, this will insert new images and replace existing images (if detectIDs match)

    This DOES NOT validate the images or check that they are appropriate for the db_name.

    Progress is reported per 100 inserts

    Insert speed does not depend on database size, but does depend on disk speed (local vs network, SSD vs HD, etc)

    :param db_name: name (with path) of the SQLite output db
    :param img_dir: location of all images to be imported  (is NOT recursive)
    :param img_regex: generally simple wildcard, like "*.png" or "*nei.png" or "*.jpg", etc
    :return:
    """

    import glob
    import re
    import time
    # import gzip #only get about 5% compression on pngs at much higher cpu cost
    import datetime

    def build_schema(conn):
        # build up sql commands to issue
        try:
            sql_create_report_image_table = """ CREATE TABLE IF NOT EXISTS report (
                                                    detectid BIGINT PRIMARY KEY,
                                                    image BLOB NOT NULL
                                                ); """

            # create report table
            cursor = conn.cursor()
            cursor.execute(sql_create_report_image_table)
            cursor.close()
            conn.commit()


            #create index: not necessary; is autocreated with BIGING Primary Key
            # sql_create_index = """ CREATE UNIQUE INDEX idx_detectid ON report (detectid); """
            # cursor = conn.cursor(sql_create_index)
            # cursor.execute()
            # cursor.close()
            # conn.commit()

            return True
        except Exception as e:
            print(e)
            return False

    def read_image(fn):
        blob_data = None
        detectid = None
        try: #assume name like <detectid><optional chars>.png
             detectid = int(re.findall('\d+',op.basename(fn).split(".")[0])[0])
             with open(fn, 'rb') as file:
                 blob_data = file.read()
        except Exception as e:
            print("Exception in read_image (bad detectid):", fn, e)

        return blob_data, detectid

    def insert_image(conn, detectid, data):
        if (detectid is None) or (data is None):
            return

        try:
            cursor = conn.cursor()
            sql_insert_blob = """ INSERT INTO report
                                              (detectid, image) VALUES (?, ?)"""
            try:
                cursor.execute(sql_insert_blob, (detectid, data))
                conn.commit()
            except Exception as e:
                if type(e) == sqlite3.IntegrityError:
                    # could be image already exists, in which case overwrite
                    try:
                        sql_update_blob = """ UPDATE report
                                              SET image = ? WHERE detectid = ?"""
                        cursor.execute(sql_update_blob, (data, detectid))
                        conn.commit()

                    except Exception as e:
                        print("second exception in insert_image:", detectid, e)
                else:
                    print("exception in insert_image:", detectid, e)

            cursor.close()
        except Exception as e:
            print(e)
            try:
                cursor.close()
            except:
                pass

    def import_images_from_path(conn,img_dir,img_regex):
        ct = 0
        total_inserts = 0
        estimated_total = 0
        modulo = 100 #print a status statement every modulo inserts

        filelist = glob.glob(op.join(img_dir, img_regex))  # just the reports, not the neighbors or mini
        estimated_total = len(filelist)

        print(f"Inserting {estimated_total} images ... ")

        if estimated_total < 1:
            return #nothing to do

        start_time = int(round(time.time() * 1000))

        for f in filelist:
            try:
                blob, detectid = read_image(f)
                insert_image(conn, detectid, blob)
                ct += 1
            except Exception as e:
                print("exception in import_images_from_path:", e)

            if ct >= modulo:
                try:
                    time_diff = int(round(time.time() * 1000)) - start_time
                    total_inserts += ct
                    per_insert_time = time_diff / ct / 1000.
                    print(f"{db_name}: Inserted {ct} ({per_insert_time:#0.3f}s per insert). Total {total_inserts}/{estimated_total} "
                          f"({float(total_inserts / estimated_total) * 100.:#0.1f}%). "
                          f"Remaining time ({datetime.timedelta(seconds=round(per_insert_time * (estimated_total - total_inserts)))}) ...")
                    start_time = int(round(time.time() * 1000))  # reset the timer
                    ct = 0
                except:
                    print("print progress failed....")

        #final log (remainder after last block of inserts
        try:
            time_diff = int(round(time.time() * 1000)) - start_time
            total_inserts += ct
            per_insert_time = time_diff / ct / 1000.
            print(f"Inserted {ct} ({per_insert_time:#0.3f}s per insert). Total {total_inserts}/{estimated_total} "
                  f"({float(total_inserts / estimated_total) * 100.:#0.1f}%). ")
        except:
            print("print progress failed....")

    #
    # main execution part
    #
    try:

        if not op.isfile(db_name):
            conn = sqlite3.connect(db_name) #not read only
            if type(conn) != sqlite3.Connection:
                print("Failed to create db connection")
                return False
            elif not build_schema(conn):
                print("Failed to build schema")
                return False
        else:
            conn = get_db_connection(db_name,readonly=False)
            if type(conn) != sqlite3.Connection:
                print("Failed to create db connection")
                return False

        import_images_from_path(conn,img_dir,img_regex)

    except Exception as e:
        print(e)

# 20200529 DD
# will revisit later if this idea becomes useful
#
# def build_local_imaging_dbs():
#     """
#     super simplified for testing
#     build the imaging databases using default settings in the current directory
#     :return:
#     """
#     import glob
#     #get the list of images and report types
#     rpt_min, rpt_max = None, None
#     rpt_list = sorted(glob.glob("*[0-9].png"))
#     if (rpt_list is not None) and (len(rpt_list) > 0):
#         rpt_min = int(rpt_list[0].rstrip(".png"))
#         rpt_max = int(rpt_list[-1].rstrip(".png"))
#
#     nei_min, nei_max = None, None
#     nei_list = sorted(glob.glob("*[0-9]*nei.png"))
#     if (nei_list is not None) and (len(nei_list) > 0):
#         nei_min = rpt_list[0].rstrip("nei.png")
#         nei_max = rpt_list[-1].rstrip("nei.png")
#         nei_min = int(nei_min.replace("_","")) #might not have an "_"
#         nei_max = int(nei_max.replace("_", ""))
#
#     mini_min, mini_max = None, None
#     mini_list = sorted(glob.glob("*[0-9]*mini.png"))
#     if (mini_list is not None) and (len(mini_list) > 0):
#         mini_min = rpt_list[0].rstrip("mini.png")
#         mini_max = rpt_list[-1].rstrip("mini.png")
#         mini_min = int(mini_min.replace("_","")) #might not have an "_"
#         mini_max = int(mini_max.replace("_", ""))
#
#     #organize by
#


class ConnMgr():
    """
    Primitive container for managing SQLite connection (to avoid repeated path search and connection building)
    Just for reading
    """

    def __init__(self):
        self.conn_dict = {} #key = detectid_prefix + type (i.e. "10003" or "10004nei" or "20007mini"

    def __del__(self):
        self.close_conns()

    def get_connection(self,detectid,report_type="report"):
        conn = None
        try:
            if not report_type in REPORT_TYPES:
                return None

            detect_prefix = int(np.int64(detectid) / 1e5)
            dkey = str(detect_prefix)+report_type
            if dkey in self.conn_dict.keys():
                conn = self.conn_dict[dkey]
            else:
                try:
                    #all ConnMgr connections are read-only (uri=True)
                    conn = get_db_connection(get_elixer_report_db_path(detectid,report_type),readonly=True)
                    if type(conn) != sqlite3.Connection:
                        conn = None
                    else:
                        self.conn_dict[dkey] = conn
                except Exception as e:
                    print(e)
        except Exception as e:
            print(e)

        return conn

    def close_conns(self):
        for key in self.conn_dict.keys():
            try:
                self.conn_dict[key].close()
            except:
                pass

        self.conn_dict.clear()


    def fetch_image(self,detectid,report_type="report"):
        """
        wrapper just to make code cleaner

        :param detectid:
        :param report_type:
        :return:
        """

        img = None
        try:
            conn = self.get_connection(detectid,report_type)
            if type(conn) == sqlite3.Connection:
                img = fetch_elixer_report_image(conn,detectid)

        except Exception as e:
            print(e)
            raise

        return img
