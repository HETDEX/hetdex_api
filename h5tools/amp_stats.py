"""

Initial Creation 2024-09-02

@author: Dustin Davis

Collection of code to support the determination of which amps, IFUs, shots are "bad"

Basic statistics with an amp+exposure as the functional unit with aggregation to the shot level.


"""

import tables
import os.path as op
import numpy as np
import glob

import astropy.stats.biweight as biweight
from astropy.table import Table, Row, join
#from astropy.io import ascii
import copy
import pickle
from hetdex_api.config import HDRconfig

import traceback



from tqdm import tqdm

#########################################
# helpers
#########################################

AMP_LIST = ['LL', 'LU', 'RL', 'RU']

def flx2cts(calfib, calfib_counts, calfib_ffsky):
    # return calfib_ffsky as counts

    zero = calfib_counts == 0
    safe = copy.copy(calfib)
    safe[calfib == 0] = np.nan

    out = calfib_ffsky * calfib_counts / safe
    out[zero] = 0
    return out


def datevshot2int(datevshot):
    try:
        return np.int64(datevshot.replace("v", ''))
    except:
        return datevshot


def shotid2datevshot(shotid):
    try:
        return str(shotid)[0:8] + "v" + str(shotid)[8:]
    except:
        return shotid


def stats_get_shot_h5(config=None, shotid=None, fqfn=None, append=False):
    """
    Find, open and return handle to a shot h5 file for use with the other stats calls

    This is almost the equivalent of hetdex_api.shot::open_shot_file(shotid, survey=hdr_version), but this
    also allows the user to specify a full path for a particular shot file without using the survey

    Shot file MUST exist
    Can use either the fully qualifed path + filename OR the shotid and an HDRConfig object

    Allow Exception to percolate up so caller can see the problem,
        if file does not exist as requested or cannot open with requested access

    Input: 1. filepath and fileanme (for direct open)
           2. -OR- HDRConfig and shotid

    """

    # validate
    if append:
        mode = "r+"
    else:
        mode = "r"

    if fqfn is None:
        # use config and shotid
        fn = str(shotid)  # whether as in int or string
        if fn[8].lower() != 'v':  # either r, s or just an integer
            if fn[8].isdigit():
                fn = fn[0:8] + "v" + fn[8:]  # was an integer so, add in the 'v'
            elif fn[8].lower() in ['r', 's']:  # assume a datevshot style format
                fn = list(fn)
                fn[8] = 'v'  # make explicitly datevshot
                fn = ''.join(fn)
            else:  # might not be properly formatted, let the exception pass up if there is a problem
                pass

        fn += ".h5"

        fqfn = op.join(config.data_dir, fn)

    h5 = tables.open_file(fqfn, mode=mode)

    return h5


def stats_shot_dict_to_table(shot_dict):
    """

    Parameters
    ----------
    shot_dict - single dictionary or an array or list of dictionaries

    Returns
    -------
    astropy table version of the dictionary(ies)

    """


    try:

        # tun nested dictionaries into lists at the amp+exposure level
        if isinstance(shot_dict, dict):  # assume a single entry
            shot_dict_array = [shot_dict]
        else:
            shot_dict_array = shot_dict

        shotid = []  # shot_dict["shotid"]
        mf_list = []
        exp_list = []
        im_median = []
        maskfraction = []
        chi2fib_avg = []
        frac_c2 = []
        sky_sub_rms = []
        sky_sub_rms_rel = []
        sky_sub_rms_median_exp = []
        dither_relflux = []
        norm = []
        frac_0 = []
        n_lo = []
        #n_lo_ss = []
        N_cont = []


        #original amp.dat
        #kNorm = []
        #kN_c = []
        #kNlo = [] #same as n_lo
        Scale = []
        Avg = []
        Avg_orig = []
        kchi = []
        #Frac_c2 = []
        #Frac0 = []

        for sd in shot_dict_array:

            for ifu in sd["ifu_dict_array"]:
                for amp in ifu["amp_dict_array"]:
                    for exp in amp["exposures"]:
                        shotid.append(sd["shotid"])
                        mf_list.append(exp['multiframe'])
                        exp_list.append(exp['expid'])

                        # Erin Stuff
                        im_median.append(exp['im_median'])
                        maskfraction.append(exp['maskfraction'])
                        chi2fib_avg.append(exp['chi2fib_avg'])
                        frac_c2.append(exp['frac_c2'])

                        sky_sub_rms.append(exp['sky_sub_rms'])
                        sky_sub_rms_rel.append(exp['sky_sub_rms_rel'])
                        try:
                            #not all versions have this in the exp_dict
                            sky_sub_rms_median_exp.append(exp['sky_sub_rms_median']) #the one actually used in the calculation
                        except:
                            sky_sub_rms_median_exp.append(np.nan)
                            #sky_sub_rms_median_exp.append(np.nanmedian(sd['sky_sub_rms_median_exp']))
                        dither_relflux.append(exp['dither_relflux'])
                        norm.append(exp['norm'])
                        frac_0.append(exp['frac_0'])
                        n_lo.append(exp['n_lo'])
                        #n_lo_ss.append(exp['n_lo_ss'])
                        N_cont.append(exp['N_cont'])


                        #tests (like original amp.dat)
                        #kNorm.append(exp['kNorm'])
                        #kN_c.append(exp['kN_c'])
                        #kNlo.append(exp['kNlo']) #same as n_lo
                        Scale.append(exp['Scale'])
                        Avg.append(exp['Avg'])
                        Avg_orig.append(exp['Avg_orig'])
                        kchi.append(exp['kchi'])


        output_tb = Table(
            [
                shotid * np.ones(len(exp_list), dtype=int),
                mf_list,
                np.array(exp_list).astype(np.int32), # !! cannot be int8 as it is necessarily xlated to BOOL by fits i/o
                im_median,
                maskfraction,
                Avg,  # Karl Average
                Scale, #Karl Scale
                chi2fib_avg,
                frac_c2,
                frac_0,
                np.array(n_lo).astype(np.int32),
                Avg_orig,
                sky_sub_rms,
                sky_sub_rms_rel,
                sky_sub_rms_median_exp,
                dither_relflux,
                norm,
                #kNorm ,
                #kN_c ,
                #kNlo, #same as n_lo
                kchi,
                np.array(N_cont).astype(np.int32),
                #n_lo_ss


            ],
            names=["shotid", "multiframe", "expnum", "im_median", "MaskFraction", "Avg",'Scale','chi2fib_med', 'frac_c2', 'frac_0',
                   'n_lo', 'Avg_orig','sky_sub_rms', 'sky_sub_rms_rel','sky_sub_rms_median','dither_relflux','norm',
                   'kchi','N_cont'
                   ]

        )

        return output_tb

    except Exception as e:
        # todo: error handling
        print("stats_shot_dict_to_table ", print(traceback.format_exc()))


###########################################
# File handling
###########################################

def stats_save_as(shot_dict,outfile,format="ascii",overwrite=True,oldstyle=False,header=True):
    """

    Parameters
    ----------
    shot_dict - single dictionary or list or array of dictionries
    outfile - optional output filename (include path if not ".")
    format - astropy table format
    overwrite - True (default) or False
    oldstyle - if True, mimic the old amp.dat style with the columns in the same order (and newer columns follow)
                (note: oldstyle, per amp.dat, DOES NOT HAVE COLUMN HEADERS)

    Returns
    -------
    the astropy table that was written

    """


    default_bad = np.nan
    tab = stats_shot_dict_to_table(shot_dict)
    tab.sort(['shotid'])

    if not oldstyle:
        tab.write(outfile,overwrite=overwrite,format=format)
    else: #try to mimic the older amp.dat stype, no column headings, short precision, different shotformat
        if not overwrite:
            if op.exists(outfile):
                print(f"{outfile} exists and overwrite set to False.")
                return None

        #column order from amp.dat
        #['shotid','multiframe','Factor','N_c', 'Avg', 'Scale', 'W0', 'W1', 'n_lo', 'Avg_orig', 'chi2fib_med',
        # 'frac_c2', 'frac_0'])

        # <TableColumns names=('shotid','multiframe','expnum','im_median','MaskFraction','chi2fib_med','frac_c2',
        # 'frac_0','n_lo','sky_sub_rms','sky_sub_rms_rel','dither_relflux'
        with open(outfile,"w") as f:
            if header:
                # !!! make sure the order matches the write at the end !!!
                try: #note some extra tabs for formatting

                    #less efficient but easier to read
                    f.write(f"#dateshot\t")
                    f.write(f"multiframe\t")
                    f.write(f"Factor\t")
                    f.write(f"N_c\t")
                    f.write(f"Avg\t")
                    f.write(f"Scale\t")
                    f.write(f"W0\t")
                    f.write(f"W1\t")
                    f.write(f"n_lo\t")
                    f.write(f"Avg_orig\t")
                    f.write(f"chi2fib_med\t")
                    f.write(f"frac_c2\t")
                    f.write(f"frac_0\t")
                    f.write(f"im_median\t")
                    f.write(f"MaskFraction\t")
                    f.write(f"sky_sub_rms\t")
                    f.write(f"sky_sub_rms_rel\t")
                    f.write(f"sky_sub_rms_median_exp\t")
                    f.write(f"dither_relflux\t")
                    f.write(f"norm\t")
                    #f.write(f"kN_c\t")
                    f.write(f"kchi\t")
                    f.write(f"N_cont\t")
                    #f.write(f"n_lo_ss\t")
                    f.write("\n")

                except:
                    print("stats_save_as ", print(traceback.format_exc()))



            for row in tab:
                try:
                    dvse = f"d{str(row['shotid'])[0:8]}s{str(row['shotid'])[8:]}exp{str(row['expnum']).zfill(2)}"
                    mf = row['multiframe'][6:]
                except:
                    #if could not build the datevshot string or the multiframe, move on to the next one
                    print("stats_save_as ", print(traceback.format_exc()))
                    continue

                #if the column does not exist, use the default bad value
                try:
                    Factor = f"{0:0.3f}"
                except:
                    Factor = f"{default_bad}"

                try:
                    N_c = f"{0:d}"
                except:
                    N_c = f"{default_bad}"

                try:
                    Avg = f"{0:0.2f}"
                except:
                    Avg = f"{default_bad}"

                try:
                    Scale = f"{0:0.2f}"
                except:
                    Scale = f"{default_bad}"

                try:
                    W0 = f"{0:0.2f}"
                except:
                    W0 = f"{default_bad}"

                try:
                    W1 = f"{0:0.2f}"
                except:
                    W1 = f"{default_bad}"

                try:
                    n_lo = f"{row['n_lo']:d}"
                except:
                    n_lo = f"{default_bad}"

                # try:
                #     n_lo_ss = f"{row['n_lo_ss']:d}"
                # except:
                #     n_lo_ss = f"{default_bad}"

                try:
                    N_cont = f"{row['N_cont']:d}"
                except:
                    N_cont = f"{default_bad}"

                try:
                    Avg_orig = f"{0:0.2f}"
                except:
                    Avg_orig = f"{default_bad}"

                try:
                    chi2fib_med = f"{row['chi2fib_med']:0.2f}"
                except:
                    chi2fib_med = f"{default_bad}"

                try:
                    frac_c2 = f"{row['frac_c2']:0.4f}"
                except:
                    frac_c2 = f"{default_bad}"

                try:
                    frac_0 = f"{row['frac_0']:0.4f}"
                except:
                    frac_0 = f"{default_bad}"


                #newer columns (that do not have an older amp.dat equivalent)



                try:
                    im_median = f"{row['im_median']:0.4f}"
                except:
                    im_median = f"{default_bad}"

                try:
                    MaskFraction = f"{row['MaskFraction']:0.4f}"
                except:
                    MaskFraction = f"{default_bad}"

                try:
                    sky_sub_rms = f"{row['sky_sub_rms']:0.4f}"
                except:
                    sky_sub_rms = f"{default_bad}"

                try:
                    sky_sub_rms_rel = f"{row['sky_sub_rms_rel']:0.4f}"
                except:
                    sky_sub_rms_rel = f"{default_bad}"

                try:
                    sky_sub_rms_median_exp = f"{row['sky_sub_rms_median_exp']:0.4f}"
                except:
                    sky_sub_rms_median_exp = f"{default_bad}"

                try:
                    #this one is really at a shot level, but has to be repeated at each amp
                    dither_relflux = f"{row['dither_relflux']:0.4f}"
                except:
                    dither_relflux = f"{default_bad}"

                try:
                    norm = f"{row['norm']:0.4f}"
                except:
                    norm = f"{default_bad}"



                # testing amp.dat-like columns (that DO have an older amp.dat equivalent)

                # try:
                #     kNorm = f"{row['kNorm']:0.4f}"
                # except:
                #     kNorm = f"{default_bad}"

                # try:
                #     kN_c = f"{row['kN_c']:0.4f}"
                # except:
                #     kN_c = f"{default_bad}"

                try:
                    Avg = f"{row['Avg']:0.4f}"
                except:
                    Avg = f"{default_bad}"

                try:
                    Scale = f"{row['Scale']:0.4f}"
                except:
                    Scale = f"{default_bad}"

                # try:
                #     kW0 = f"{row['kW0']:0.4f}"
                # except:
                #     kW0 = f"{default_bad}"
                #
                # try:
                #     kW1 = f"{row['kW1']:0.4f}"
                # except:
                #     kW1 = f"{default_bad}"

                try:
                    Avg_orig = f"{row['Avg_orig']:0.4f}"
                except:
                    Avg_orig = f"{default_bad}"

                try:
                    kchi = f"{row['kchi']:0.4f}"
                except:
                    kchi = f"{default_bad}"

                # try:
                #     kNlo = f"{row['kNlo']:0.4f}"
                # except:
                #     kNlo = f"{default_bad}"

                try:
                    f.write(f"{dvse}\t")
                    f.write(f"{mf}\t")
                    f.write(f"{Factor}\t")
                    f.write(f"{N_c}\t")
                    f.write(f"{Avg}\t")
                    f.write(f"{Scale}\t")
                    f.write(f"{W0}\t")
                    f.write(f"{W1}\t")
                    f.write(f"{n_lo}\t")
                    f.write(f"{Avg_orig}\t")
                    f.write(f"{chi2fib_med}\t")
                    f.write(f"{frac_c2}\t")
                    f.write(f"{frac_0}\t")
                    f.write(f"{im_median}\t")
                    f.write(f"{MaskFraction}\t")
                    f.write(f"{sky_sub_rms}\t")
                    f.write(f"{sky_sub_rms_rel}\t")
                    f.write(f"{sky_sub_rms_median_exp}\t")
                    f.write(f"{dither_relflux}\t")
                    f.write(f"{norm}\t")
                    #f.write(f"{kN_c}\t")
                    f.write(f"{kchi}")
                    f.write(f"{N_cont}")
                    #f.write(f"{n_lo_ss}\t")
                    f.write(f"\n")
                except:
                    print("stats_save_as ", print(traceback.format_exc()))
                    continue


    return tab



# def stats_save_as_fits(shot_dict):
#     """
#     just save the shot_dict as a table (in ascii format so is easy to cat files together) named for the shotid
#     """
#     tab = stats_shot_dict_to_table(shot_dict)
#     tab.write(shotid2datevshot(shot_dict['shotid'])+"_ampstats.dat",overwrite=True,format="ascii")



def save_shot_stats_pickle(shot_dict):
    """

    Parameters
    ----------
    shot_dict - single shotdict or list

    Returns
    -------
    None
    """

    try:
        # just save as a pickle

        if isinstance(shot_dict, dict):  # assume a single entry
            shot_dict_array = [shot_dict]
        else:
            shot_dict_array = shot_dict

        for sd in shot_dict_array:
            fn = f"{sd['shotid']}_stats.pickle"
            with open(fn, "wb+") as f:
                pickle.dump(sd, f)

    except Exception as e:
        # todo: error handling
        print(f"save_shot_stats_pickle  ", print(traceback.format_exc()))


def get_shotids_in_path(path="./"):
    """
    load a list of shotids that correspond to pickle files in the specified directory
    Parameters
    ----------
    path

    Returns
    -------
    list of shotids

    """

    shotids = []
    try:
        fns = sorted(glob.glob(op.join(path,"*_stats.pickle")))
        for fn in fns:
            try:
                shotid = int(op.basename(fn).split("_")[0])
                shotids.append(shotid)
            except:
                print(f"Invalid name: {fn}")

    except:
        print(f"get_shotids_in_path ", print(traceback.format_exc()))

    return shotids

def load_shot_stats_pickle(shotid,path="./"):
    """

    Parameters
    ----------
    shotid - single integer style shotid or an array of shotids
    path = path underwhich to look for the pickle files

    Returns
    -------
        single dictionary or array of dictionaries to match the passed in shot
    """
    isarray = True
    shot_dict_array = []
    try:

        if isinstance(shotid, (int,np.integer)):  # assume a single entry
            shotid = [shotid]
            isarray = False #a single value was passed in

        for s in shotid:
            try:
                fn = f"{s}_stats.pickle"
                with open(op.join(path,fn), "rb") as f:
                    shot_dict = pickle.load(f)
                    shot_dict_array.append(shot_dict)
            except:
                pass
                #print("load_shot_stats_pickle ", print(traceback.format_exc()))

        if isarray or len(shot_dict_array) == 0:
            return shot_dict_array
        else: #should not be able to have isarray == False and len(shot_dict_array) > 1
            return shot_dict_array[0]

    except Exception as e:
        # todo: error handling
        print("load_shot_stats_pickle ", print(traceback.format_exc()))
        return None






def stats_update_shot(h5, shot_dict=None, shot_dict_tab=None):
    """
    Takes a stats dictionary and updates the associated shot

    new group, organized by multiframe + exposure, extra columns are the statistics/metrics and the interpreted results
    """


    #DD 20240905 This is currently out of date
    #DD 202509011 Updated to match last AmpStats table ... intended for post-HETDEX shot.h5 use
    #print("!!!!! This needs to be updated with new columns  !!!!!")

    try:
        if shot_dict_tab is not None:
            tab = shot_dict_tab
        elif shot_dict is not None:
            tab = stats_shot_dict_to_table(shot_dict)
        else:
            #this is a problem
            print("Invalid parameters. Neither shot_dict nor shot_dict_tab provided.")
            return
    except:
        print("Invalid parameters. Neither shot_dict nor shot_dict_tab provided.")
        return

    try:
        # create the AmpStats table if it does not exist
        # Note: if exists, we want to UPDATE the rows and append if new entry
        # NOTE!!! this does NOT create new columns
        _ = h5.root.__getattr__('AmpStats')
        create_tab = False
        print("AmpStats table already exists. Will update.")
    except:
        create_tab = True
        print("AmpStats table does not exist. Will create.")

        # print("Creating AmpStats")

        class AmpStats(tables.IsDescription):
            # shotid ... do not need shotid
            multiframe = tables.StringCol(itemsize=20, pos=0)
            expnum = tables.Int32Col(pos=1)  # !! cannot be int8 as it is necessarily xlated to BOOL by fits i/o
            status = tables.Int32Col(pos=2)  # a status indicator, TBD ... could be a value or a bitmapped mask (-1 bad, 0 unchecked, 1 good?)
            im_median = tables.Float32Col()
            mask_fraction = tables.Float32Col()
            avg = tables.Float32Col()
            scale = tables.Float32Col()
            chi2fib_med = tables.Float32Col()
            frac_c2 = tables.Float32Col()
            frac_0 = tables.Float32Col()

            n_lo = tables.Int32Col()
            avg_orig = tables.Float32Col()
            sky_sub_rms = tables.Float32Col()
            sky_sub_rms_rel = tables.Float32Col()
            sky_sub_rms_median = tables.Float32Col()
            dither_relflux = tables.Float32Col()
            norm = tables.Float32Col()
            kchi = tables.Float32Col()
            n_cont = tables.Int32Col()
            #date ... do not need date
            flag = tables.Int32Col()
            #flag_manual = tables.Int32Col()
            #flag_manual_desc ... cannot have that here

        if create_tab:
            h5.create_table(
                h5.root,
                "AmpStats",
                AmpStats,
                "Statistics for each Amp per exposure for this shot",
                expectedrows=1000,  # 78 IFU x 4 Amps x 3 Exposures = 936
            )

    astb = h5.root.AmpStats


    if create_tab:
        # all new rows
        for entry in tab:
            row = astb.row
            row['multiframe'] = entry['multiframe']
            row['expnum'] = entry['expnum']
            row['status'] = 0  # TBD
            row['im_median'] = entry['im_median']
            row['mask_fraction'] = entry['MaskFraction']
            row['avg'] = entry['Avg']
            row['scale'] = entry['Scale']
            row['chi2fib_med'] = entry['chi2fib_med']
            row['frac_c2'] = entry['frac_c2']
            row['frac_0'] = entry['frac_0']
            row['n_lo'] = entry['n_lo']
            row['avg_orig'] = entry['Avg_orig']
            row['sky_sub_rms'] = entry['sky_sub_rms']
            row['sky_sub_rms_rel'] = entry['sky_sub_rms_rel']
            row['sky_sub_rms_median'] = entry['sky_sub_rms_median']
            row['dither_relflux'] = entry['dither_relflux']
            row['norm'] = entry['norm']
            row['kchi'] = entry['kchi']
            row['n_cont'] = entry['N_cont']
            try:
                row['flag'] = entry['flag']
            except:
                row['flag'] = -1 #unset
                print("stats_update_shot (1)", print(traceback.format_exc()))

            row.append()

        astb.flush()
        astb.cols.multiframe.create_csindex()

    else:
        # update matched rows and ADD new ones if multiframe+exp does not exist
        new_rows = False
        for entry in tab:
            q_mf = entry['multiframe']
            q_exp = entry['expnum']
            # find the counterpart
            row = astb.read_where("(multiframe==q_mf) & (expnum==q_exp)")

            if len(row) == 0:  # new row
                print("Adding new row")
                new_rows = True
                row = astb.row
                row['multiframe'] = entry['multiframe']
                row['expnum'] = entry['expnum']
                row['status'] = 0  # TBD
                row['im_median'] = entry['im_median']
                row['mask_fraction'] = entry['MaskFraction']
                row['avg'] = entry['Avg']
                row['scale'] = entry['Scale']
                row['chi2fib_med'] = entry['chi2fib_med']
                row['frac_c2'] = entry['frac_c2']
                row['frac_0'] = entry['frac_0']
                row['n_lo'] = entry['n_lo']
                row['avg_orig'] = entry['Avg_orig']
                row['sky_sub_rms'] = entry['sky_sub_rms']
                row['sky_sub_rms_rel'] = entry['sky_sub_rms_rel']
                row['sky_sub_rms_median'] = entry['sky_sub_rms_median']
                row['dither_relflux'] = entry['dither_relflux']
                row['norm'] = entry['norm']
                row['kchi'] = entry['kchi']
                row['n_cont'] = entry['N_cont']
                try:
                    row['flag'] = entry['flag']
                except:
                    print("stats_update_shot (2)", print(traceback.format_exc()))
                    row['flag'] = -1  # unset


                row.append()
            elif len(row) == 1:  # update
                # need to get row iterator now, so the old (numpy) row is replaced
                # with a row from the iterator
                for row in astb.where("(multiframe==q_mf) & (expnum==q_exp)"):
                    row['status'] = 0  # TBD
                    row['im_median'] = entry['im_median']
                    row['mask_fraction'] = entry['MaskFraction']
                    row['avg'] = entry['Avg']
                    row['scale'] = entry['Scale']
                    row['chi2fib_med'] = entry['chi2fib_med']
                    row['frac_c2'] = entry['frac_c2']
                    row['frac_0'] = entry['frac_0']
                    row['n_lo'] = entry['n_lo']
                    row['avg_orig'] = entry['Avg_orig']
                    row['sky_sub_rms'] = entry['sky_sub_rms']
                    row['sky_sub_rms_rel'] = entry['sky_sub_rms_rel']
                    row['sky_sub_rms_median'] = entry['sky_sub_rms_median']
                    row['dither_relflux'] = entry['dither_relflux']
                    row['norm'] = entry['norm']
                    row['kchi'] = entry['kchi']
                    row['n_cont'] = entry['N_cont']
                    try:
                        row['flag'] = entry['flag']
                    except:
                        row['flag'] = -1  # unset
                        print("stats_update_shot (3)", print(traceback.format_exc()))


                    row.update()
            else:  # problem
                pass  # todo, error control, too many rows

        astb.flush()
        if new_rows:  # need to clear and reset the index
            astb.cols.multiframe.remove_index()
            astb.cols.multiframe.create_csindex()


#########################################
# dictionaries
#########################################

def stats_amp_dict(shotid, multiframe, num_exposures=3, expid=None):
    """
    return a new stats dictionary, this is the functional unit

    input: optional the shot information or shot h5 to pre-populate dictionary info

    """

    def amp_exp_dict(shotid, multiframe, expid):
        # todo: statistics go here (with empty initializations)

        exp_dict = {
            "shotid": shotid,
            "multiframe": multiframe,
            "expid": expid,  # 1,2,3
            # todo: various statistics
        }

        return exp_dict

    try:

        # validate:
        amp = multiframe[-2:]

        if amp.upper() not in AMP_LIST:
            return {"status": -1, "reason": f"invalid amp specifier: ({amp})"}

        # at this point assume the shot info is okay .. can fail later when stats are done
        # format "multi_<specid>_<ifuslotid>_<ifuid>_<amp>"  i.e. multi_415_048_067_RL
        # or w/o "multi"
        toks = multiframe.split("_")
        if toks[0] == "multi":
            i = 1
        else:
            i = 0
            multiframe = f"multi_{multiframe}"  # always want the prefix "multi_"

        spec = toks[i]  # leave as a string
        slot = toks[i + 1]
        ifu = toks[i + 2]
        amp = toks[i + 3]

        amp_dict = {
            "shotid": datevshot2int(shotid),
            "multiframe": multiframe,
            "specid": spec,  # aka camid or camera ID
            "slotid": slot,  # aka IFUSlotID
            "ifuid": ifu,
            "amp": amp,

            "exposures": [amp_exp_dict(shotid, multiframe, expid + 1) for expid in range(num_exposures)],
            # todo: various roll up statistics over all shots (if any)
        }

        return amp_dict

    except Exception as e:
        # todo: error handling
        print("stats_amp_dict", print(traceback.format_exc()))


def stats_ifu_dict(shotid, multiframe=None, num_exposures=3, expid=None, amp_dict_array=None):
    """
    return a new stats dictionary

    input: optional the shot information or shot h5 to pre-populate dictionary info

    """

    try:

        if multiframe is not None:
            amp_dict_array = []  # actually a list

            if multiframe[-2:].upper() in AMP_LIST:
                stripped_mf = multiframe[:-3]  # include the "_"
            else:
                stripped_mf = multiframe

            for ampid in AMP_LIST:
                mf = stripped_mf + f"_{ampid}"
                if expid is not None:
                    amp_dict = stats_amp_dict(shotid, mf, 1)
                    amp_dict['exposures'][0]['expid'] = expid
                    amp_dict_array.append(amp_dict)
                else:
                    amp_dict_array.append(stats_amp_dict(shotid, mf, 3))  # assume exposures 1,2,3
        else:  # assume the amp_dict_array is good, else the exception will trigger
            pass

        ifu_dict = {
            "shotid": amp_dict_array[0]['shotid'],
            # "multiframe": amp_dict_array[0]['multiframe'],
            "specid": amp_dict_array[0]['specid'],  # aka camid or camera ID
            "slotid": amp_dict_array[0]['slotid'],  # aka IFUSlotID
            "ifuid": amp_dict_array[0]['ifuid'],
            # just a list of exposures numbers, actual data is in amp_dict_array, assume all are the same
            "expid_array": np.array([a['expid'] for a in [amp['exposures'] for amp in amp_dict_array][0]]),

            "amp_dict_array": np.array(amp_dict_array),
        }

        return ifu_dict

    except Exception as e:
        # todo: error handling
        print("stats_ifu_dict", print(traceback.format_exc()))


def stats_shot_dict(h5, expid=None, ifu_dict_array=None, multiframes=None, expid_array=None):
    """
    return a new stats dictionary

    if expid supplied, just load that exposure
    if ifU_dict_array is supplied, just load those IFUs and exposures
    otherswise, load it all

    input: optional the shot information or shot h5 to pre-populate dictionary info

    """

    try:

        shotid = datevshot2int(op.basename(h5.filename).split('.')[0])
        status = 0
        # notice, this is a little backward from the IFU and Amp loading
        if ifu_dict_array is None:
            # load from the h5 file
            # there are fewer Images entries, so read multiframe from there
            # [:-3] strips off the _<AMP> since we only care about the full IFU here
            if multiframes is None:
                multiframes = np.unique([mf.decode()[:-3] for mf in h5.root.Data.Images.read(field="multiframe")])
            if len(multiframes) == 0:
                status = -1
            # normally about 78x4 = 312 multiframes

            ifu_dict_array = []

            if expid is None:  # load them all
                if expid_array is None:
                    expid_array = h5.root.Shot.read(field="expnum")[0]
                num_exposures = len(expid_array)

                # stats_ifu_dict(shotid, multiframe=None,num_exposures=3,expid=None,amp_dict_array=None):
                for mf in multiframes:
                    ifu_dict = stats_ifu_dict(shotid, multiframe=mf, num_exposures=num_exposures, expid=None)
                    ifu_dict_array.append(ifu_dict)
            else:
                expid_array = [expid]
                for mf in multiframes:
                    ifu_dict = stats_ifu_dict(shotid, multiframe=mf, num_exposures=1, expid=expid)
                    ifu_dict_array.append(ifu_dict)

        else:
            expid_array = ifu_dict_array[0]["expid_array"]

        # now we have an ifu_dict_array

        shot_dict = {
            "shotid": shotid,
            "status": status,
            "expid_array": expid_array,
            "ifu_dict_array": ifu_dict_array,

        }

        return shot_dict
    except Exception as e:
        # todo: error handling
        print("stats_shot_dict", print(traceback.format_exc()))






###############################
# statistics collection
###############################


def stats_amp(h5, multiframe=None, expid=None, amp_dict=None, fibers_table=None, images_table=None, shot_table=None):
    """
    This is the primary functional unit; most statistics/metrics are computed for each
      individual amp+exposure.

    Inputs: 1. needs existing h5 handle so can just open once and pass around
            2. amp identifier (all 3 or just any one?)
            3. exposure identifier

    For the given h5 file, either use the multiframe with or without the expid -OR-
                           use the amp_dict and populate it
    In the first case, if multiframe is specified but expid is NOT, then use all exposures

    Returns populated amp_dict

    """

    try:

        shotid = datevshot2int(op.basename(h5.filename).split('.')[0])

        if multiframe is not None:
            if expid is not None:
                # def stats_amp_dict(shotid, multiframe, num_exposures=3):
                amp_dict = stats_amp_dict(shotid, multiframe, 1)
                amp_dict['exposures'][0]['expid'] = expid
            else:
                amp_dict = stats_amp_dict(shotid, multiframe, 3)

        # now we have a base amp_dict (or there is an exception)

        #################################
        # collect amp+exp statistics
        #################################

        #######################################
        # based on Erin stats from
        #  /home1/05350/ecooper/work/stampede2/hdr3/get_im_values.py
        #######################################

        target_mf = amp_dict['multiframe']

        for exp_dict in amp_dict['exposures']:
            try:

                target_expnum = exp_dict['expid']
                query = f"((multiframe==target_mf) & (expnum==target_expnum))"

                if images_table is None:
                    #print("amp level read")
                    image = h5.root.Data.Images.read_where(query,field="image")
                else:
                    image_sel = np.array(images_table['multiframe']==target_mf) & np.array(images_table['expnum']==target_expnum)
                    image = images_table['image'][image_sel]
                #really only need "image" and "error" (don't use "clean_image" so there is some waste

                if len(image) != 1:
                    del image
                    #from original amp.dat
                    #exp_dict['kNorm'] = np.nan  #NOT COMPUTED YET
                    #exp_dict['kN_c'] = -1 #integer, so -1 is unset
                    exp_dict['Avg'] = np.nan
                    exp_dict['Scale'] = np.nan
                    #exp_dict['kW0'] = np.nan
                    #exp_dict['kW1'] = np.nan
                    exp_dict['Avg_orig'] = np.nan

                    exp_dict['kchi'] = np.nan
                    #exp_dict['kNlo'] = np.nan #same as n_lo


                    #in both
                    exp_dict['chi2fib_avg'] = np.nan
                    exp_dict['frac_c2'] = np.nan
                    exp_dict['frac_0'] = np.nan
                    exp_dict['n_lo'] = -1 # this is an integer, 0 is lowest possible so -1 is unset
                   # exp_dict['n_lo_ss'] = -1
                    exp_dict['N_cont'] = -1

                    #newer only
                    exp_dict['maskfraction'] = np.nan
                    exp_dict['im_median'] = np.nan
                    exp_dict['sky_sub_rms'] = np.nan
                    exp_dict['sky_sub_rms_median'] = np.nan #gets added on shot rollup
                    exp_dict['sky_sub_rms_rel'] = np.nan  # gets added on shot rollup

                    # exp_dict['sky_sub_bws'] = -999.0
                    exp_dict['dither_relflux'] = np.nan #this is used for Norm
                                                        #which is the max / min
                    exp_dict['norm'] = np.nan
                else:

                    if images_table is None:
                        #image does == 1, but is 1x 1032x 1032 so just want the single
                        image = image[0]
                        #now also need the error for the image
                        #print("amp level read")
                        error = h5.root.Data.Images.read_where(query, field="error")[0]
                        #NOTICE: image and error are just too big to read in ALL of them up front,
                        #        so still need to read them just for the amp+exmposure here
                    else: #note we aleady have image
                        error = images_table['error'][image_sel][0]
                        image = image[0]

                    ###############################
                    # used in various calcs below
                    ###############################

                    if fibers_table is None:
                        clean_tab = True
                        #print("amp level read")
                        tab = Table(h5.root.Data.Fibers.read_where(query))  # '(multiframe == mf) & (expnum == expi)'))
                    else:
                        clean_tab = False
                        tab = fibers_table

                    # calfib = h5.root.Data.Fibers.read_where(query, field="calfib")
                    # calfib_counts = h5.root.Data.Fibers.read_where(query, field="calfib_counts")
                    # calfib_ffsky = h5.root.Data.Fibers.read_where(query, field="calfib_ffsky")
                    #
                    # calfib_ffsky_counts = flx2cts(calfib, calfib_counts, calfib_ffsky)
                    #
                    # #used as approximates for some older amp.dat columns
                    # chi2 = h5.root.Data.Fibers.read_where(query, field="chi2")
                    # spectrum = h5.root.Data.Fibers.read_where(query, field="spectrum")
                    # sky_subtracted = h5.root.Data.Fibers.read_where(query, field="sky_subtracted")
                    #
                    #


                    tab_sel = np.array(tab['multiframe']==target_mf) & np.array(tab['expnum']==target_expnum)

                    calfib = np.array(tab['calfib'][tab_sel])
                    #calfibe = np.array(tab['calfibe'][tab_sel])
                    calfib_counts = np.array(tab['calfib_counts'][tab_sel])
                    calfib_ffsky = np.array(tab['calfib_ffsky'][tab_sel])

                    calfib_ffsky_counts = flx2cts(calfib, calfib_counts, calfib_ffsky)

                    #used as approximates for some older amp.dat columns
                    chi2 = np.array(tab['chi2'][tab_sel])
                    spectrum = np.array(tab['spectrum'][tab_sel])
                    sky_subtracted = np.array(tab['sky_subtracted'][tab_sel])
                    #fiber_to_fiber = np.array(tab['fiber_to_fiber'][tab_sel])


                    ############################################
                    #try to get close to original amp.dat cols
                    ############################################

                    #exp_dict['kNorm'] = np.nan                      #NOT COMPUTED YET
                    #exp_dict['kN_c'] = -1 #integer, so -1 is unset #NOT COMPUTED YET
                    exp_dict['Avg'] = np.nan
                    exp_dict['Scale'] = np.nan
                    #exp_dict['kW0'] = np.nan                        #NOT COMPUTED YET, unimportant, unsued
                    #exp_dict['kW1'] = np.nan                        #NOT COMPUTED YET, unimportant, unsued
                    exp_dict['Avg_orig'] = np.nan

                    exp_dict['kchi'] = np.nan
                    #exp_dict['kNlo'] = np.nan #same as n_lo

                    # #Nlo  #same as n_lo later
                    # M = calfib_ffsky_counts[:,299:800]
                    # M0 = M[M != 0]
                    # if np.size(M0) > 1:
                    #     sddev = biweight.biweight_scale(M0)
                    #
                    #     ct = 0
                    #     for i in range(len(M)):
                    #         f = M[i][299:800]
                    #         if np.count_nonzero(f) != 0:
                    #             if biweight.biweight_location(f[f!=0]) <(-2 * sddev):
                    #                 ct += 1
                    #         #else:
                    #         #    pass
                    #
                    #     exp_dict['kNlo'] = ct
                    # else:
                    #     exp_dict['kNlo'] = np.nan

                    #Scale
                    #a2a is 1036 long, but spectrum is 1032
                    #a2a = np.loadtxt("/home/jovyan/Hobby-Eberly-Telesco/lib_calib/202407/i067aLLata.dat",usecols=[1],dtype=float)
                    #M = (a2a*calfib_ffsky_counts[sel])[9:100,299:800]
                    M = calfib_ffsky_counts[9:100,299:800]
                    M0 = M[M != 0]
                    if np.size(M0) > 1:
                        sddev = biweight.biweight_scale(M0) #may need amp2amp, fiber2fiber normalization
                        exp_dict['Scale'] = sddev
                        exp_dict['Avg'] = biweight.biweight_location(M0)
                    else:
                        exp_dict['Scale'] = np.nan
                        exp_dict['Avg'] = np.nan


                    #Avg_orig (close but computation is not exact on the same array (post-corrections for amp2amp, fiber2fiber))
                    M = spectrum[9:100,299:800]
                    M0 = M[M != 0]
                    if np.size(M0) > 1:
                        exp_dict['Avg_orig'] = biweight.biweight_location(M0)
                    else:
                        exp_dict['Avg_orig'] = np.nan

                    #chi arrays 100 to 1000, 10:100
                    #    index   99   999
                    M = chi2[9:100,99:1000]
                    M0 = M[M > 0] #NOTICE >0 not !=0
                    if np.size(M0) > 1:
                        exp_dict['kchi'] = biweight.biweight_location(M0)
                    else:
                        exp_dict['kchi'] = np.nan

                    ######################################################
                    # new/replacement statistics apart from amp.dat
                    #######################################################

                    #im_media, maskfraction, chi2fib, frac_c2 based on:
                    #/home1/05350/ecooper/work/stampede2/hdr3/get_im_values.py

                    # similar to what original amp.dat called Scale, but over a different range
                    # Scale is from ffsky counts instead of the "image"
                    #      and for fibers 10 to 100 (inclusive), wavebins 300 to 800 (inclusive)
                    #      and explicitly excludes 0 values
                    #      and uses biweight instead of median
                    #image = row["image"][0]
                    exp_dict['im_median'] = np.nanmedian(image[300:800, 300:800])

                    # similar to original amp.dat Frac0
                    # Frac0 is from ffsky counts instead of the "image"
                    #      and for fibers 10 to 100 (inclusive), wavebins 300 to 800 (inclusive)
                    #      and is the fraction of those wavelength bin values that == 0
                    #error = row["error"][0]

                    sel_zero = np.array(error == 0)
                    exp_dict['maskfraction'] = np.sum(sel_zero) / np.size(sel_zero)
                    del sel_zero
                    #sel_reg = (tab['fibnum'] >= 10) & (tab['fibnum'] <= 100)
                    #chi2fib = tab['chi2'][sel_reg]
                    #chi2_arr = np.array(chi2fib[:, 100:1000])
                    #chi2_arr = np.array(chi2[9:100,299:800])
                    chi2_arr = np.array(chi2[9:100, 100:1000])

                    chi2_gtzero = chi2_arr[chi2_arr > 0]
                    # similar to the orginal amp.dat chi, with a +/-1 count/index difference on fibers and wavelengths
                    #   both uses fibers 10 to 100
                    #   both use the same wavelength bins (100:1000)
                    # Original chi vs this one are very close (within less than about 1%)
                    if np.size(chi2_gtzero) > 1:
                        exp_dict['chi2fib_avg'] = biweight.biweight_location(chi2_gtzero)

                        # combined these are similar to the Frac_c2, with a +/-1 count/index difference on fibers and wavelengths
                        #   both uses fibers 10 to 100
                        #   both use the same wavelength bins (100:1000)
                        # Original Frac_c2 vs this one are very close (within less than about 1%)
                        exp_dict['frac_c2'] = np.sum(chi2_arr > 2.) / np.sum(chi2_gtzero)
                    else:
                        exp_dict['chi2fib_avg'] = np.nan
                        exp_dict['frac_c2'] = np.nan

                    #####################################
                    # similar to original amp.dat stuff
                    #####################################



                    # because we are dividing out (in the shot level stats), it does not matter
                    #   if this is in flux or counts, however, for individual amp+exposures,
                    #   the thresholds are historically in counts
                    M = calfib_ffsky_counts[9:100, 299:800]

                    # a little off but okay (we think a better calculation than the original in amp.dat) ??
                    # Frac0
                    denom = np.size(M)
                    if denom != 0:
                        exp_dict['frac_0'] = np.size(M[M == 0]) / denom
                    else:
                        exp_dict['frac_0'] = np.nan

                    # sky_sub_rms  (also used with sky_sub_rms_rel)
                    M0 = M[M != 0]
                    if np.size(M0) > 1:
                        exp_dict['sky_sub_rms'] = biweight.biweight_scale(M0)
                    else:
                        exp_dict['sky_sub_rms'] = np.nan



                    if False:
                        #####################
                        # N_cont (compact)
                        #####################
                        M = calfib_ffsky_counts  #uses all 112 x 1032 (or 1036 in this case) # [:, 299:800]
                        #M = sky_subtracted #[:,66:-66] * fiber_to_fiber[:,66:-66]
                        M[M==0] = np.nan
                        M[calfibe==0] = np.nan
                        if np.any(M) and np.any(~np.isnan(M)):
                            sddev = biweight.biweight_scale(M,ignore_nan=True) / 30.0 #i.e. np.sqrt(900) #yes, 900 not 1032 or 1036, per Karl ... empirircal
                            avg = biweight.biweight_location(M,ignore_nan=True)

                            bwl_f = biweight.biweight_location(M,axis=1,ignore_nan=True)
                            hi_ct = np.count_nonzero(bwl_f > (avg + 7.0 * sddev))

                            exp_dict['N_cont'] = hi_ct
                        else:
                            exp_dict['N_cont'] = -1

                    else:
                        ###############################
                        #N_cont iterative
                        ###############################
                        #M = calfib_ffsky_counts  # uses all 112 x 1032 (or 1036 in this case) # [:, 299:800]
                        M = sky_subtracted #[:,66:-66] * fiber_to_fiber[:,66:-66]
                        M0 = M[M != 0]
                        #Msel = np.array(M != 0) & np.array(calfibe != 0)
                        #M0 = M[Msel] #slightly better with this 0 filter
                        #M0 = M[calfibe !=0]
                        #M0 = M0[M0 != 0]

                        if np.size(M0) > 1:
                            sddev = biweight.biweight_scale(M0,ignore_nan=True) / 30.0  # i.e. np.sqrt(900) #yes, 900 not 1032 or 1036, per Karl ... empirircal
                            avg = biweight.biweight_location(M0, ignore_nan=True)

                            hi_ct = 0
                            for i in range(len(M)): #work down all the 112 fibers
                                f = M[i]
                                #e = calfibe[i]
                                #f = f[e!=0]
                                f = f[f!=0]
                                if np.count_nonzero(f) != 0:

                                    bwl_f = biweight.biweight_location(f)
                                    if bwl_f > (avg + (7.0 * sddev)):
                                    #if bwl_f > (0.0 + (7.0 * sddev)):  #DEFINITELY NOT THIS ONE
                                        hi_ct += 1
                                #else: #all zero? mark it?
                                #    hi_ct += 1

                            exp_dict['N_cont'] = hi_ct
                        else:
                            exp_dict['N_cont'] = -1



                    if False:
                        # Nlo (different M selection than frac_0)
                        ##############
                        # n_lo
                        ##############
                        # Nlo (different M selection than frac_0)
                        # again, this is different than the original amp.dat, but may be a better calculation
                       # M = calfib_ffsky_counts[:, 299:800] #uses clipped interior bit 300 to 800 inclusisve per Karl, so 299:800 for Python
                        M = calfib_ffsky_counts[:,299:800] #* fiber_to_fiber[:,299:800]
                        M[M == 0] = np.nan
                        M[calfibe[:,299:800] == 0] = np.nan
                        #M0 = M[M != 0]
                        #if np.size(M0) > 1:
                        if np.any(M) and np.any(~np.isnan(M)):
                            sddev = biweight.biweight_scale(M,ignore_nan=True)
                            avg = biweight.biweight_location(M,ignore_nan=True)

                            bwl_f = biweight.biweight_location(M, axis=1, ignore_nan=True)
                            lo_ct = np.count_nonzero(bwl_f < (avg - 2.0 * sddev))

                            # lo_ct = 0
                            # for i in range(len(M)): #work down all the 112 fibers
                            #     f = M[i]
                            #     if np.count_nonzero(f) != 0:
                            #         bwl_f = biweight.biweight_location(f[f != 0])
                            #         if bwl_f < (avg - (2.0 * sddev)):
                            #             lo_ct += 1
                            #     # else: #all zero? mark it?
                            #     #     lo_ct += 1

                            exp_dict['n_lo'] = lo_ct
                        else:
                            exp_dict['n_lo'] = -1
                    else:
                        #iterative version
                        ##############
                        # n_lo
                        ##############
                        # Nlo (different M selection than frac_0)
                        # again, this is different than the original amp.dat, but may be a better calculation
                        # M = calfib_ffsky_counts[:, 299:800] #uses clipped interior bit 300 to 800 inclusisve per Karl, so 299:800 for Python
                        #M = calfib_ffsky_counts[:, 299:800]  # * fiber_to_fiber[:,299:800]
                        M = sky_subtracted[:, 299:800]
                        M0 = M[M!=0]
                        #M0 = M[calfibe[:, 299:800] != 0]
                        #M0 = M0[M0 != 0]
                        # M0 = M[M != 0]
                        if np.size(M0) > 1:
                            sddev = biweight.biweight_scale(M0, ignore_nan=True)
                            avg = biweight.biweight_location(M0, ignore_nan=True)

                            lo_ct = 0
                            for i in range(len(M)): #work down all the 112 fibers
                                f = M[i]
                                #e = calfibe[:, 299:800][i]
                                #f = f[e != 0]
                                f = f[f != 0]
                                if np.count_nonzero(f) != 0:
                                    bwl_f = biweight.biweight_location(f[f != 0])
                                    if bwl_f < (avg - (2.0 * sddev)):
                                        lo_ct += 1
                                # else: #all zero? mark it?
                                #     lo_ct += 1

                            exp_dict['n_lo'] = lo_ct
                        else:
                            exp_dict['n_lo'] = -1

                    ##############
                    # n_lo (alternate)
                    ##############
                    # Nlo (different M selection than frac_0)
                    # again, this is different than the original amp.dat, but may be a better calculation
                    #M = calfib_ffsky_counts[:, 299:800]
                    # M = sky_subtracted[:, 299:800] #uses clipped interior bit 300 to 800 inclusisve per Karl, so 299:800 for Python
                    # M0 = M[M != 0]
                    # if np.size(M0) > 1:
                    #     sddev = biweight.biweight_scale(M0)
                    #     avg = biweight.biweight_location(M0)
                    #     lo_ct = 0
                    #     for i in range(len(M)): #work down all the 112 fibers
                    #         f = M[i]
                    #         if np.count_nonzero(f) != 0:
                    #             bwl_f = biweight.biweight_location(f[f != 0])
                    #             if bwl_f < (avg - (2.0 * sddev)):
                    #                 lo_ct += 1
                    #         #else: # the counts are ALL zero? while not technically "low"
                    #               # by this definition, it IS a problem (and really, all zero is "low")
                    #             #lo_ct += 1
                    #     exp_dict['n_lo_ss'] = lo_ct
                    # else:
                    #     exp_dict['n_lo_ss'] = -1

                    ########################################
                    # used at shot level (may have feedback)
                    ########################################
                    dither_flux = np.nan
                    try:
                        # local read
                        if shot_table is None:
                            dither_flux = h5.root.Shot.read(field="relflux_virus")
                        else:
                            dither_flux = shot_table['relflux_virus']
                        exp_dict['dither_relflux'] = dither_flux[0][int(target_expnum) - 1]
                        if np.any(np.isnan(dither_flux)):
                            exp_dict['norm'] = np.nan
                        else:
                            exp_dict['norm'] = np.nanmax(dither_flux)/np.nanmin(dither_flux)
                    except:
                        exp_dict['norm'] = np.nan
                        exp_dict['dither_relflux'] = np.nan
                        print(f"stats_amp statistics, dither_flux: {dither_flux}", print(traceback.format_exc()))

                    exp_dict['sky_sub_rms_rel'] = np.nan  # gets a value at shot level feedback/rollup

                    # cleanup
                    if clean_tab:
                        del tab

                    del calfib_ffsky_counts

            except Exception as e:
                # todo: handle fail for this amp+exposure
                print("stats_amp statisitcs (Erin)", print(traceback.format_exc()))



        return amp_dict
    except Exception as e:
        # todo: error handling
        print("stats_amp statisitcs", print(traceback.format_exc()))


def stats_ifu(h5, multiframe=None, expid=None, ifu_dict=None,fibers_table=None,images_table=None,shot_table=None):
    """

    Rollup of calls to stats_amp(). Normally 12 calls, one for each of 4 amps x 3 exposures

    operate on all amps for an IFU
    note: we IGNORE the amp on the multiframe if provided
    """
    try:

        shotid = datevshot2int(op.basename(h5.filename).split('.')[0])

        if multiframe is not None:
            if expid is not None:
                # def stats_amp_dict(shotid, multiframe, num_exposures=3):
                ifu_dict = stats_ifu_dict(shotid, multiframe, num_exposures=1, expid=expid)
            else:
                ifu_dict = stats_ifu_dict(shotid, multiframe, num_exposures=3)
        else:
            pass  # assume ifu_dict is already populated, otherwise will hit exception

        # now we have a base ifu_dict (or there is an exception)
        # so iterate over the amp_dicts and exposures
        # def stats_amp(h5,multiframe=None,expid=None,amp_dict=None):

        for amp_dict in ifu_dict['amp_dict_array']:
            # print(amp_dict)
            stats_amp(h5, multiframe=None, expid=None, amp_dict=amp_dict,
                      fibers_table=fibers_table,images_table=images_table,shot_table=shot_table)

        return ifu_dict
    except Exception as e:
        # todo: error handling
        print("stats_ifu statisitcs", print(traceback.format_exc()))



def make_stats_for_shot(shotid=None, survey=None,fqfn=None, save=True, preload=True):
    """

    Open a SINGLE shot file and build the full set of amp+exp and full shot statistics.
    This can be done via the shotid and a lookup -OR- by specifying a path to a specific shot h5 file

    NOTICE: there is NO option to append/update the new Group to the hdf5 file here! That needs to be used with
            caution and so requires an explicit set of function calls to do it.

    Parameters
    ----------
    shotid - SINGLE shotid in integer format
    survey - if not specified, uses whatever is configured for LATEST. If specified, is string: e.g "hdr4" or "hdr5", etc
    fqfn - path to an h5 file .. if specified us THIS h5 file. Ignores shotid and survey.
    save - if TRUE (default) save the stats dictionary as a pickle

    Returns
    -------

    """

    try:

        #get the config
        if fqfn is None:
            if survey:
                HETDEX_API_CONFIG = HDRconfig(survey=survey)
            else:
                HETDEX_API_CONFIG = HDRconfig() #use "latest"

            #load the h5
            h5 = stats_get_shot_h5(config=HETDEX_API_CONFIG, shotid=shotid, fqfn=None, append=False)
        else:
            h5 = stats_get_shot_h5(config=None, shotid=None, fqfn=fqfn, append=False)
            shotid = h5.root.Shot.read(field="shotid")[0]

        if h5 is not None:

            if preload:
                print(f"{shotid} preload read ...")
                t_fib, t_img, t_shot = stats_make_shot_tables(h5)
            else:
                t_fib, t_img, t_shot = None, None, None

            print(f"{shotid} building stats ...")
            shot_dict = stats_shot(h5,expid=None, shot_dict=None, rollup=True, fibers_table=t_fib,images_table=t_img, shot_table=t_shot)

            if save and shot_dict is not None:
                save_shot_stats_pickle(shot_dict)
                print(f"{shotid} done.")
            else:
                print(f"[{shotid}] Error. Could not assemble stats dictionary.")
            h5.close()
        else:
            print(f"[{shotid}] Error. Could not find shot h5 file.")
            shot_dict = None

        return shot_dict
    except:
        print("make_stats_for_shot", print(traceback.format_exc()))
        return None

def stats_shot(h5, expid=None, shot_dict=None, rollup=True,fibers_table=None,images_table=None,shot_table=None):
    """
    Rollup of calls to stats_ifu(), which calls stats_amp()
    Normally 78 ifus x 4 amps x 3 expsoures = 936 total calls to stats_amp()

    """

    try:
        if shot_dict is None:
            if fibers_table is not None:
                multiframes = np.unique([mf.decode()[:-3] for mf in np.array(fibers_table["multiframe"])])
            else:
                multiframes = None

            if shot_table is not None:
                expid_array = shot_table['expnum'][0] #usually [1,2,3]
            else:
                expid_array  = None

            shot_dict = stats_shot_dict(h5, expid,multiframes=multiframes,expid_array=expid_array)


        # operate on all IFUs
        #for ifu_dict in tqdm(shot_dict['ifu_dict_array']):
        for ifu_dict in shot_dict['ifu_dict_array']:
            ifu_dict = stats_ifu(h5, ifu_dict=ifu_dict,fibers_table=fibers_table,images_table=images_table,shot_table=shot_table)

        #
        # shot level statistics, over all amps and exposures
        #

        if rollup:
            shot_dict = stats_shot_rollup(h5, shot_dict)

        return shot_dict
    except Exception as e:
        # todo: error handling
        print("stats_shot statisitcs", print(traceback.format_exc()))


def stats_make_shot_tables(h5):
    """
    select the key fields from the h5 file and make a shot level table of the entire shot

    Parameters
    ----------
    h5

    Returns
    -------
    astropy tables or None, None
    """

    t_fibers, t_images, t_shot = None, None, None

    try:
        # expid_array = h5.root.Shot.read(field="expnum")
        # dither_flux = h5.root.Shot.read(field="relflux_virus")
        #
        # t_shot = Table([expid_array, dither_flux], names=["expnum", "relflux_virus"])


        #just read the whole thing. ... it is not that big and may have future useful  bits
        t_shot = Table(h5.root.Shot.read())

    except:
        t_shot = None
        print("stats_make_shot_tables", print(traceback.format_exc()))

    try:
        #t_fibers = Table(h5.root.Data.Fibers.read())
        #only keep the columns we need
        #there is extra data here, but fewer read calls ... is that a good trade off for TACC?
        #t_fibers.keep_columns(
        #    ['multiframe', 'expnum', 'fibnum','calfib', 'calfibe', 'calfib_counts', 'calfib_ffsky', 'chi2', 'spectrum',
        #    'sky_subtracted','fiber_to_fiber])


        #alternate (one column at time)
        multiframe = h5.root.Data.Fibers.read(field="multiframe")
        expnum = h5.root.Data.Fibers.read(field="expnum")
        fibnum = h5.root.Data.Fibers.read(field="fibnum")
        calfib = h5.root.Data.Fibers.read(field="calfib")
        calfib_counts = h5.root.Data.Fibers.read(field="calfib_counts")
        calfib_ffsky = h5.root.Data.Fibers.read(field="calfib_ffsky")
        chi2 = h5.root.Data.Fibers.read(field="chi2")
        spectrum = h5.root.Data.Fibers.read(field="spectrum")
        sky_subtracted = h5.root.Data.Fibers.read(field="sky_subtracted")

        #fiber_to_fiber = h5.root.Data.Fibers.read(field="fiber_to_fiber")
        #calfibe = h5.root.Data.Fibers.read(field="calfibe")

        t_fibers = Table([multiframe, expnum, fibnum, calfib,calfib_counts,calfib_ffsky,chi2,spectrum,sky_subtracted],
                  names=["multiframe", "expnum", "fibnum", "calfib","calfib_counts","calfib_ffsky","chi2","spectrum",
                         "sky_subtracted"])
        #,"fiber_to_fiber","calfibe"]
        # fiber_to_fiber,calfibe],


        t_fibers.add_index('multiframe')
        t_fibers.add_index('expnum')


        del multiframe
        del expnum
        del fibnum
        del calfib
        del calfib_counts
        del calfib_ffsky
        del chi2
        del spectrum
        del sky_subtracted

        #del fiber_to_fiber
        #del calfibe

    except Exception as e:
        # todo: error handling
        print("stats_make_shot_tables", print(traceback.format_exc()))
        t_fibers = None

    try:
        #this might be too much
        #cannot read all at once, so have to read the 4 columns we need individually
        images_mf = h5.root.Data.Images.read(field="multiframe")
        images_en = h5.root.Data.Images.read(field="expnum")
        images_im = h5.root.Data.Images.read(field="image")
        images_er = h5.root.Data.Images.read(field="error")

        t_images = Table([images_mf, images_en, images_im, images_er], names=["multiframe", "expnum", "image", "error"])
        t_images.add_index('multiframe')
        t_images.add_index('expnum')

        del images_mf
        del images_en
        del images_im
        del images_er

    except Exception as e:
        # todo: error handling
        print("stats_make_shot_tables", print(traceback.format_exc()))
        t_images = None

    return t_fibers, t_images, t_shot

def stats_shot_rollup(h5, shot_dict):
    # todo: shot specific or shot aggregate statistics
    try:
        # todo: shot specific or shot aggregate statistics
        T = stats_shot_dict_to_table(shot_dict)
        shot_dict['exposures'] = np.unique(T['expnum'])
        shot_dict['sky_sub_rms_median_exp'] = []
        dither_relflux_array = []  # note: don't really need to compute this HERE, but it will be used in the qc part later

        for exp in shot_dict['exposures']:
            sel = np.array(T['expnum'] == exp)
            if np.count_nonzero(sel) == 0:
                continue

            # dither norms (should be the SAME VALUE for each amp for a given exposure)
            # so a mean should be redundant, but just to be safe, do it anyway
            if np.all(np.isnan(T['dither_relflux'][sel])):
                dither_relflux_array.append(np.nan)
            else:
                dither_relflux_array.append(np.nanmean(T['dither_relflux'][sel]))

            x = np.nan_to_num(T['sky_sub_rms'][sel], nan=-999.0)
            data_sel = np.array(x != -999.)
            shot_dict['sky_sub_rms_median_exp'].append(np.nanmedian(x[data_sel]))

        shot_dict['sky_sub_rms_median_exp'] = np.array(shot_dict['sky_sub_rms_median_exp'])
        dither_relflux_array = np.array(dither_relflux_array)
        shot_dict['dither_relflux'] = dither_relflux_array

        # now, retroactively update the amps
        for ifu in shot_dict['ifu_dict_array']:
            for amp in ifu['amp_dict_array']:
                for exp in amp['exposures']:
                    sel = shot_dict['exposures'] == exp['expid']
                    if np.count_nonzero(sel) != 1:
                        print(f"Fail! Bad rms selection.", exp, shot_dict['exposures'])
                        exp['sky_sub_rms_rel'] = np.nan
                        # exp['sky_sub_bws_rel'] = -999.0
                        continue
                    try:
                        #different from Erin's original method in that the median used is the median for THIS exposure
                        #rather than the median over all 3 exposures
                        exp['sky_sub_rms_median'] = shot_dict['sky_sub_rms_median_exp'][sel][0] #record the one we used
                        exp['sky_sub_rms_rel'] = exp['sky_sub_rms'] / exp['sky_sub_rms_median']
                        # exp['sky_sub_bws_rel'] = exp['sky_sub_bws'] / shot_dict['sky_sub_bws_median_exp'][sel][0]
                    except Exception as e:
                        exp['sky_sub_rms_rel'] = np.nan
                        # exp['sky_sub_bws_rel'] = -999.0
                        print("stats_shot statistics, sky_sub", print(traceback.format_exc()))
        return shot_dict

    except Exception as e:
        # todo: error handling
        print("stats_shot statistics", print(traceback.format_exc()))



#######################################
# QC stuff
#######################################

def is_masked(x): #masked or NaN
    #if ANY entry is "masked" then the whole is "masked" and you
    #only get back a SINGLE value (T/F) when running np.ma.is_masked
    #so have to use list comprehension, which is comparitively slow
    #NOTE: however, this gives the same results as just:  np.array(np.isnan(x))
    #  which turns masked into NaN and yields True in those cases
    #return np.isnan(x) | np.array([np.ma.is_masked(y) for y in x])
    return np.array(np.isnan(x)) #| np.array([np.ma.is_masked(y) for y in x])

def stats_qc(data,extend=False,total_exp_time=None):
    """
    todo: what constitutes a bad amp ...
    evaluate and set flag in the amp_dict

    Parameters
    ----------
    data - an stats dictionary or row from a table for a single amp+exposure
    extend - if TRUE, add a flag column to the Table representation of data and return as the Table
             (NOTICE: this IMPLICITY requires the return of an astropy Table
    Returns
    -------
    single value OR array of values for status

    """

    exp_time_norm = 1200.0 #seconds (20 min) ... normalization for some values or thresholds scaled to the total time
                           #nominal HETEX is 1080 secs but can be a bit longer under poorer conditions
    single = True #used later to return a single or iterable objects
    if isinstance(data,Row):
        amp_stats = Table(data)
    elif isinstance(data,Table):
        #this is for one row, which one?
        single = False
        amp_stats = data
    elif isinstance(data,dict):
        #could be a single amp+exp dict or a whole shot
        if 'ifu_dict_array' in data.keys(): #this is a whole shot
            amp_stats = stats_shot_dict_to_table(data)
            single = False
        elif 'amp_dict_array' in data.keys(): #this is an IFU
            print("Not ready for IFU only") # not read for this one yet
            single = False
        elif 'exposures' in data.keys(): #this is an amp with exposures
            print("Not ready for amp with exposures")  # not read for this one yet
            single = False
        else:
            amp_stats = Table([data])
    else: #assume a list or array of dicts, but this could fail
        amp_stats = Table(data)
        single = False




    #need date for subselection
    if 'shotid' in amp_stats.colnames:
        #this is the dateshot as integer format
        amp_stats['date'] = [int(str(r['shotid'])[0:8]) for r in amp_stats]
    elif 'dateshot' in amp_stats.colnames:
        #this should look like "d20230904s011exp01"
        amp_stats['date'] = [int(r['dateshot'][1:9]) for r in amp_stats]

    flags = np.zeros(len(amp_stats))


    #cloned from Erin's bad_amp_hunting notebook
    #note: 'background' is 'Avg'
    #'MaskFraction' may be 'maskfraction', but is the newer definition (not MaskFraction_karl)
    if 'maskfraction' in amp_stats.colnames:
        amp_stats.rename_column('maskfraction','MaskFraction')
    #"N_cont' is not currently computed ... not sure what to do with it
    #'nfib_bad' is not currently computed ... may be n_lo


    #todo: should we be allowing NaNs? if we could not compute the value and it is NaN, doesn't that
    #      already mean it is a problem??

    if True: #allow Nan or masked

        # "Original"
        # sel1 = ((amp_stats['Avg'].astype(float) > -10) *
        #         (amp_stats['Avg'].astype(float) < 100)) | (is_masked(amp_stats['Avg']))
        # sel2 = (amp_stats['sky_sub_rms_rel'] < 1.5) | (is_masked(amp_stats['sky_sub_rms_rel']))
        # sel3 = (amp_stats['sky_sub_rms'] > 0.2) | (is_masked(amp_stats['sky_sub_rms_rel']))
        # sel4 = (amp_stats['im_median'] > 0.05) | (is_masked(amp_stats['im_median']))
        # sel5 = (amp_stats['MaskFraction'] < 0.25) | (is_masked(amp_stats['MaskFraction']))
        # # sel6 = (amp_stats['N_cont'] < 35) | (is_masked(amp_stats['N_cont']))
        # sel6 = np.full(len(amp_stats), True)
        # # sel7 = (amp_stats['nfib_bad'] <= 1) | (is_masked(amp_stats['nfib_bad']))
        # sel7 = np.full(len(amp_stats), True)
        #
        # sel8 = (amp_stats['norm'] > 0.5) | (np.isnan(amp_stats['norm']))
        #
        # sel9_hdr3 = ((amp_stats['frac_c2'] < 0.5) | (is_masked(amp_stats['frac_c2']))) * (amp_stats['date'] < 20210901)
        # # next is hdr4 and ABOVE (some changes from HDR3 to HDR4 necessitate a different selection here)
        # sel9_hdr4 = ((amp_stats['frac_c2'] < 0.1) | (is_masked(amp_stats['frac_c2']))) * (amp_stats['date'] >= 20210901)
        # sel9 = sel9_hdr3 | sel9_hdr4

        if total_exp_time is None or total_exp_time < exp_time_norm:
            sel1 = ((amp_stats['Avg'].astype(float) > -5.0) *
                    (amp_stats['Avg'].astype(float) < 100.0)) | (is_masked(amp_stats['Avg']))
        else:
            avg_max = 100.00 * np.sqrt(total_exp_time/exp_time_norm)
            sel1 = ((amp_stats['Avg'].astype(float) > -5.0) *
                    (amp_stats['Avg'].astype(float) < avg_max)) | (is_masked(amp_stats['Avg']))
        #down to a max of 20 to 25 helps a little


        #should have an lower limit too? this is the sky_sub_rms / median(all the sky_sub_rms for the shot)
        #originally was 1.5 max, but 1.7 is working better to not include good amps
        #sel2 = ( (amp_stats['sky_sub_rms_rel'] >= 0.0) & (amp_stats['sky_sub_rms_rel'] < 1.7)) | (is_masked(amp_stats['sky_sub_rms_rel']))
        sel2 = (amp_stats['sky_sub_rms_rel'] < 1.7) | (is_masked(amp_stats['sky_sub_rms_rel']))

        #should maybe have an upper limit? .. this is a biweight scale (so stddev like) ...
        # maybe an upper limit based on the Avg? but  again consider if there is a bright object ... could be a huge
        # range in counts that is legit ... SO NO. there should NOT be an upper limit ...tested, does not work
        #upper limit does not seem to help any ... only cuts into "good"
        #sel3 = ((amp_stats['sky_sub_rms'] > 0.2) & (amp_stats['sky_sub_rms'] < 40.0)) | (is_masked(amp_stats['sky_sub_rms']))
        #not finding this to be helpful ... the sky_sub_rms_rel works well though
        #sel3 = (amp_stats['sky_sub_rms'] > 0.2) | (is_masked(amp_stats['sky_sub_rms'])) #original

        #a low value is a good check, but a high value may not be... consider if we are on a large, bright galaxy,
        # then the median would be naturally high
        sel4 = (amp_stats['im_median'] > 0.05) | (is_masked(amp_stats['im_median']))

        #could consider lowering ... 25% default may be too high, 20% seems a little better
        sel5 = (amp_stats['MaskFraction'] < 0.20) | (is_masked(amp_stats['MaskFraction']))

        sel6 = (amp_stats['N_cont'] < 37)  #| (is_masked(amp_stats['N_cont'])) #original is 35 but 37 is a better balance
                                           #between flagging more amps and adding maybe good ones
        #sel6 = np.full(len(amp_stats), True)
        #sel7 = (amp_stats['nfib_bad'] <= 1) | (is_masked(amp_stats['nfib_bad']))
        sel7 = (amp_stats['n_lo'] <= 1)    #| (is_masked(amp_stats['n_lo']))
        #sel7 = np.full(len(amp_stats), True)

        #since norm is max(dither flux norms) / min(dither flux norms) MUST always be at least 1.0 (if the max == min)
        #though this seems unusual that max == min would be the case
        #if max/min is too large, there is likely something wrong. x3.0 seems to be about that spot
        sel8 = ((amp_stats['norm'] > 1.0) & (amp_stats['norm'] < 3.0)) | (np.isnan(amp_stats['norm']))

        if total_exp_time is None or total_exp_time <= exp_time_norm:
            sel9_hdr3 = ((amp_stats['frac_c2'] < 0.5) | (is_masked(amp_stats['frac_c2']))) * (amp_stats['date'] < 20210901)
            # next is hdr4 and ABOVE (some changes from HDR3 to HDR4 necessitate a different selection here)
            #sel9 for hdr4+ is pretty insensitive until you get to < 0.08 or so (that is < .1 vs < 0.5 is almost no difference)
            sel9_hdr4 = ((amp_stats['frac_c2'] < 0.1) | (is_masked(amp_stats['frac_c2']))) * (amp_stats['date'] >= 20210901)

        else:
            sel9_hdr3 = ((amp_stats['frac_c2'] < 0.5) | (is_masked(amp_stats['frac_c2']))) * (
                        amp_stats['date'] < 20210901)
            # next is hdr4 and ABOVE (some changes from HDR3 to HDR4 necessitate a different selection here)
            # sel9 for hdr4+ is pretty insensitive until you get to < 0.08 or so (that is < .1 vs < 0.5 is almost no difference)
            max_frac_c2 = 0.1 * np.sqrt(total_exp_time / exp_time_norm)
            sel9_hdr4 = ((amp_stats['frac_c2'] < max_frac_c2) | (is_masked(amp_stats['frac_c2']))) * (
                        amp_stats['date'] >= 20210901)

        sel9 = sel9_hdr3 | sel9_hdr4

        #for a sub selection of wavelegths and fibers, the fraction that are zero valued
        #(so frac_0 == 1.0 means ALL are zero)
        #this helps and picks up a few more, but I think it HAS to be here
        # does not matter if the threshold is 1.0 or down to 0.20 or maybe even lower
        # this one REALLY needs the masked check
        sel11 = (amp_stats['frac_0'] < 1.0) | (is_masked(amp_stats['frac_0']))

        #lower limit of about 0.5 seems okay, less than that maybe implies abonormally uniform data
        #upper limit somewhere 1.03 to 1.5 is good .. below 1.03  gets too aggressive
        if total_exp_time is None or total_exp_time <= exp_time_norm:
            sel12 = ((amp_stats['chi2fib_med'] > 0.5) & (amp_stats['chi2fib_med'] < 1.05)) | \
                    (is_masked(amp_stats['chi2fib_med']))
        else:
            max_chi2fib = 1.05 * total_exp_time/exp_time_norm #np.sqrt(total_exp_time/exp_time_norm)
            sel12 = ((amp_stats['chi2fib_med'] > 0.5) & (amp_stats['chi2fib_med'] < 1.05)) | \
                    (is_masked(amp_stats['chi2fib_med']))

        # Not useful (at least not compared to others ... or may be redundant with others)
        #n_lo does not seem to be useful ... SHOULD be a replacement for nfib_bad (sel7), but does not seem to work that way
        # produces WAY too many hits?
        #sel10 = (amp_stats['n_lo'] <= 1) | (is_masked(amp_stats['n_lo']))





        #Scale is not useful
        #sel13 =  ((amp_stats['Scale'] > 0) & (amp_stats['Scale'] < 30.0)) | (is_masked(amp_stats['Scale']))

        #Avg_orig
        # sel14 = ((amp_stats['Avg_orig'].astype(float) > 0.0) *
        #         (amp_stats['Avg_orig'].astype(float) < 5000.0)) | (is_masked(amp_stats['Avg_orig']))

        #basically redundant to chi2fib_med
        #sel15 = ((amp_stats['kchi'] > 0.5) & (amp_stats['kchi'] < 1.2)) | (is_masked(amp_stats['kchi']))

        #dither_relflux .. not useful ...probably redundant
        # sel16 = ((amp_stats['dither_relflux'] > 0.5) & (amp_stats['dither_relflux'] < 1.5)) | (
        #     is_masked(amp_stats['dither_relflux']))

    else: #Prohibit Nan or masked (counts as BAD)
        sel1 = ((amp_stats['Avg'].astype(float) > -10) *
                (amp_stats['Avg'].astype(float) < 100)) & (~is_masked(amp_stats['Avg']))
        sel2 = (amp_stats['sky_sub_rms_rel'] < 1.5) & (~is_masked(amp_stats['sky_sub_rms_rel']))
        sel3 = (amp_stats['sky_sub_rms'] > 0.2)  & (~is_masked(amp_stats['sky_sub_rms_rel']))
        sel4 = (amp_stats['im_median'] > 0.05 ) & (~is_masked(amp_stats['im_median']))
        sel5 = (amp_stats['MaskFraction'] < 0.25) & (~is_masked(amp_stats['MaskFraction']))
        # sel6 = (amp_stats['N_cont'] < 35) & (is_masked(amp_stats['N_cont']))
        sel6 = np.full(len(amp_stats), True)
        # sel7 = (amp_stats['nfib_bad'] <= 1) & (~is_masked(amp_stats['nfib_bad']))
        sel7 = np.full(len(amp_stats), True)
        sel8 = (amp_stats['norm'] > 0.5) & (~is_masked(amp_stats['norm']))

        sel9_hdr3 = (amp_stats['frac_c2'] < 0.5) & (~is_masked(amp_stats['frac_c2']))  * (amp_stats['date'] < 20210901)
        # next is hdr4 and ABOVE (some changes from HDR3 to HDR4 necessitate a different selection here)
        sel9_hdr4 = (amp_stats['frac_c2'] < 0.1) & (~is_masked(amp_stats['frac_c2'])) * (amp_stats['date'] >= 20210901)
        sel9 = sel9_hdr3 | sel9_hdr4

    # else: #ignore NaN
    #     sel1 = ((amp_stats['Avg'].astype(float) > -10) * (amp_stats['Avg'].astype(float) < 100))
    #     sel2 = amp_stats['sky_sub_rms_rel'] < 1.5
    #     sel3 = amp_stats['sky_sub_rms'] > 0.2
    #     sel4 = amp_stats['im_median'] > 0.05
    #     sel5 = amp_stats['MaskFraction'] < 0.25
    #     # sel6 = amp_stats['N_cont'] < 35)
    #     sel6 = np.full(len(amp_stats), True)
    #     # sel7 = amp_stats['nfib_bad'] <= 1
    #     sel7 = np.full(len(amp_stats), True)
    #     sel8 = amp_stats['norm'] > 0.5
    #
    #
    #     sel9_hdr3 = (amp_stats['frac_c2'] < 0.5)  * (amp_stats['date'] < 20210901)
    #     #next is hdr4 and ABOVE (some changes from HDR3 to HDR4 necessitate a different selection here)
    #     sel9_hdr4 = (amp_stats['frac_c2'] < 0.1)  * (amp_stats['date'] >= 20210901)
    #
    #     sel9 = sel9_hdr3 | sel9_hdr4

    #GOOD selection
    sel = sel1 & sel2 &  sel4 & sel5 & sel6 & sel7 & sel8 & sel9 & sel11 & sel12
    #sel = sel1 & sel2 & sel3 & sel4 & sel5         & sel7 & sel8 & sel9 & sel11 & sel12
    #sel = sel1 & sel2 & sel3 & sel6 & sel7 & sel8

    # sel_wiggles = amp_stats['ft_flag'] == 1
    #todo: re-activate manual flagging

    #sel_manual = amp_stats['flag_manual'] == 1
    #sel = sel & sel_manual
    #amp_stats['flag'] = (sel * sel_manual ).astype(int)

    flags[sel] = 1

    if extend:
        amp_stats['flag'] = flags.astype(np.int32) # !! cannot be int8 as it is necessarily xlated to BOOL by fits i/o

        if 'flag_manual' not in amp_stats.colnames:
            flag_col_idx = list(amp_stats.colnames).index('flag')
            amp_stats.add_column( np.full(len(amp_stats),np.int32(-1)),name="flag_manual",index=flag_col_idx+1)

            if 'flag_manual_desc' not in amp_stats.colnames:
                amp_stats.add_column(np.full( len(amp_stats), str(" ")*256), name="flag_manual_desc", index=flag_col_idx + 2)

        return amp_stats
    elif single:
        del amp_stats
        return flags[0]
    else:
        del amp_stats
        return flags


def stats_update_flag_manual(db, shotid, multiframe=None,expnum=None,flag_manual=-1, flag_manual_desc=None, savefmt=None,interactive=True):
    """
    Updates the flag_manual and/or flag_manual_desc column(s) for the row(s) matching the input
      shotid, (optional) multiframe, and (optional) expnum

    Creates the flag_manual and flag_manual_desc columns if they do not already exist

    Parameters
    ----------
    db - required, the amp data table to update or full path to table to update. NOTE: will NOT save if db is a table in memory.
    shotid - required
    multiframe [optional] if not specfied ALL mutliframes for a shot are updated
    expnum [optional] if not specified ALL exposures are updated; can be set independently of multiframe
    flag_manual -1 = unset, 0 = bad, 1 = good
    flag_manual_desc - string up to 256 characters; suggestion - have brief resason why the flag is set, who set it and when (date)
    savefmt = if None, do not save, otherwise save using the astropy format specified (str). Only applies if db is a string (filename)
    interactive = if True (default) prompt the user to confirm

    Returns the number of rows updated or -1 if an exception,
       !!!NOTICE!!! This modifies the table in memory  BUT DOES NOT write out to disk
    -------

    """

    try:
        rows_to_update = 0

        if isinstance(db, str): #assume this is a filename
            try:
                tab = Table.read(db)
            except Exception as e:
                if interactive:
                    print(f"Unable to open file {db}. Cannot update.")
                    print(traceback.format_exc())
                return -1
        elif isinstance(db,Table):
            tab = db
        else:
            if interactive:
                print("Unknown database type. Cannot update.")
            return -1


        if 'flag' not in tab.colnames:
            if interactive:
                print("Input is missing flagging. Will apply now ... ")
            tab = stats_qc(tab,extend=True)


        if 'flag_manual' not in tab.colnames:
            flag_col_idx = list(tab.colnames).index('flag')
            # !! cannot be int8 as it is necessarily xlated to BOOL by fits i/o
            tab.add_column( np.full(len(tab),np.int32(-1)),name="flag_manual",index=flag_col_idx+1)

            if 'flag_manual_desc' not in tab.colnames:
                tab.add_column(np.full( len(tab), str(" ")*256), name="flag_manual_desc", index=flag_col_idx + 2)

        sel = np.array(tab['shotid']==shotid)
        if multiframe is not None:
            sel = sel &  np.array(tab['multiframe']==multiframe)
        if expnum is not None:
            sel = sel & np.array(tab['expnum']==expnum)

        rows_to_update = np.count_nonzero(sel)
        if rows_to_update > 0:
            if interactive:
                i = input(f"Records to update = {rows_to_update}. Proceed (y/n)?")
                if len(i) > 0 and i.upper() == "Y":
                    proceed = True
                else:
                    proceed = False
            else:
                proceed = True

            if proceed:
                tab['flag_manual'][sel] = flag_manual
                if flag_manual_desc is not None:
                    tab['flag_manual_desc'][sel] = flag_manual_desc

                if savefmt is not None and isinstance(db,str):
                    tab.write(db,overwrite=True,format=savefmt)
            else:
                print("Update cancelled.")
                rows_to_update = -1
        elif interactive:
            print("No records match the input. No update made.")
            print(f"shotid {shotid}, multiframe {multiframe}, expnum {expnum}")

    except Exception as e:
        # todo: error handling
        print("stats_update_manual_flag", print(traceback.format_exc()))
        rows_to_update = -1

    return rows_to_update

def stats_make_simple_amp_flag_file(db,outname="amp_flag",overwrite=False):
    """
    Simplifies amp statistics to just shotid + multiframe + flag
    Saves as two output files (fits and ascii format tables)

    Parameters
    ----------
    db = astropy table or filename of the amp statistics to simplify
    outname = basename of the files to save. The ".fits" and ".tab" extensions are added in this code
    overwrite = set to True if you want to replace files that already exist

    Returns
    -------
    the simplified table or None if the function fails

    """

    try:
        if isinstance(db, str): #assume this is a filename
            try:
                tab = Table.read(db)
            except Exception as e:
                print(f"Unable to open file {db}. Cannot output simple file.")
                print(traceback.format_exc())
                return -1
        elif isinstance(db,Table):
            tab = copy.copy(db)
        else:
            print("Unknown database type. Cannot output simple file.")
            return -1

        if "flag" not in tab.colnames:
            print("flag column not present. Cannot build simple flag file.")
            return -1

        flag_manual_check = "flag_manual" in tab.colnames
        tab.keep_columns(["shotid","multiframe","expnum","flag","flag_manual"])

        #combine over all exposures ... if any exposure has the flag, then the whole amp is flagged (for all exposures)
        #a flag of 0 is "bad"; all others (1 or even unset, -1) are assumed "good"
        exposures = np.unique(tab['expnum'])

        Tsimple = None

        for expnum in exposures:
            if Tsimple is None:
                Tsimple = tab[tab['expnum'] == expnum]
                Tsimple.remove_column('expnum')
                # make sure we have only 1 or 0 (anything not == 0 goes to 1)
                # here, 1 is "good"
                sn0 = np.array(Tsimple['flag'] != 0)
                Tsimple['flag'][sn0] = 1
                # now flip s|t 0 is good and -1 is bad
                Tsimple['flag'][~sn0] = -1
                Tsimple['flag'][sn0] = 0

                if flag_manual_check:
                    sn0 = np.array(Tsimple['flag_manual'] != 0)
                    Tsimple['flag_manual'][sn0] = 1
                    # now flip s|t 0 is good and -1 is bad
                    Tsimple['flag_manual'][~sn0] = -1
                    Tsimple['flag_manual'][sn0] = 0

                    # if either flag or flag_manual is set to bad, make it bad (HERE not 0 is bad)
                    Tsimple['flag'] = np.array(Tsimple['flag']) + np.array(Tsimple['flag_manual'])
                    Tsimple.remove_column('flag_manual')


            else:
                T2b = tab[tab['expnum'] == expnum]
                T2b.remove_column('expnum')
                # make sure we have only 1 or 0 (anything not == 0 goes to 1)
                # here, 1 is "good"
                sn0 = np.array(T2b['flag'] != 0)
                T2b['flag'][sn0] = 1
                # now flip s|t 0 is good and -1 is bad
                T2b['flag'][~sn0] = -1
                T2b['flag'][sn0] = 0

                if flag_manual_check:
                    sn0 = np.array(T2b['flag_manual'] != 0)
                    T2b['flag_manual'][sn0] = 1
                    # now flip s|t 0 is good and -1 is bad
                    T2b['flag_manual'][~sn0] = -1
                    T2b['flag_manual'][sn0] = 0

                    # if either flag or flag_manual is set to bad, make it bad (HERE not 0 is bad)
                    T2b['flag'] = np.array(T2b['flag']) + np.array(T2b['flag_manual'])
                    T2b.remove_column('flag_manual')
                    T2b.rename_column('flag', 'next_flag')

                # join and add up flag column
                Tsimple = join(Tsimple, T2b, keys=["shotid", "multiframe"])
                Tsimple['flag'] = np.array(Tsimple['flag']) + np.array(Tsimple['next_flag'])
                Tsimple.remove_column('next_flag')

        # lastly make < 0 become "bad" and == 0 while what was 0 becomes 1
        s = np.array(Tsimple['flag'] == 0)
        Tsimple['flag'][s] = 1  # now good
        Tsimple['flag'][~s] = 0  # now bad

        if outname is not None:
            Tsimple.write(outname+".fits",overwrite=overwrite,format="fits")
            Tsimple.write(outname+".tab",overwrite=overwrite,format='ascii')


        return Tsimple
    except Exception as e:
        print("stats_make_simple_amp_flag_file", print(traceback.format_exc()))
        return None

