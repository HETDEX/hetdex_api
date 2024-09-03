"""

Initial Creation 2024-09-02

@author: Dustin Davis

Collection of code to support the determination of which amps, IFUs, shots are "bad"

Basic statistics with an amp+exposure as the functional unit with aggregation to the shot level.


"""

import tables
import os.path as op
import numpy as np

import astropy.stats.biweight as biweight
from astropy.table import Table
from astropy.io import ascii
import copy
import pickle
from hetdex_api.config import HDRconfig

import traceback

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
    useful for comparing maybe?


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
        dither_relflux = []
        frac_0 = []
        n_lo = []

        #         Nlo = []
        #         Scale = []
        #         Avg = []
        #         Avg_orig = []
        #         chi = []
        #         Frac_c2 = []
        #         Frac0 = []

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
                        dither_relflux.append(exp['dither_relflux'])
                        frac_0.append(exp['frac_0'])
                        n_lo.append(exp['n_lo'])


        output_tb = Table(
            [
                shotid * np.ones(len(exp_list), dtype=int),
                mf_list,
                exp_list,
                im_median,
                maskfraction,
                chi2fib_avg,
                frac_c2,
                frac_0,
                n_lo,
                sky_sub_rms,
                sky_sub_rms_rel,
                dither_relflux,


            ],
            names=["shotid", "multiframe", "expnum", "im_median", "MaskFraction", 'chi2fib_med', 'frac_c2', 'frac_0',
                   'n_lo', 'sky_sub_rms', 'sky_sub_rms_rel', 'dither_relflux']

        )

        return output_tb

    except Exception as e:
        # todo: error handling
        print("stats_shot_dict_to_table ", print(traceback.format_exc()))


###########################################
# File handling
###########################################

def stats_save_as_fits(shot_dict):
    """
    just save the shot_dict as a table (in ascii format so is easy to cat files together) named for the shotid
    """
    tab = stats_shot_dict_to_table(shot_dict)
    tab.write(shotid2datevshot(shot_dict['shotid'])+"_ampstats.dat",overwrite=True,format="ascii")



def save_shot_stats_pickle(shot_dict):
    fn = None
    try:
        # just save as a pickle
        fn = f"{shot_dict['shotid']}_stats.pickle"
        with open(fn, "wb+") as f:
            pickle.dump(shot_dict, f)

    except Exception as e:
        # todo: error handling
        print(f"save_shot_stats_pickle ({fn}) ", print(traceback.format_exc()))


def load_shot_stats_pickle(shot):
    fn = None
    try:
        fn = f"{shot}_stats.pickle"
        with open(fn, "rb") as f:
            shot_dict = pickle.load(f)
        return shot_dict
    except Exception as e:
        # todo: error handling
        print("save_shot_stats_flat ", print(traceback.format_exc()))
        return None

def stats_update_shot(h5, shot_dict):
    """
    Takes a stats dictionary and updates the associated shot

    new group, organized by multiframe + exposure, extra columns are the statistics/metrics and the interpreted results
    """

    try:
        # create the AmpStats table if it does not exist
        # Note: if exists, we want to UPDATE the rows and append if new entry
        # NOTE!!! this does NOT create new columns
        _ = h5.root.__getattr__('AmpStats')
        create_tab = False
        # print("AmpStats already exists")
    except:
        create_tab = True

        # print("Creating AmpStats")

        class AmpStats(tables.IsDescription):
            multiframe = tables.StringCol(itemsize=20, pos=0)
            expnum = tables.Int32Col(pos=1)  # could be an int8, but keep it this was to be consistent with other tables
            status = tables.Int32Col(
                pos=2)  # a status indicator, TBD ... could be a value or a bitmapped mask (-1 bad, 0 unchecked, 1 good?)
            im_median = tables.Float32Col()
            mask_fraction = tables.Float32Col()
            chi2fib_avg = tables.Float32Col()  # e.g. chi2fib_med, but is not a median, is a biweight
            frac_c2 = tables.Float32Col()

        if create_tab:
            h5.create_table(
                h5.root,
                "AmpStats",
                AmpStats,
                "Statistics for each Amp per exposure for this shot",
                expectedrows=1000,  # 78 IFU x 4 Amps x 3 Exposures = 936
            )

    astb = h5.root.AmpStats
    tab = stats_shot_dict_to_table(shot_dict)

    if create_tab:
        # all new rows
        for entry in tab:
            row = astb.row
            row['multiframe'] = entry['multiframe']
            row['expnum'] = entry['expnum']
            row['status'] = 0  # TBD
            row['im_median'] = entry['im_median']
            row['mask_fraction'] = entry['MaskFraction']
            row['chi2fib_avg'] = entry['chi2fib_med']
            row['frac_c2'] = entry['frac_c2']

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
                row['chi2fib_avg'] = entry['chi2fib_med']
                row['frac_c2'] = entry['frac_c2']

                row.append()
            elif len(row) == 1:  # update
                # need to get row iterator now, so the old (numpy) row is replaced
                # with a row from the iterator
                for row in astb.where("(multiframe==q_mf) & (expnum==q_exp)"):
                    row['status'] = 0  # TBD
                    row['im_median'] = entry['im_median']
                    row['mask_fraction'] = entry['MaskFraction']
                    row['chi2fib_avg'] = entry['chi2fib_med']
                    row['frac_c2'] = entry['frac_c2']

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


def stats_shot_dict(h5, expid=None, ifu_dict_array=None):
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
            multiframes = np.unique([mf.decode()[:-3] for mf in h5.root.Data.Images.read(field="multiframe")])
            if len(multiframes) == 0:
                status = -1
            # normally about 78x4 = 312 multiframes

            ifu_dict_array = []

            if expid is None:  # load them all
                expid_array = h5.root.Shot.read(field="expnum")
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


def stats_amp(h5, multiframe=None, expid=None, amp_dict=None):
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

                row = h5.root.Data.Images.read_where(query)

                if len(row) != 1:
                    # print("Failure, Erin's stuff", target_mf, target_expnum )
                    exp_dict['maskfraction'] = np.nan
                    exp_dict['im_median'] = np.nan
                    exp_dict['chi2fib_avg'] = np.nan
                    exp_dict['frac_c2'] = np.nan
                    exp_dict['sky_sub_rms'] = np.nan
                    exp_dict['sky_sub_rms_rel'] = np.nan  # gets added on shot rollup
                    exp_dict['n_lo'] = -1  # this is an integer, 0 is lowest possible so -1 is unset
                    exp_dict['frac_0'] = np.nan
                    # exp_dict['sky_sub_bws'] = -999.0
                    exp_dict['dither_relflux'] = np.nan
                else:

                    # todo: need to add in:
                    # background
                    # sky_sub_rms
                    # sky_sub_rms_rel
                    # N_cont
                    # nfib_bad
                    # norm

                    # similar to what original amp.dat called Scale, but over a different range
                    # Scale is from ffsky counts instead of the "image"
                    #      and for fibers 10 to 100 (inclusive), wavebins 300 to 800 (inclusive)
                    #      and explicitly excludes 0 values
                    #      and uses biweight instead of median
                    image = row["image"][0]
                    exp_dict['im_median'] = np.nanmedian(image[300:800, 300:800])

                    # similar to original amp.dat Frac0
                    # Frac0 is from ffsky counts instead of the "image"
                    #      and for fibers 10 to 100 (inclusive), wavebins 300 to 800 (inclusive)
                    #      and is the fraction of those wavelength bin values that == 0
                    error = row["error"][0]
                    sel_zero = error == 0
                    exp_dict['maskfraction'] = np.sum(sel_zero) / np.size(sel_zero)

                    tab = Table(h5.root.Data.Fibers.read_where(query))  # '(multiframe == mf) & (expnum == expi)'))

                    sel_reg = (tab['fibnum'] >= 10) & (tab['fibnum'] <= 100)
                    chi2fib = tab['chi2'][sel_reg]
                    chi2_arr = np.array(chi2fib[:, 100:1000])

                    # similar to the orginal amp.dat chi, with a +/-1 count/index difference on fibers and wavelengths
                    #   both uses fibers 10 to 100
                    #   both use the same wavelength bins (100:1000)
                    # Original chi vs this one are very close (within less than about 1%)
                    exp_dict['chi2fib_avg'] = biweight.biweight_location(chi2_arr[chi2_arr > 0])

                    # combined these are similar to the Frac_c2, with a +/-1 count/index difference on fibers and wavelengths
                    #   both uses fibers 10 to 100
                    #   both use the same wavelength bins (100:1000)
                    # Original Frac_c2 vs this one are very close (within less than about 1%)
                    exp_dict['frac_c2'] = np.sum(chi2_arr > 2) / np.sum(chi2_arr > 0)

                    #####################################
                    # similar to original amp.dat stuff
                    #####################################

                    # sky_sub_rms ~ Scale, but we use RMS, not stddev or biweight location?
                    calfib = h5.root.Data.Fibers.read_where(query, field="calfib")
                    calfib_counts = h5.root.Data.Fibers.read_where(query, field="calfib_counts")
                    calfib_ffsky = h5.root.Data.Fibers.read_where(query, field="calfib_ffsky")
                    calfib_ffsky_counts = flx2cts(calfib, calfib_counts, calfib_ffsky)

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

                    # Nlo (different M selection than frac_0)
                    # again, this is different than the original amp.dat, but may be a better calculation
                    M = calfib_ffsky_counts[:, 299:800]
                    if np.count_nonzero(M) > 0:
                        sddev = biweight.biweight_scale(M[M != 0])

                        ct = 0
                        for i in range(len(M)):
                            f = M[i][299:800]
                            if np.count_nonzero(f) > 0:
                                if biweight.biweight_location(f[f != 0]) < (-2 * sddev):
                                    ct += 1

                        exp_dict['n_lo'] = ct
                    else:
                        exp_dict['n_lo'] = -1

                    ########################################
                    # used at shot level (may have feedback)
                    ########################################
                    dither_flux = np.nan
                    try:
                        dither_flux = h5.root.Shot.read(field="relflux_virus")
                        exp_dict['dither_relflux'] = dither_flux[0][int(target_expnum) - 1]
                    except:
                        print(f"stats_amp statisitcs, dither_flux: {dither_flux}", print(traceback.format_exc()))

                    exp_dict['sky_sub_rms_rel'] = np.nan  # gets a value at shot level feedback/rollup

                    # cleanup
                    del chi2_arr
                    del tab
                    del calfib
                    del calfib_counts
                    del calfib_ffsky
                    del calfib_ffsky_counts


            except Exception as e:
                # todo: handle fail for this amp+exposure
                print("stats_amp statisitcs (Erin)", print(traceback.format_exc()))



        return amp_dict
    except Exception as e:
        # todo: error handling
        print("stats_amp statisitcs", print(traceback.format_exc()))


def stats_ifu(h5, multiframe=None, expid=None, ifu_dict=None):
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
            stats_amp(h5, multiframe=None, expid=None, amp_dict=amp_dict)

        return ifu_dict
    except Exception as e:
        # todo: error handling
        print("stats_ifu statisitcs", print(traceback.format_exc()))


def stats_shot(h5, expid=None, shot_dict=None, rollup=True):
    """
    Rollup of calls to stats_ifu(), which calls stats_amp()
    Normally 78 ifus x 4 amps x 3 expsoures = 936 total calls to stats_amp()

    """

    try:
        if shot_dict is None:
            shot_dict = stats_shot_dict(h5, expid)

        # operate on all IFUs
        for ifu_dict in shot_dict['ifu_dict_array']:
            ifu_dict = stats_ifu(h5, ifu_dict=ifu_dict)

        #
        # shot level statistics, over all amps and exposures
        #

        if rollup:
            shot_dict = stats_shot_rollup(h5, shot_dict)

        return shot_dict
    except Exception as e:
        # todo: error handling
        print("stats_shot statisitcs", print(traceback.format_exc()))


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

            # dither norms (should be the SAME VALUE for each amp for a given exposure)
            # so a mean should be redundant, but just to be safe, do it anyway
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
                        exp['sky_sub_rms_rel'] = -999.0
                        # exp['sky_sub_bws_rel'] = -999.0
                        continue
                    try:
                        exp['sky_sub_rms_rel'] = exp['sky_sub_rms'] / shot_dict['sky_sub_rms_median_exp'][sel][0]
                        # exp['sky_sub_bws_rel'] = exp['sky_sub_bws'] / shot_dict['sky_sub_bws_median_exp'][sel][0]
                    except Exception as e:
                        exp['sky_sub_rms_rel'] = -999.0
                        # exp['sky_sub_bws_rel'] = -999.0
                        print("stats_shot statisitcs, sky_sub", print(traceback.format_exc()))
        return shot_dict

    except Exception as e:
        # todo: error handling
        print("stats_shot statisitcs", print(traceback.format_exc()))


#######################################
# QC stuff
#######################################

def stats_qc_amp(amp_dict):
    """
    todo: what consitutes a bad amp ...
    evaluate and set flag in the amp_dict

    """

    #     sel1 = ((amp_stats['background'].astype(float) > -10) * (amp_stats['background'].astype(float) < 100) ) | (np.isnan(amp_stats['background']))
    #     sel2 = (amp_stats['sky_sub_rms_rel'] < 1.5) | (np.isnan(amp_stats['sky_sub_rms_rel']))
    #     sel3 = (amp_stats['sky_sub_rms'] > 0.2)  | (np.isnan(amp_stats['sky_sub_rms_rel']))
    #     sel4 = (amp_stats['im_median'] > 0.05 ) | (np.isnan(amp_stats['im_median']))
    #     sel5 = (amp_stats['MaskFraction'] < 0.25) | (np.isnan(amp_stats['MaskFraction']))
    #     sel6 = (amp_stats['N_cont'] < 35) | (np.isnan(amp_stats['N_cont']))
    #     sel7 = (amp_stats['nfib_bad'] <= 1) | (np.isnan(amp_stats['nfib_bad']))
    #     sel8 = (amp_stats['norm'] > 0.5) | (np.isnan(amp_stats['norm']))

    #     sel9_hdr3 = ((amp_stats['frac_c2'] < 0.5) | (np.isnan(amp_stats['frac_c2'])) ) * (amp_stats['date'] < 20210901)
    #     sel9_hdr4 = ((amp_stats['frac_c2'] < 0.1) | (np.isnan(amp_stats['frac_c2'])) ) * (amp_stats['date'] >= 20210901)

    #     sel9 = sel9_hdr3 | sel9_hdr4

    #     sel = sel1 & sel2 & sel3 & sel4 & sel5 & sel6 & sel7 & sel8 & sel9
    #     #sel = sel1 & sel2 & sel3 & sel6 & sel7 & sel8

    #     #sel_wiggles = amp_stats['ft_flag'] == 1
    #     sel_manual = amp_stats['flag_manual'] == 1
    #     amp_stats['flag'] = (sel * sel_manual ).astype(int)
    #     amp_stats['flag'] = (sel * sel_manual ).astype(int)

    pass


def stats_qc_ifu(ifu_dict):
    """
    todo: what consitutes a bad IFU ...
    operate on all amps in the ifu, calling into stats_qc_amp
    evaluate and set flag in the amp_dict

    """
    pass


def stats_qc_shot(ifu_dict):
    """
    todo: what consitutes a bad Shot ...
    operate on all IFUs in the shot, calling into stats_qc_ifu (which calls stats_qc_amp)
    evaluate and set flag in the amp_dict

    """
    pass