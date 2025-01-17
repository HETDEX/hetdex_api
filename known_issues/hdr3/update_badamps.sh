#!/usr/bin/env python

#
# ampstats.fits is the file of record (but only for HDR5+)
# amp_flag.fits / .tab are the simplified versions
#

from astropy.table import Table
from h5tools import amp_stats as AS
from hetdex_api.config import HDRconfig
import os.path as op
from tqdm import tqdm
import sys
import numpy as np
import hetdex_api

badamps_path = op.join("/".join(hetdex_api.__path__[0].split("/")[:-1]),"known_issues/hdr3/") #note all hdr(x) point to hdr3

bad_flag_value = 0


recfn0 = "amp_flag"#.fits"
recfn0_full = "ampstats"#.fits"

badfn1 = "badamps.list"
badfn2 = "badamps_single.list"


T0 = Table.read(f"{recfn0}.fits")
T0_full = Table.read(f"{recfn0_full}.fits")

if op.isfile(badfn1):
    print(f"Using local {badfn1}")
    T1 = Table.read(badfn1,format="ascii")
else:
    print(f"Using config {op.join(badamps_path,badfn1)}")
    T1 = Table.read(op.join(badamps_path,badfn1),format="ascii")
if op.isfile(badfn2):
    print(f"Using local {badfn2}")
    T2 = Table.read(badfn2,format="ascii")
else:
    print(f"Using config {op.join(badamps_path,badfn2)}")
    T2 = Table.read(op.join(badamps_path,badfn2),format="ascii")

#add date column to recfn0 so can match against the multi (this is slow, but straightforward)
print("adding proxy columns ...")
T0['date'] = [int(str(s)[:-3]) for s in T0['shotid']]


# first the ranges (T1)
print("badamps.list")
for row in tqdm(T1):
    #will need to update both the full and simple versions

    date_start = row['date_start']
    date_end = row['date_end']
    mf = row["multiframe"].encode()

    #check the read in start shot, stop shot, multiframe
    #simple version
    sel = np.array(T0['multiframe']==mf) & np.array(T0['date'] >= date_start) & np.array(T0['date'] <= date_end)
    if np.count_nonzero(sel) > 0:
        #update
        T0['flag'][sel] = bad_flag_value
    else:
        print(f"Not found: {str(row['ifuslot']).zfill(3)} {row['ampid']} {mf.decode()} {date_start} {date_end}")
        #print(row)


    #now the full version
    sel = np.array(T0_full['multiframe']==mf) & np.array(T0_full['date'] >= date_start) & np.array(T0_full['date'] < date_end)

    if np.count_nonzero(sel) > 0:
        #update
        T0_full['flag_manual'][sel] = bad_flag_value
        T0_full['flag_manual_desc'][sel] = b'badamps.list'


#now the singles
print("badamps_single.list")
for row in tqdm(T2):
    #will need to update both the full and simple versions??
    shotid = row['shotid']
    mf = row["multiframe"].encode()

    #check the read in start shot, stop shot, multiframe
    #simple version
    sel = np.array(T0['multiframe']==mf) & np.array(T0['shotid'] == shotid)
    if np.count_nonzero(sel) > 0:
        #update
        T0['flag'][sel] = bad_flag_value
    else:
        print(f"Not found: {shotid} {mf}")

    #now the full version
    sel = np.array(T0_full['multiframe']==mf) & np.array(T0_full['shotid'] == shotid)

    if np.count_nonzero(sel) > 0:
        #update
        T0_full['flag_manual'][sel] = bad_flag_value
        T0_full['flag_manual_desc'][sel] = b'badamps.list'

T0.remove_column('date')
T0.write(f"{recfn0}.updated.fits",overwrite=True,format="fits")
T0.write(f"{recfn0}.updated.tab",overwrite=True,format="ascii")

#fits only format
T0_full.write(f"{recfn0_full}.updated.fits",overwrite=True,format="fits")


print("If you are happy with the update, replace amp_flags.fits, amp_flags.tab, and ampstats.fits")
print("with amp_flags.updated.fits, amp_flags.updated.tab, and ampstats.updated.fits, respectively.")
