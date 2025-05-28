import numpy as np
from astropy.table import Table, hstack
from hetdex_api.survey import  FiberIndex, Survey
from multiprocessing import Pool
import time

version='5.0.1'

def get_fracs( shotid):

    fib_tab = FibIndex.return_shot( shotid)
        
    n_ifu = len( np.unique( fib_tab['ifuslot']))
    n_fib = len( fib_tab)
    fracs = dict()
    for col in fib_tab.colnames:
        if 'flag' in col:
            fracs[col] = np.float32( ( np.sum( fib_tab[col] == 1) / n_fib))
            
    return shotid, np.unique( fib_tab['datevobs'])[0], n_ifu, n_fib, fracs

def get_fracs_mp(shotid):
    try:
        res = get_fracs(shotid)
    except:
        res = shotid
    return res

FibIndex = FiberIndex('hdr5')
S = Survey('hdr5')

sel_good = S.remove_shots() # this removes bad shots, tp<0.08 frames and other frames

shotlist = S.shotid[ sel_good ]

stab = S.return_astropy_table()
stab.write( 'survey_use_{}.txt'.format(version), format='ascii', overwrite=True)

t0 = time.time()
P = Pool(20)
res = P.map( get_fracs_mp, shotlist)
P.close()
t1 = time.time()

shots = []
datevobs = []
n_ifu = []
n_fib = []
flag = []
flag_badamp = []
flag_badfib = []
flag_meteor = []
flag_satellite = []
flag_largegal = []
flag_shot = []
flag_throughput = []

for r in res:
    if isinstance(r, np.int64):
        update_shot = r
        print('updating for', r)
        r = get_fracs(update_shot)
    shots.append( r[0])
    datevobs.append( r[1])
    n_ifu.append( r[2])
    n_fib.append( r[3])
    frac = r[4]
    flag.append( frac['flag'] )
    flag_badamp.append( frac['flag_badamp'] )
    flag_badfib.append( frac['flag_badfib'] )
    flag_meteor.append( frac['flag_meteor'] )
    flag_satellite.append( frac['flag_satellite'] )  
    flag_largegal.append( frac['flag_largegal'] )
    flag_shot.append( frac['flag_shot'] )
    flag_throughput.append( frac['flag_throughput'] )
    
cols = ['shotid',
        'datevobs',
        'n_ifu',
        'n_fib',
        'frac', 
        'frac_badamp',
        'frac_badfib',
        'frac_meteor',
        'frac_satellite',
        'frac_largegal',
        'frac_shot',
        'frac_throughput']

tab = Table(
    [
        shots,
        datevobs,
        n_ifu,
        n_fib,
        flag,
        flag_badamp,
        flag_badfib,
        flag_meteor,
        flag_satellite,
        flag_largegal,
        flag_shot,
        flag_throughput
    ],
    names = cols,
    dtype=[ np.int64,
            str,
            np.int32,
            np.int32,
            np.float32,
            np.float32,
            np.float32,
            np.float32,
            np.float32,
            np.float32,
            np.float32,
            np.float32,
    ]
)
tab.write( 'survey_frac_{}.txt'.format(version), format='ascii', overwrite=True)

FibIndex.close()
S.close()

print('Done in {:3.2f} min.'.format( (t1-t0)/60))
