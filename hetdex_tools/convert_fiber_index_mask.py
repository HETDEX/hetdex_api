import numpy as np
from astropy.table import Table, hstack
from hetdex_api.survey import  FiberIndex, Survey
from multiprocessing import Pool
import time


def make_ascii_fiber_table( shotid):

    try:
        fib_tab = FibIndex.return_shot( shotid)
        
        multiname = [ "{}_{}".format( row['multiframe'], str( row['fibnum']).zfill(3)) for row in fib_tab]
        
        fib_tab.add_column( multiname, name='multiname', index=0)

        exp = [ "exp{}".format( str( row['expnum']).zfill(2)) for row in fib_tab]
    
        fib_tab.add_column( exp, name='exp', index=1)
    
        for col in fib_tab.colnames:        
            if 'flag' in col:
                fib_tab[col].dtype = np.int8
    
        fib_tab.write('fiber_index_mask/{}/fibmask_{}_{}.txt'.format(version, fib_tab['datevobs'][0], version), format='ascii', overwrite=True)

    except:
        print('Failed for {}'.format(shotid))
        
    return

FibIndex = FiberIndex('hdr5')
S = Survey('hdr5')

version='5.0.0'

sel_good = S.remove_shots()
shotlist = S.shotid[sel_good]

t0 = time.time()
P = Pool(32)
res = P.map( make_ascii_fiber_table, shotlist)
P.close()
t1 = time.time()

FibIndex.close()

print('Done in {:3.2f} min.'.format( (t1-t0)/60))
