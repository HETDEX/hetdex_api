import sys
import os.path as op
import tables as tb

from hetdex_api.config import HDRconfig

config = HDRconfig('hdr2.1')

date = sys.argv[1]
obs = sys.argv[2]

shotfile1 = op.join( config.data_dir, "{}v{}.h5".format(date, str(obs).zfill(3)))
fileh1 = tb.open_file(shotfile1)

shotfile2 = "{}v{}.h5".format(date, str(obs).zfill(3))
fileh2 = tb.open_file(shotfile2, 'a')

print('Copying Astrometry group from hdr2.1 for {}'.format(shotfile2))

fileh1.copy_node('/Astrometry', newparent=fileh2.root, recursive=True, overwrite=True)

fileh1.close()
fileh2.close()
