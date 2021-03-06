###!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 11:52:01 2019

@author: gregz
"""

import os.path as op
import numpy as np
import subprocess
import datetime

from hetdex_api import config

#filename = sys.argv[1]
filename = op.join(config.red_dir, 'hdr1.scilist')
object_table = [line.rstrip('\n').split() for line in open(filename)]

filename = op.join(config.red_dir, 'hdr1.callist')
object_table2 = [line.rstrip('\n').split() for line in open(filename)]
object_table = object_table + object_table2

N = len(object_table) / 20 + 1
object_chunks = np.array_split(object_table, N)
spath = op.join(config.software_dir, 'scripting', 'rwrangler_shotfiles.slurm')
G = open(op.join(config.software_dir, 'scripting/calls_to_run/at_calls.txt'),
         'w')
atcalls = []
D = datetime.datetime.now()
for i, object_chunk in enumerate(object_chunks):
    rname = 'rwrangler_shotfiles_%i' % (i+1)
    name = op.join(config.software_dir, 'scripting', 'calls_to_run',
                   rname)
    sname = name + '.slurm'
    f = open(name, 'w')
    s = []
    d = op.join(config.software_dir, 'HETDEX_API/mkhdf5')
    for _object in object_chunk:
        line = ('%s %s %s /work/03946/hetdex/hdr1/reduction/data' %
                (d, _object[0], _object[1]))
        s.append(line)
    f.write('\n'.join(s))
    f.close()
    sedcall = ('sed "s/rwrangler_shotfiles/%s/g" '
              '%s > %s' % (rname, spath, sname))
    subprocess.call([sedcall], shell=True)
    d1 = D + datetime.timedelta(0, 600.*(i+1), 0)
    d2 = d1.strftime('%H:%M %B %d')
    atcalls.append('echo "source ~hetdex/.bashrc; sbatch %s" | at %s' % (sname, d2))
G.write('\n'.join(atcalls))
