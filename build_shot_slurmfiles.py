#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 11:52:01 2019

@author: gregz
"""

import os.path as op
import numpy as np
import config
import subprocess

#filename = sys.argv[1]
filename = op.join(config.red_dir, 'throughput/hdr1.scilist')
object_table = [line.rstrip('\n').split() for line in open(filename)]
object_table = [[_object[0], _object[1]] for _object in object_table
                if _object[0][:4] == '2019']

N = len(object_table) / 20 + 1
object_chunks = np.array_split(object_table, N)
spath = op.join(config.software_dir, 'scripting', 'rwrangler_shotlists.slurm')
for i, object_chunk in enumerate(object_chunks):
    rname = 'rwrangler_shotlists_%i' % (i+1)
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
    sedcall = ('sed "s/rwrangler_shotlists/%s/g" '
              '%s > %s' % (rname, spath, sname))
    proc = subprocess.Popen(sedcall.split())