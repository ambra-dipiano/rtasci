# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import sys
import os
import pandas as pd
import argparse

# ---------------------------------------------------------------------------- !

parser = argparse.ArgumentParser()
parser.add_argument('--tt', type=float, default=1000000, help='total trials')
parser.add_argument('--tn', type=float, default=50000, help='trials per node')
parser.add_argument('--delay', type=float, default=50, help='delay')
parser.add_argument('--off', type=str, default='gw', help='offset')
parser.add_argument('--flux', type=float, default=1, help='flux scaling factor')
parser.add_argument('--env', type=str, default='ctools', help='environment to activate')
parser.add_argument('--pipe', type=str, default='pipeline', help='which pipeline to run')
args = parser.parse_args()

path = '/data01/homes/cta/gammapy_integration/cta-sag-sci/RTAscience'
for i in range(int(args.tt/args.tn)):
    outname = f'{path}/jobs/seed{i*args.tn+1:06d}-{(i+1)*args.tn:06d}.sh'
    # write bash
    sh = outname.replace('cfg/', 'jobs/') + '.sh'
    with open(sh, 'w+') as f:
        f. write('#!/bin/bash\n')
        f.write(f'\nsource activate {args.env}')
        f.write('\n\texport DATA=/data01/homes/cta/gammapy_integration/DATA/')
        if args.pipe == 'pipeline':
            f.write(f'\n\tpython {path}/thesis/pipeline.py --count {i*args.tn} --trials {args.tn} --delay {args.delay} --flux {args.flux} --off {args.off}\n')
        elif args.pipe == 'wilks':
            f.write(f'\n\tpython {path}/thesis/wilks.py --count {i*args.tn} --trials {args.tn} --off {args.off}\n')           
        else:
            raise ValueError('Argument "pipe" is not valid.')
    os.system(f'sbatch {sh}')