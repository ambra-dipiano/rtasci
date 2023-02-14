# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import argparse
from rtasci.lib.RTAIrfs import RTAIrfs

# files ---!
parser = argparse.ArgumentParser()
parser.add_argument("--caldb", nargs="+", default=['prod3b', 'prod3b-v2'])
args = parser.parse_args()

# initialise ---!
for db in args.caldb:
  print(f'processing {db} degradation')
  irf = os.listdir(os.environ.get('CTOOLS') + '/share/caldb/data/cta/' + db + '/bcf/')
  for fits in irf:
    irfObj = RTAIrfs()
    irfObj.irf = fits
    irfObj.caldb = db
    irfObj.degradeIrf()

print('caldb degradation completed')