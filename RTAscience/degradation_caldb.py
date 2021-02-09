# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

from RTAscience.lib.RTAIrfs import RTAIrfs
import os

# files ---!
caldb = ['prod3b', 'prod3b-v2']

# initialise ---!
for db in caldb:
  print(f'processing {db} degradation')
  irf = os.listdir(os.environ.get('CTOOLS') + '/share/caldb/data/cta/' + db + '/bcf/')
  for fits in irf:
    irfObj = RTAIrfs()
    irfObj.irf = fits
    irfObj.caldb = db
    irfObj.degradeIrf()

print('caldb degradation completed')