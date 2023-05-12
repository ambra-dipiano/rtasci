# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import os
import sys

start = int(sys.argv[1])
if len(sys.argv) > 2:
    runs = int(sys.argv[2])
else:
    runs = 20

for i in range(runs):
    os.system(f'scancel {start+i}')


