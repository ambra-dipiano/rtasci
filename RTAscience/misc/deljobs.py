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

start = float(sys.argv[1])

runs = 20
for i in range(runs):
    os.system(f'scancel {start+i}')