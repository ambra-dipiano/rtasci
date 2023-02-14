# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import numpy as np
from rtasci.lib.RTAStats import *
from os.path import expandvars, join

x = np.random.chisquare(1, int(1e7))
#np.insert(x, 0, 0)

path = expandvars('$DATA/outputs/LEO')
save_data_on_file(x, filename=join(path, 'chi2sample_s1e7_r0-36_n100.csv'), hdr=None)

