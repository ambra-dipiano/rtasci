# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import sys
from os.path import expandvars
from rtasci.lib.RTACtoolsSimulation import RTACtoolsSimulation as sim

events = sim()
events.input = expandvars(sys.argv[1])
events.t = [0, 1000]
events.sortObsEvents(key='TIME')

