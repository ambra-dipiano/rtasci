# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Leonardo Baroncelli <leonardo.baroncelli@inaf.it>
# *******************************************************************************

import os
import argparse
from os.path import isdir, join
import numpy as np
from rtasci.lib.RTAStats import ts_wilks, p_values, ts_wilks_cumulative

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", nargs='+', type=str, required=True)
    parser.add_argument("-o", "--outname", type=str, required=True)
    args = parser.parse_args()

    outDir="./tmp"
    if not isdir(outDir):
        os.mkdir(outDir)
    data = [np.genfromtxt(file, usecols=(0), skip_header=0, dtype=float) for file in args.files]

    ts_wilks(data, nbin=100, width=None, filename = join(outDir, f"{args.outname}_wilks.png"), overlay=False, write_data=True)
    ts_wilks_cumulative(data, nbin=100, width=None, filename = join(outDir, f"{args.outname}_cumulative.png"), overlay=False, write_data=True)
    p_values(data, nbin=100, width=None, filename = join(outDir, f"{args.outname}_pvalues.png"), overlay=False, sigma5=False, write_data=True)
