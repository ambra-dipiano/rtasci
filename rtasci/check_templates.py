# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import argparse
import numpy as np
import matplotlib.pyplot as plt
from rtasci.lib.RTAVisualise import get_template_lc, get_template_spectra

parser = argparse.ArgumentParser(description='Simulate empty fields.')
parser.add_argument('--runid', type=str, required=True, help="template RUNID")
parser.add_argument('-p', '--path', type=str, default='/Users/iasfbo/Desktop/CTA/DATA/templates/grb_afterglow/GammaCatalogV1.0/', help='absolute path where cat is installed')
args = parser.parse_args()

time, flux = get_template_lc(runid=args.runid+'.fits', erange=[0.04, 1], path=args.path)
print(np.shape(time), np.shape(flux))

plt.plot(time, flux)
plt.xlabel('time')
plt.ylabel('flux')
plt.yscale('log')
plt.xscale('log')
plt.savefig(f'test_lightcurve_{args.runid}.png')
plt.close()

energy, flux = get_template_spectra(runid=args.runid+'.fits', erange=[0.04, 1], path=args.path)
print(np.shape(energy), np.shape(flux))

plt.plot(energy, flux)
plt.xlabel('energy')
plt.ylabel('flux')
plt.yscale('log')
plt.xscale('log')
plt.savefig(f'test_spectra_{args.runid}.png')
plt.close()