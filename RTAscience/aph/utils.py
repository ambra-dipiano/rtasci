# Copyright 2019,2020 Simone Tampieri
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import csv
import math
from astropy.coordinates import SkyCoord, Angle
from RTAscience.aph.photometry import Photometrics
from RTAscience.aph.irf import EffectiveArea

def photometrics_counts(events_list, events_type, pointing, true_coords, region_rad=0.2):
  phm = Photometrics({events_type: events_list})
  reflected_regions = phm.reflected_regions(pointing, true_coords, region_rad)
  on_count = phm.region_counter(true_coords, region_rad)
  off_count = 0
  for r in reflected_regions:
    off_count += phm.region_counter(r, r['rad'])
  alpha = 1 / len(reflected_regions)
  return {'on': on_count, 'off': off_count, 'alpha': alpha, 'excess': on_count - alpha * off_count}

def li_ma (n_on, n_off, alpha):
    if n_on <= 0 or n_off <= 0 or alpha == 0:
        return None
    fc = (1 + alpha) / alpha
    fb = n_on / (n_on + n_off)
    f  = fc * fb
    gc = 1 + alpha
    gb = n_off / (n_on + n_off)
    g  = gc * gb
    first  = n_on * math.log(f)
    second = n_off * math.log(g)
    fullb   = first + second
    return math.sqrt(2) * math.sqrt(fullb)

def read_timeslices_tsv(filename):
    ts = []
    with open(filename, mode='r', newline='\n') as fh:
        reader = csv.reader(fh, delimiter='\t')
        headers = next(reader)
        for row in reader:
            ts.append(dict(zip(headers, row)))
    return ts

def get_angle(input_angle):
    ang = None
    if isinstance(input_angle, Angle):
        ang = input_angle
    elif isinstance(input_angle, float):
        ang = Angle(input_angle, unit='deg')
    else:
        raise Exception('The input parameter must be an Angle or a float for decimal degree.')
    return ang

def get_skycoord(pnt_coord):
    coord = None
    if isinstance(pnt_coord, SkyCoord):
        coord = pnt_coord
    elif isinstance(pnt_coord, tuple):
        coord = SkyCoord(ra=pnt_coord[0], dec=pnt_coord[1], unit='deg', frame='icrs')
    elif isinstance(pnt_coord, dict) and 'ra' in pnt_coord and 'dec' in pnt_coord:
        coord = SkyCoord(ra=pnt_coord['ra'], dec=pnt_coord['dec'], unit='deg', frame='icrs')
    else:
        raise Exception('The input parameter must be a SkyCoord, a { "ra": 12.3, "dec": 45.6 } dictionary or (12.3, 45.6) tuple.')
    return coord

def counting(phm, src, rad, off_regions, e_min=None, e_max=None, t_min=None, t_max=None, draconian=False):
    on_count = phm.region_counter(src, rad, emin=e_min, emax=e_max, tmin=t_min, tmax=t_max)
    off_count = 0
    for r in off_regions:
        off_count += phm.region_counter(r, r['rad'], emin=e_min, emax=e_max, tmin=t_min, tmax=t_max)

    alpha = 1 / len(off_regions)
    excess = on_count - alpha * off_count
    signif = li_ma(on_count, off_count, alpha)

    # !!! here we can implement checks
    err_note = None
    if on_count < 10 or off_count < 10:
        err_note = 'Not compliant with Li & Ma requirements.'
        if draconian:
            raise Exception(err_note)

    return on_count, off_count, alpha, excess, signif, err_note

def find_off_regions(phm, algo, src, pnt, rad, verbose=False, save=None):
    if not len(pnt) == 2 and len(src) == 2:
        raise Exception('need source and pointing coordinates and a region radius to do aperture photometry')

    off_regions = None
    if algo.lower() == 'cross':
        off_regions = phm.cross_regions(pnt, src, rad)
    elif algo.lower() == 'reflection':
        off_regions = phm.reflected_regions(pnt, src, rad)
    else:
        raise Exception('invalid background regions algorithm')

    if verbose > 1:
        print('off regions algorithm:', algo)
        for i, o in enumerate(off_regions):
            print('      off regions #{:02d}:'.format(i), o)

    if save:
        phm.write_region(off_regions, save, color='red', dash=True, width=2)

    return off_regions

def aeff_eval(args, src, pnt):
    if not(args.energy_min and args.energy_max and args.pixel_size and args.power_law_index):
        raise Exception('need energy min and max, a pixel size to eval the flux')

    aeff = EffectiveArea(irf_filename=args.irf_file)
    # these IRFs return value in m², so we need convert
    # the source data struct need a 'rad'
    source_reg_aeff = aeff.weighted_value_for_region(src, pnt, [args.energy_min, args.energy_max], args.pixel_size, args.power_law_index) * 1e4 # cm2
    return source_reg_aeff

class ObjectConfig(object):
    def __init__(self, d):
        self.__dict__ = d
        self.__dict__['begin_time'] = float(self.__dict__['begin_time'])
        self.__dict__['end_time'] = float(self.__dict__['end_time'])
        self.__dict__['source_ra'] = float(self.__dict__['source_ra'])
        self.__dict__['source_dec'] = float(self.__dict__['source_dec'])
        self.__dict__['region_radius'] = float(self.__dict__['region_radius'])
        self.__dict__['verbose'] = float(self.__dict__['verbose'])
        self.__dict__['energy_min'] = float(self.__dict__['energy_min'])
        self.__dict__['energy_max'] = float(self.__dict__['energy_max'])
        self.__dict__['pixel_size'] = float(self.__dict__['pixel_size'])
        self.__dict__['power_law_index'] = float(self.__dict__['power_law_index'])
        #self.__dict__['irf_file'] = float(self.__dict__['irf_file'])
