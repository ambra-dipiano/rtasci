# *******************************************************************************
# Copyright (C) 2021 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# Simone Tampieri <simone.tampieri@inaf.it>
# *******************************************************************************

from genericpath import isfile
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from rtasci.aph import utils
from regions import CircleSkyRegion
from regions import write_ds9
import astropy.units as u
import numpy as np
import logging
import os
logging.basicConfig(level=logging.WARN)

# Bintable columns:
# 0  name = 'EVENT_ID'; format = '1J'; bscale = 1; bzero = 2147483648
#    name = 'TIME'; format = '1D'; unit = 's'
#    name = 'RA'; format = '1E'; unit = 'deg'
#    name = 'DEC'; format = '1E'; unit = 'deg'
#    name = 'ENERGY'; format = '1E'; unit = 'TeV'
# 5  name = 'DETX'; format = '1E'; unit = 'deg'
# 6  name = 'DETY'; format = '1E'; unit = 'deg'
#    name = 'MC_ID'; format = '1J'

class Photometrics():
    def __init__(self, args):
        self.events_data = None
        self.events_filename = None
        self.mandatory_fields = ['RA', 'DEC', 'ENERGY']
        if 'events_filename' in args:
            self.events_filename = args['events_filename']
            self.events_data = self.load_data_from_fits_file(self.events_filename)
        elif 'events_list' in args:
            self.events_data = args['events_list']
        self.events_list_checks()
        logging.info('Events data type: {}'.format(type(self.events_data)))

    def events_list_checks(self):
        """Data con be a FITS_rec or a np.recarray
        see here: https://docs.astropy.org/en/stable/io/fits/usage/table.html
        """
        if self.events_data is None:
            raise Exception('Events data is empy. Need a events list.')
        if isinstance(self.events_data, fits.fitsrec.FITS_rec):
            for f in self.mandatory_fields:
                if f not in self.events_data.columns.names:
                    raise Exception("Events data has no '{}' col".format(f))
        elif isinstance(self.events_data, np.recarray):
            for f in self.mandatory_fields:
                if f not in self.events_data.dtype.names:
                    raise Exception("Events data has no '{}' col".format(f))
        else:
            raise Exception("Events data must be FITS_rec or np.recarray")

    @staticmethod
    def load_data_from_fits_file(filename):
        """Load events extension data from a fits file.

        Parameters
        ----------
        filename: str

        Returns
        -------
        events_bintable data
        """
        with fits.open(filename, mode='readonly') as hdul:
            data = hdul['EVENTS'].data
        return data

    def region_counter(self, input_center, input_radius, emin=None, emax=None, tmin=None, tmax=None):
        """Counts photons in an input area"""
        region_center = utils.get_skycoord(input_center)
        region_radius = utils.get_angle(input_radius)

        # filtering...
        condlist = np.full(len(self.events_data.field('ENERGY')), True)
        # ... w/ energy boundaries
        if emin is not None:
            condlist &= self.events_data.field('ENERGY') >= emin
        if emax is not None:
            condlist &= self.events_data.field('ENERGY') <= emax
        # FIXME: TIME needs a better implementation
        # atm it consider users that knows the time format in the input fits
        if tmin is not None:
            condlist &= self.events_data.field('TIME') >= tmin
        if tmax is not None:
            condlist &= self.events_data.field('TIME') <= tmax

        events_list = np.extract(condlist, self.events_data)
        # events coordinates from the selected events list
        events_coords = SkyCoord(events_list.field('RA'), events_list.field('DEC'), unit='deg', frame='icrs')
        distances = region_center.separation(events_coords)
        return np.count_nonzero(distances < region_radius)

    @classmethod
    def reflected_regions(cls, input_pointing_center, input_region_center, input_region_radius):
        """Find regions with reflected algorithm.

        Parameters
        ----------
        input_pointing_center: SkyCoord or dict
        input_region_center: SkyCoord or dict
        input_region_radius: Angle or float

        Returns
        -------
        array of regions
        """
        pointing_center = utils.get_skycoord(input_pointing_center)
        region_center = utils.get_skycoord(input_region_center)
        region_radius = utils.get_angle(input_region_radius)

        # Angular separation of reflected regions. 1.05 factor is to have a margin
        region_diameter = 1.05 * 2.0 * region_radius
        radius = pointing_center.separation(region_center)
        # the numbers_of_reflected regions is the number of center that can stay
        # on the circumference, NOT the really computated number of regions.
        # the number is floor down
        numbers_of_reflected_regions = int(2 * np.pi * radius / region_diameter)
        # Indeed, we skip the source region and the two near (see below), so we
        # need at least 4 centers to get one off region.
        if numbers_of_reflected_regions < 4:
            raise Exception('the combination of region radius and coordinates does not allow to compute reflected regions.')
        regions_offset_angle = Angle(360, unit='deg') / numbers_of_reflected_regions

        regions = []
        # starting from the source region 0, we skip region 1 and region N, so 2..N-1
        starting_pos_angle = pointing_center.position_angle(region_center)
        for i in range(2, numbers_of_reflected_regions-1):
            theta = starting_pos_angle + i * regions_offset_angle
            coord_pos = pointing_center.directional_offset_by(theta, radius)
            regions.append({ 'ra': coord_pos.ra.deg, 'dec': coord_pos.dec.deg, 'rad': region_radius.deg })
        return regions

    @classmethod
    def wobble_regions(cls, *args):
        return cls.cross_regions(*args)

    @classmethod
    def cross_regions(cls, input_pointing_center, input_region_center, input_region_radius):
        """Return the three background regions starting from pointing and source one.

        Parameters
        ----------
        input_pointing_center: SkyCoord or dict
        input_region_center: SkyCoord or dict
        input_region_radius: Angle or float

        Returns
        -------
        array of regions
        """
        # FIXME Wobble algorithm has no check about distance and region radius.
        pointing_center = utils.get_skycoord(input_pointing_center)
        region_center = utils.get_skycoord(input_region_center)
        region_radius = utils.get_angle(input_region_radius)
        radius = pointing_center.separation(region_center)
        starting_pos_angle = pointing_center.position_angle(region_center)
        regions = []
        for i in range(1,4):
            theta = starting_pos_angle + i * Angle(90, unit='deg')
            coord_pos = pointing_center.directional_offset_by(theta, radius)
            regions.append({ 'ra': coord_pos.ra.deg, 'dec': coord_pos.dec.deg, 'rad': region_radius.deg })
        return regions

    @classmethod
    def write_region(cls, coords, filename, **kwargs):
        try:
            iter(coords)
        except TypeError as te:
            raise Exception('Coords must be iterable')
        circles = []
        for coord in coords:
            center = utils.get_skycoord(coord)
            rad = utils.get_angle(coord['rad'])
            circles.append(CircleSkyRegion(center=center, radius=rad, visual=kwargs))
        if isfile(filename):
            os.remove(filename)
        write_ds9(circles, filename)

