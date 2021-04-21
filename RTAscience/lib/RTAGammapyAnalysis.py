# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.data import EventList, GTI, Observation, Observations
from gammapy.irf import load_cta_irfs
from gammapy.modeling import Fit
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel, PointSpatialModel


def gammapy_config(cfg, target, radius=0.2, rbins=20, etrue=[0.02, 300], tbins=30, bkg_method='reflected', maxoffset=2.5, fitflux=False, fbins=30, source='GRB', level='info', stack=False, exclusion=None, save=False):
    config = AnalysisConfig()
    config.general.log = {'level': level}
    config.datasets.type = cfg.get('type')
    config.datasets.stack = stack
    # 1d analysis
    if cfg.get('type').lower() == "1d":
        # ON region and make sure that PSF leakage is corrected
        config.datasets.on_region = dict(frame="icrs", lon=f"{target[0]} deg", lat=f"{target[1]} deg", radius=f"{radius} deg")
        config.datasets.geom.selection.offset_max = f"{maxoffset} deg"
    config.datasets.containment_correction = True
    # background
    config.datasets.background=dict(method=bkg_method, exclusion=exclusion)
    #config.datasets.safe_mask.methods = ["edisp-bias"]
    # energy binning for the spectra
    config.datasets.geom.axes.energy = dict(min=f"{cfg.get('emin')} TeV", max=f"{cfg.get('emax')} TeV", nbins=rbins)
    config.datasets.geom.axes.energy_true = dict(min=f"{etrue[0]} TeV", max=f"{etrue[1]} TeV", nbins=tbins)
    # fit and flux points (currently not used by the analysis)
    if fitflux:
        config.fit.fit_range = dict(min=f"{cfg.get('emin')} TeV", max=f"{cfg.get('emax')} TeV")
        config.flux_points.energy = dict(min=f"{cfg.get('emin')} TeV", max=f"{cfg.get('emax')} TeV", nbins=fbins)
        config.flux_points.source = source
    # write
    if save:
        config.write("config1d.yaml", overwrite=True)
    return config

def set_model(target, source='Crab', freeze_spc=['index'], freeze_spt=['lon_0', 'lat_0'], default=True):
    target = SkyCoord(target[0], target[1], unit='deg', frame='icrs')
    spatial_model = PointSpatialModel(lon_0=target.ra, lat_0=target.dec, frame="icrs")
    if default:
        spectral_model = PowerLawSpectralModel(index=2.4, amplitude="5.7e-16 cm-2 s-1 MeV-1", reference="1e6 MeV")
    else:
        raise ValueError('Option default=False not implemented yet.')
    for prm in freeze_spt:
        spatial_model.parameters[prm].frozen = True
    for prm in freeze_spc:
        spectral_model.parameters[prm].frozen = True
    sky_model = SkyModel(spatial_model=spatial_model, spectral_model=spectral_model, name=source)
    return sky_model, spectral_model, spatial_model