# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

from astropy.coordinates import SkyCoord
from gammapy.analysis import AnalysisConfig
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel, PointSpatialModel


def gammapy_config(cfg, obs, target=None, pointing=None, radius=0.2, rbins=20, etrue=[0.02, 200], tbins=30, maxoffset=2.5, fitflux=False, fbins=30, source='GRB', level='info', stack=False, exclusion=None, safe_mask=['aeff-default', 'offset-max'], save=False, blind=False):
    """Wrapper for gammapy configuration class."""
    if target is None and pointing is None:
        raise ValueError("Either target or pointing must be not null.")
    config = AnalysisConfig()
    config.general.log = {'level': level}
    config.datasets.type = cfg.get('type')
    config.datasets.stack = stack
    config.observations.datastore = ''
    config.observations.obs_file = obs
    # 3d analysis or blind
    if cfg.get('type').lower() == "3d" or cfg.get('blind') == True:
        #print(f"BLIND {int(cfg.get('roi')*2*cfg.get('skyroifrac'))}")
        # type
        config.datasets.type = "3d"
        # geometry of the map for 3d
        config.datasets.geom.wcs.skydir = {'lon': pointing.ra, 'lat': pointing.dec, 'frame': 'icrs'}  
        config.datasets.geom.wcs.fov = {'width': f"{int(cfg.get('roi')*2*cfg.get('skyroifrac'))} deg", 'height': f"{int(cfg.get('roi')*2*cfg.get('skyroifrac'))} deg"}
        config.datasets.geom.wcs.binsize = f"{cfg.get('skypix')} deg"
        if not cfg.get('blind'):
            config.datasets.background=dict(method="fov_background", exclusion=exclusion)
    # 1d analysis
    elif cfg.get('type').lower() == "1d":
        # type
        config.datasets.type = "1d"
        # ON region and make sure that PSF leakage is corrected
        config.datasets.on_region = dict(frame="icrs", lon=f"{target[0]} deg", lat=f"{target[1]} deg", radius=f"{radius} deg")
        # background
        config.datasets.background=dict(method="reflected", exclusion=exclusion)#, parameters={'max_region_numberint': 23})
    # what maps to compute
    config.datasets.map_selection = ['counts', 'exposure', 'background', 'psf', 'edisp']
    # safe mask and IRF
    config.datasets.safe_mask.methods = safe_mask
    config.datasets.containment_correction = True
    # roi
    config.datasets.geom.selection.offset_max = f"{maxoffset} deg"
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

def set_model(target, source='Crab', freeze_spc=['index'], freeze_spt=['lon_0', 'lat_0'], default=True, index=2.4):
    """Wrapper to create gammapy model."""
    target = SkyCoord(target[0], target[1], unit='deg', frame='icrs')
    spatial_model = PointSpatialModel(lon_0=target.ra, lat_0=target.dec, frame="icrs")
    if default:
        spectral_model = PowerLawSpectralModel(index=index, amplitude="5.7e-16 cm-2 s-1 MeV-1", reference="1e6 MeV")
    else:
        raise ValueError('Option default=False not implemented yet.')
    for prm in freeze_spt:
        spatial_model.parameters[prm].frozen = True
    for prm in freeze_spc:
        spectral_model.parameters[prm].frozen = True
    sky_model = SkyModel(spatial_model=spatial_model, spectral_model=spectral_model, name=source)
    return sky_model, spectral_model, spatial_model