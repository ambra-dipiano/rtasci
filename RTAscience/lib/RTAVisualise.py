# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import matplotlib.pyplot as plt
import pyregion
import seaborn as sns
from astropy.io import fits
from matplotlib.colors import SymLogNorm
from matplotlib import rc
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Rectangle
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.visualization as avis


# handle DS9 regiond (wip) ---!
def handleReg(reg, col='black'):
  r = pyregion.open(reg)
  r[0].attr[1]['color'] = col
  return r

# plot sky map ---!
def showSkymap(file, reg='none', col='black', suffix='none', title='skymap',
               xlabel='R.A. (deg)', ylabel='Dec (deg)', fontsize=12, show=True, tex=False):
  with fits.open(file) as hdul:
    wcs = WCS(hdul[0].header)
    data = hdul[0].data
    hdr = hdul[0].header

  plt.rc('text', usetex=False) if tex else None
  ax = plt.subplot(111)
  # load region ---!
  if reg != 'none' :
    r = pyregion.open(reg).as_imagecoord(hdr)
    for i in range(len(r)):
      r[i].attr[1]['color'] = col
      patch_list, text_list = r.get_mpl_patches_texts()
      for p in patch_list:
        ax.add_patch(p)
      for t in text_list:
        ax.add_artist(t)
  # plot with projection ---!
  plt.subplot(projection=wcs)
  plt.imshow(data, cmap='jet', norm=SymLogNorm(1), interpolation='gaussian')
  plt.grid(color='white', ls='solid')
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title, fontsize=fontsize)
  plt.colorbar().set_label('cts')
  # save fig ---!
  if suffix != 'none' :
    plt.savefig(file.replace('.fits', '_%s.png' % suffix))
  else :
    plt.savefig(file.replace('.fits', '.png'))
  # show fig ---!
  plt.show() if show else None
  plt.close()
  return

# plot residual count map ---!
def showResmap(file, reg='none', col='black', suffix='none', title='map redisuals',
               xlabel='R.A. (deg)', ylabel='Dec (deg)', fontsize=12, show=True, tex=False):
  with fits.open(file) as hdul:
    data = hdul[0].data
    hdr = hdul[0].header

  plt.rc('text', usetex=False) if tex else None
  ax = plt.subplot(111)
  # load region ---!
  if reg != 'none':
    r = pyregion.open(reg).as_imagecoord(hdr)
    for i in range(len(r)):
      r[i].attr[1]['color'] = col
      patch_list, text_list = r.get_mpl_patches_texts()
      for p in patch_list:
        ax.add_patch(p)
      for t in text_list:
        ax.add_artist(t)
  # plot ---!
  plt.imshow(data, cmap='bwr', interpolation='gaussian', filterrad=0.04)
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title, fontsize=fontsize)
  plt.colorbar().set_label('Significance $\sigma$')
  # save fig ---!
  if suffix != 'none':
    plt.savefig(file.replace('.fits', '_%s.png' % suffix))
  else:
    plt.savefig(file.replace('.fits', '.png'))
  #show fig ---!
  plt.show() if show else None
  plt.close()
  return

# show spectral residuals ---!
def showResiduals(file, yscale='log', title='spectral residuals', figsize=(10,8),
                  xlabel='energy (TeV)', ylabel=('counts', 'residuals'), fontsize=12, show=True, tex=False):
  with fits.open(file) as hdul:
    data = hdul[1].data
    # store data ---!
    Emin = data.field(0)
    Emax = data.field(1)
    cts = data.field(2)
    model = data.field(3)
    res = data.field(4)
  # binning ---!
  en_bins = Emax - 0.5 * (Emax - Emin)

  plt.figure(figsize=figsize)
  plt.rc('text', usetex=False) if tex else None
  if yscale.lower() == 'lin' :
    ax1 = plt.subplot(211, xscale='log')
  else :
    ax1 = plt.subplot(211, yscale='log', xscale='log')
  # plot ---!
  plt.plot(en_bins, cts, 'ro')
  plt.step(Emax, cts, 'r-', where='pre', label='cts')
  plt.step(Emin, model, 'g-', where='post', label='model')
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel[0], fontsize=fontsize)
  plt.title(title, fontsize=fontsize)
  plt.legend(loc=0)
  plt.grid(True)
  # plot ---!
  plt.subplot(212, sharex=ax1)
  plt.plot(en_bins, res, 'b+', label='cts residuals')
  plt.axhline(y=0, c='r', lw='1', ls='--')
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel[1], fontsize=fontsize)
  plt.grid(True)
  # adjust ---!
  plt.subplots_adjust(hspace=0.2)
  # save fig ---!
  plt.savefig(file.replace('.fits', '.png'))
  # show fig ---!
  plt.show() if show == True else None
  plt.close()
  return

# flux butterfly diagram ----!
def showButterfly(file, flux_pnts=0.0, fluxEn_pnts=0.0, suffix='none', title='flux', fontsize=12,
                  xlabel='Energy (TeV)', ylabel='E $\cdot \\frac{dN}{dE}$ (erg/$cm^2$/s)', show=True, tex=False):
  data = np.loadtxt(file, delimiter=' ')
  energy = data[:, 0]
  intensity = data[:, 1]
  lerr = data[:, 2]
  uerr = data[:, 3]

  # from intensity (ph/cm^2/s/MeV) to flux in (erg/cm^2/s) ---!
  (flux, f_lerr, f_uerr) = ((intensity, lerr, uerr) * (energy ** 2)) / (6.4215 * 1e5)
  if flux_pnts != 0.0  and fluxEn_pnts != 0.0:
    f_pnts = (np.array(flux_pnts) * (np.array(fluxEn_pnts) ** 2)) / (6.4215 * 1e5)
  else:
    print('without flux points')

  plt.rc('text', usetex=False) if tex else None
  ax = plt.subplot(111, yscale='log', xscale='log')
  # plot ---!
  plt.plot(energy / 1e6, flux, 'b-', alpha=1, label='best fit')
  plt.fill_between(energy / 1e6, f_lerr, f_uerr, facecolor='blue', alpha=0.3, label='errors')
  if flux_pnts != 0.0 and fluxEn_pnts != 0.0:
    plt.scatter(np.array(fluxEn_pnts) / 1e6, f_pnts, marker='o', color='red', label='data flux pnts')
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title, fontsize=fontsize)
  plt.legend(loc=0)
  plt.grid(True)
  # save fig ---!
  if suffix != 'none':
    plt.savefig(file.replace('.txt', '_%s.png' % suffix))
  else:
    plt.savefig(file.replace('.txt', '.png'))
  # show fig ---!
  plt.show() if show == True else None
  plt.close()
  return

# plot spectrum ---!
def showSpectrum(file, figsize=(8,15), fontsize=12, title=('spectrum with errors', 'spectrum with errors', 'log spectrum'),
                 xlabel='Energy (TeV)', ylabel='Flux (erg/$cm^2$/s)', show=True, tex=False):
  with fits.open(file) as hdul:
    data = hdul[1].data
    # store data ---!
    en = data.field(0)
    en_down = data.field(1)
    en_up = data.field(2)
    flux = data.field(3)
    err_flux = data.field(4)
  # adjusting errors ---!
  en_err = (en_up - en_down) / 2

  plt.rc('text', usetex=False) if tex else None
  fig = plt.figure(figsize=figsize)
  ax1 = plt.subplot(311, xscale='log')
  # plot ---!
  plt.errorbar(en, flux, yerr=err_flux, xerr=en_err, fmt='ro', label='data')
  plt.step(en, flux, 'r-', where='mid')
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title[0], fontsize=fontsize)
  plt.legend(loc=0)
  plt.grid(True)

  ax2 = plt.subplot(312, sharex=ax1)
  plt.plot(en, flux, 'ro', label='data')
  plt.step(en, flux, 'r-', where='mid')
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title[1], fontsize=fontsize)
  plt.legend(loc=0)
  plt.grid(True)

  ax3 = plt.subplot(313, sharex=ax1, yscale='log')
  plt.plot(en, flux, 'ro', label='data')
  plt.step(en, flux, 'r-', where='mid')
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title[2], fontsize=fontsize)
  plt.legend(loc=0)
  plt.grid(True)
  # adjust ---!
  plt.subplots_adjust(hspace=0.5)
  # save fig ---!
  extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(file.replace('.fits', '_errors.png'), bbox_inches=extent.expanded(1.3, 1.3))
  extent = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(file.replace('.fits', '.png'), bbox_inches=extent.expanded(1.3, 1.3))
  extent = ax3.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(file.replace('.fits', '_log.png'), bbox_inches=extent.expanded(1.3, 1.3))
  # show fig --
  plt.show() if show else None
  plt.close()
  return

# plot cslightcrv output ---!
def showLightCurve(file, figsize=(15,15), axisLim ='auto', title='lightcurve', yscale=('lin','log'), xscale=('lin','log'),
                   show = True, tex=False):

  fig = plt.figure(figsize=figsize)
  plt.rc('text', usetex=False) if tex else None
  ax1 = plt.subplot(211, yscale=yscale[0], xscale=xscale[0])
  ax2 = plt.subplot(212, yscale=yscale[1], xscale=xscale[1])

  for i in range(len(file)):
    with fits.open(file) as hdul:
     data = hdul[1].data
     t_mjd = data.field(0)  # days
     et_mjd = data.field(1)  # days
     prefact = data.field(6)  # ph/cm^2/s/MeV
     e_prefact = data.field(7)  # ph/cm^2/s/MeV
     index = data.field(8)
     e_index = data.field(9)
     TS = data.field(10)
     diff_uplim = data.field(11) # ph/cm^2/s/MeV
     flux_uplim = data.field(12) # ph/cm^2/s
     Eflux_uplim = data.field(13) # erg/cm^2/s

    pnts = []
    e_pnts = []
    t_pnts = []
    et_pnts = []
    ul_pnts = []
    eul_pnts = []
    tul_pnts = []
    etul_pnts = []
    # list flux point or upper limit ---!
    for el in range(len(data)):
      if TS[el] > 9 and 2.0*e_prefact[el] < prefact[el] :
        pnts.append(prefact[el])
        e_pnts.append(e_prefact[el])
        t_pnts.append(t_mjd[el])
        et_pnts.append(et_mjd[el])
      else :
        ul_pnts.append(diff_uplim[el])
        eul_pnts.append(0.5*diff_uplim[el])
        tul_pnts.append(t_mjd[el])
        etul_pnts.append(et_mjd[el])

    # linear ---!
    ax1.errorbar(t_pnts, pnts, xerr=et_pnts, yerr=e_pnts, fmt='o', mec='k', label='data')
    ax1.errorbar(tul_pnts, ul_pnts, xerr=[etul_pnts, etul_pnts], yerr=eul_pnts, uplims=True, fmt='bo', mec='k')
    ax1.axis(axisLim) if axisLim != 'auto' else None
    ax1.grid()
    ax1.set_xlabel('t (MJD)')
    ax1.set_ylabel('dN/dE (ph/$cm^2$/s/MeV)')
    ax1.set_title('lightcurve') if title == 'none' else plt.title(title)
    # log ---!
    ax2.errorbar(t_pnts, pnts, xerr=et_pnts, yerr=e_pnts, fmt='o', mec='k', label='data')
    ax2.errorbar(tul_pnts, ul_pnts, xerr=[etul_pnts, etul_pnts], yerr=eul_pnts, uplims=True, fmt='bo', mec='k')
    ax2.axis(axisLim) if axisLim != 'auto' else None
    ax2.grid()
    ax2.set_xlabel('t (MJD)')
    ax2.set_ylabel('dN/dE (ph/$cm^2$/s/MeV)')
    ax2.set_title('lightcurve') if title == 'none' else plt.title(title)
  ax1.legend()
  ax2.legend()
  # adjust ---!
  plt.subplots_adjust(hspace=0.5)
  # save fig ---!
  extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(file.replace('.fits', '.png'), bbox_inches=extent.expanded(1.3, 1.3))
  extent = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(file.replace('.fits', '_log.png'), bbox_inches=extent.expanded(1.3, 1.3))
  # show fig ---!
  plt.show() if show else None
  plt.close()
  return

# show TS map ---!
def showTSmap(file, reg='none', col='black', suffix='none', title='TS map', cbar='Significance TSV',
              xlabel='RA (deg)', ylabel='DEC (deg)', fontsize=12, show=True, tex=False):
  with fits.open(file) as hdul:
    data = hdul[0].data
    hdr = hdul[0].header

  plt.rc('text', usetex=False) if tex else None
  ax = plt.subplot(111)
  # load region ---!
  if reg != 'none':
    r = pyregion.open(reg).as_imagecoord(hdr)
    for i in range(len(r)):
      r[i].attr[1]['color'] = col
      patch_list, text_list = r.get_mpl_patches_texts()
      for p in patch_list:
        ax.add_patch(p)
      for t in text_list:
        ax.add_artist(t)
  # plot ---!
  plt.imshow(data, cmap='bwr')
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title[0], fontsize=fontsize)
  plt.colorbar().set_label(cbar, fontsize=fontsize)
  # save fig ---!
  if suffix != 'none':
    plt.savefig(file.replace('.fits', '_%s.png' % suffix))
  else:
    plt.savefig(file.replace('.fits', '.png'))
  #show fig ---!
  plt.show() if show == True else None
  plt.close()
  return

# butterfly flux + spectrum ---!
def showButterflySpectrum(file, spectrum, axisLim='auto', suffix='none', title='butterfly diagram', fontsize=12,
                          xlabel='Energy (TeV)', ylabel='Flux (erg/$cm^2$/s)', show=True, tex=False):
  data = np.loadtxt(file, delimiter=' ')
  energy = data[:, 0]
  intensity = data[:, 1]
  lerr = data[:, 2]
  uerr = data[:, 3]
  with fits.open(spectrum) as hdul:
    spc = hdul[1].data
    en = spc.field(0)
    en_down = spc.field(1)
    en_up = spc.field(2)
    flux_pnt = spc.field(3)
    err_flux = spc.field(4)
  # adjusting errors ---!
  en_err = (en_up - en_down) / 2

  # from intensity (ph/cm^2/s/MeV) find flux in (erg/cm^2/s) ---!
  (flux, f_lerr, f_uerr) = ((intensity, lerr, uerr) * (energy ** 2)) / (6.4215 * 1e5)

  plt.rc('text', usetex=False) if tex else None
  ax = plt.subplot(111, yscale='log', xscale='log')
  # plot ---!
  plt.plot(energy / 1e6, flux, 'b-', alpha=1, label='best fit')
  plt.fill_between(energy / 1e6, f_lerr, f_uerr, facecolor='blue', alpha=0.3)
  plt.scatter(en, flux_pnt, marker='+', c='r', label='spectrum', alpha=1)
  ax.axis(axisLim) if axisLim != 'auto' else None
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title, fontsize=fontsize)
  plt.legend(loc=0)
  plt.grid(True)
  # save fig ---!
  if suffix != 'none':
    plt.savefig(file.replace('.txt', '_%s.png' % suffix))
  else:
    plt.savefig(file.replace('.txt', '.png'))
  # show fig ---!
  plt.show() if show == True else None
  plt.close()
  return

# irf degradation via aeff ---!
def degradedIRF_3d(x, y, z, xlabel='x', ylabel='y', zlabel='z', title=None, c=('b'), zscale='linear', tex=False,
                   fontsize=14, figsize=(12,6), rotation=0, zlim=(0,1), alpha=(1), label=None, savefig=None, show=True):

  fig = plt.figure(figsize=figsize)
  plt.rc('text', usetex=False) if tex else None
  sns.set_style("whitegrid", {'axes.grid': False})
  ax = fig.add_subplot(111, projection='3d', zscale=zscale)
  plt.xticks(fontsize=fontsize, rotation=rotation)
  plt.yticks(fontsize=fontsize, rotation=rotation)

  curve = []
  for i in range(len(z)):
    ax.plot_surface(x, y, z[i], alpha=alpha[i], color=c[i], label=label[i])
    curve.append(Rectangle((0, 0), 1, 1, fc=c[i], fill=True))
  ax.set_zlim(zlim)
  ax.set_xlabel(xlabel, fontsize=fontsize, labelpad=fontsize)
  ax.set_ylabel(ylabel, fontsize=fontsize, labelpad=fontsize)
  ax.set_zlabel(zlabel, fontsize=fontsize, labelpad=fontsize)
  ax.set_title(title, fontsize=fontsize) if title != None else None
  ax.legend(curve, label, loc=0, fontsize=fontsize) if label != None else None
  plt.tight_layout()
  fig.savefig(savefig) if savefig != None else None
  plt.show() if show else None
  plt.close()
  return

# plot ebl interpolation ---!
def interp_ebl(x, y, savefig, type='linear', xlabel='x', ylabel='y', title='title',
               label=('y', 'y2'), fontsize=12, show=True, tex=False, sns_style=False):

  fig = plt.figure()
  plt.rc('text', usetex=False) if tex else None
  sns.set() if sns_style else None
  ax = plt.subplot(111, xscale='log', yscale='log')
  plt.plot(x[0], y[0], '.', label=label[0], c='g')
  plt.plot(x[1], y[1], 'o', c='k', markeredgecolor='k', markerfacecolor='none', label=label[1])
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.xlabel(xlabel, fontsize=fontsize)
  plt.title(title, fontsize=fontsize)
  plt.legend(loc=0)
  plt.tight_layout()
  fig.savefig(savefig)
  plt.show() if show else None
  plt.close()
  return

# SENSITIVITY ---!
def showSensitivity(x, y, savefig, xlabel='energy (GeV)', ylabel='sensitivity', label=('nominal', 'nominal/2'), xscale='log', yscale='log', title='', fontsize=12, marker=('.'), ratio=True, show=True, tex=False, sns_style=False):

  fig = plt.figure()
  plt.rc('text', usetex=tex) if tex else None
  sns.set() if sns_style else None

  if ratio:
    ax1 = plt.subplot(211, xscale=xscale, yscale=yscale)
  else:
    ax1 = plt.subplot(111, xscale=xscale, yscale=yscale)
  for i in range(len(y)):
    plt.plot(x[i], y[i], marker=marker[i], label=label[i])
  plt.ylabel(ylabel, fontsize=fontsize)
  plt.title(title, fontsize=fontsize)
  plt.legend(loc=0)

  if ratio:
    ax2 = plt.subplot(212, sharex=ax1, yscale='linear')
    for i in range(len(y)-1):
      plt.plot(x[0], y[0]/y[i+1])
    plt.axhline(0.5, ls='-.', c='r')
    plt.ylabel('ratio nom/deg', fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylim([0.,1.])
  plt.grid()
  plt.tight_layout()
  fig.savefig(savefig)
  plt.show() if show else None
  plt.close()
  return

# plot lightcurve from pipeline ---!
def plotLightCurve(flux, t1, uplims, t2, xerr, yerr, filename, temp_t, temp_f, c1=('b'), c2=('r'), lf=('flux'), lup=('upper limit'), alpha=0.5, fontsize=20, figsize=(20, 16), rotation=0, ylim=None, xlim=None, interp=False, tex=False, sns_style=False):

  fig = plt.figure(figsize=figsize)
  plt.rc('text', usetex=False) if tex else None
  sns.set() if sns_style else None

  ax = plt.subplot(111, yscale='log', xscale='log')
  plt.xticks(fontsize=fontsize, rotation=rotation)
  plt.yticks(fontsize=fontsize, rotation=rotation)
  # template ---!
  plt.plot(temp_t[20:50], temp_f[20:50], '-g', lw=5, label='expected (30GeV-150TeV)', alpha=1)
  # detections ---!
  for i, f in enumerate(flux):
    plt.errorbar(t1[i], f, xerr=np.array([xerr[i] for i in range(len(t1[i]))]),
                 # yerr=np.array([yerr[i] for i in range(len(t1[i]))]),
                 marker='o', c=c1[i], alpha=alpha, label=lf[i], ls='none', ms=10)
    plt.fill_between(t1[i][0],
                     f[0] - np.array([yerr[i] for i in range(len(t1[i]))]),
                     f[0] + np.array([yerr[i] for i in range(len(t1[i]))]),
                     alpha=0.3, color=c1[i], interpolate=True)
    plt.axvline(t1[i][0][-1], ls='-.', lw=2, c=c1[i])
    plt.text(t1[i][0][-1] + 100, 7e-9 - 1e-9 * (i + 1), '%d s' % t1[i][0][-1], color=c1[i], fontsize=fontsize)
  # upper limits ---!
  for i, u in enumerate(uplims):
    if len(u[0]) > 1:
      plt.errorbar(t2[i], u, xerr=np.array([xerr[i] for i in range(len(t2[i]))]),
                   marker='v', c=c2[i], alpha=alpha, label=lup[i], ls='none', uplims=True)
      plt.axhline(np.mean(u), ls='-.', lw=2, c=c2[i])
      plateau = float(np.mean(uplims[i][0])) * 1e9
      plt.text(30, float(np.mean(uplims[i][0])) * 1.1, '%0.2fe-9 ph/cm$^2$/s' % plateau,
               color=c2[i], fontsize=fontsize + 5)

  plt.plot(temp_t[20:50], temp_f[20:50], '-g', lw=5, label='expected (30GeV-150TeV)', alpha=1)
  plt.ylabel('Flux (ph/cm$^2$/s)', fontsize=fontsize)
  plt.xlabel('time (s)', fontsize=fontsize)
  plt.legend(fontsize=fontsize)
  plt.xlim(xlim) if not xlim else None
  plt.ylim(ylim) if not ylim else None
  plt.show()
  fig.savefig(filename)
  plt.close()
  return