# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib as mpl
import scipy.ndimage as sp
from matplotlib.patches import Rectangle
from scipy import stats
from scipy.stats import rayleigh, norm, chi2
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Circle
from scipy.ndimage.filters import gaussian_filter
from matplotlib.ticker import FormatStrFormatter
from matplotlib.image import NonUniformImage
from scipy.ndimage.filters import gaussian_filter

extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
extra2 = Line2D([0], [0], ls='-.', color='k', lw='1')

def hist1d(x, mean, nbin=20, hist=True, fontsize=15, color='b', xscale='linear', figsize=(15,12), rotation=0, alpha=0.5, lw=3, ls=('-', '--', '-.', ':'), title='gaussian fit', ax_thresh=None, xlabel='x', ylabel='y', leglabel='data', filename='hist1d_gauss.png', usetex=False, sns_style=False, show=True):

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111, xscale=xscale)
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	# plt.plot([],[], color='none', label='wbin=%.2fdeg' %width)
	for index, el in enumerate(x):
		if el[0] is list():
			el=el[0]
		sns.distplot(el, bins=nbin, kde=False, hist=hist, fit=norm, norm_hist=True, fit_kws={"color": color[index]}, color=color[index], hist_kws={'alpha':alpha}, label=leglabel[index])
		plt.axvline(mean[index], c=color[index], ls=ls[index], lw=lw, label='mean ~ %.1E' %mean[index]) if mean != None else None
	plt.title(title, fontsize=fontsize)
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.legend(fontsize=fontsize)

	plt.grid() if not sns_style else None
	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax


# HIST 1D GAUSSIAN DISTRIBUTION ---!
def hist1d_gauss(x, mean, loc=0, threshold=1, nbin=20, width=None, hist=True, fontsize=15, figsize=(15,12), color='b', alpha=0.5, lw=3, ls=('-', '--', '-.', ':'), title='gaussian fit', ax_thresh=0.2, xlabel='x', ylabel='y', leglabel='data', rotation=0, filename='hist1d_gauss.png', usetex=False, sns_style=False, show=True):

	if nbin == None:
		if width == None:
			print('Error: set either nbin or width')
		nbin = int(threshold/width)

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	# plt.plot([],[], color='none', label='wbin=%.2fdeg' %width)
	for index, el in enumerate(x):
		if el[0] is list():
			el=el[0]
		sns.distplot(el, bins=nbin, kde=False, hist=hist, fit=norm, norm_hist=True, fit_kws={"color": color[index]}, color=color[index], hist_kws={'alpha':alpha, 'range':[loc-threshold, loc+threshold]}, label=leglabel[index])
		plt.axvline(mean[index], c=color[index], ls=ls[index], lw=lw, label='mean ~ %.3fdeg' %mean[index]) if mean != None else None
	plt.axvline(loc, c='k', ls='-', lw=lw, label='true ~ %.3fdeg' %loc) if loc != None else None
	plt.title(title, fontsize=fontsize)
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.legend(fontsize=fontsize)
	plt.xlim([loc-ax_thresh, loc+ax_thresh])

	plt.grid() if not sns_style else None
	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax


# HIST 1D RAYLEIGH DISTRIBUTION ---!
def hist1d_rayleigh(x, mean, rayleigh_prms={'loc':0, 'scale':[1]}, threshold=1, nbin=None, width=None, hist=True, fontsize=15, figsize=(15,12), rotation=0, color='b', alpha=0.5, lw=3, ls=('-', '--', '-.', ':'),title='rayleigh fit', ax_thresh=0.2, xlabel='x', ylabel='y', leglabel='data', filename='hist1d_rayleigh.png', usetex=False, sns_style=False, show=True):

	if width == None:
		width = threshold/nbin
	if nbin == None:
		nbin = int(threshold/width)
	if nbin == None and width == None:
		print('Error: set either nbin or width')

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	# plt.plot([],[], color='none', label='wbin=%.2fdeg' %width)
	for index, el in enumerate(x):
		if el[0] is list():
			el=el[0]
		sns.distplot(el, bins=nbin, kde=False, hist=hist, fit=rayleigh, norm_hist=True, fit_kws={"color": color[index]}, color=color[index], hist_kws={'alpha':alpha, 'range':[0.0, threshold]}, label=leglabel[index])
		plt.axvline(mean[index], c=color[index], ls=ls[index], lw=lw, label='mean ~ %.3fdeg' %mean[index]) if mean != None else None
		if rayleigh_prms['scale'] != None:
			plt.axvline(rayleigh_prms['scale'][index], c=color[index], ls='-', lw=w, label='mode ~ %.3fdeg' %rayleigh_prms['scale'][index])
	plt.title(title, fontsize=fontsize)
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.legend(fontsize=fontsize)
	plt.xlim([rayleigh_prms['loc'], rayleigh_prms['loc']+ax_thresh]) if rayleigh_prms['loc'] != None else None

	plt.grid() if not sns_style else None
	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax


# RAYLEIGH CDF WITH CONFIDENCE INTERVAL ---!
def rayleigh_cdf(x, loc=0, scale=1, if_CI=True, probs=(0.6827, 0.9545, 0.9973, 0.99994), xlabel='x', title='x ~ RA(gamma) CDF', colors=('k', 'r', 'orange', 'm'), fontsize=15, figsize=(15,12), rotation=0, filename='theo_rayleigh_cdf.png', usetex=False, sns_style=False, show=False):

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	ax.plot(np.sort(x), stats.rayleigh.cdf(np.sort(x), loc=loc, scale=scale), ls='-', label='cdf')
	ax.axvline(scale, c='maroon', label='gamma')
	ax.axvline(np.std(x), c='maroon', ls=':', label='1 std =%.2f' %(np.std(x)))

	if if_CI is True:
		x_critical = []
		for i in range(len(probs)):
			x_critical.append(stats.rayleigh.ppf(q=probs[i], loc=loc, scale=scale))
			ax.axvline(x_critical[i], c=colors[i], ls='-.', label='x=%.2f, %.2f' %(x_critical[i],probs[i]*100)+'%')

	plt.ylabel('1-alpha', rotation=90, fontsize=fontsize)
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	ax.set_xlim(left=0)
	ax.set_ylim(bottom=0)
	plt.legend(loc=0)

	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax


# RAYLEIGH PDF WITH CONFIDENCE INTERVAL ---!
def rayleigh_pdf(x, loc=0, scale=1, if_CI=True, probs=(0.6827, 0.9545, 0.9973, 0.99994), xlabel='x', title='x ~ RA(gamma) CDF', colors=('k', 'r', 'orange', 'm'), fontsize=15, figsize=(15,12), rotation=0, filename='theo_rayleigh_cdf.png', usetex=False, sns_style=False, show=False):

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	ax.plot(np.sort(x), stats.rayleigh.pdf(np.sort(x), loc=loc, scale=scale), ls='-', label='cdf')
	ax.axvline(scale, c='maroon', label='gamma')
	ax.axvline(np.std(x), c='maroon', ls=':', label='1 std =%.2f' %(np.std(x)))

	if if_CI is True:
		x_critical = []
		for i in range(len(probs)):
			x_critical.append(stats.rayleigh.ppf(q=probs[i], loc=loc, scale=scale))
			ax.axvline(x_critical[i], c=colors[i], ls='-.', label='x=%.2f, %.2f' %(x_critical[i],probs[i]*100)+'%')

	plt.ylabel('counts density', rotation=90, fontsize=fontsize)
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	ax.set_xlim(left=0)
	ax.set_ylim(bottom=0)
	plt.legend(loc=0)

	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax


# 2D HISTOGRAM WITH RAYLEIGH CONFIDENCE INTERVAL ---!
def hist2d_rayleigh_CI(x, y, nbin=None, width=None, rayleigh_prms={'loc':0, 'scale':1}, xcentre=0, ycentre=0, interp=None, threshold=1, probs=(0.6827, 0.9545, 0.9973, 0.99994), colors=('k', 'r', 'orange', 'm'), ls=('-','--','-.',':'), cmap='gist_heat', lw=4, ms=2e2, ax_thresh=0.2, xlabel='x', ylabel='y', title='confidence intervals from theoretical distribution', fontsize=15 , figsize=(10,8), rotation=0, filename='hist2d_CIrayleigh.png', usetex=False, sns_style=False, show=False):

	xmean = np.mean(x)
	ymean = np.mean(y)

	if width is None:
		width = threshold/nbin
	if nbin is None:
		nbin = int(threshold/width)
	if nbin is None and width is None:
		print('Error: set either nbin or width')

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	if interp == None:
		h = plt.hist2d(x, y, bins=nbin, cmap=cmap, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
	else:
		h, xedges, yedges = np.histogram2d(x, y, bins=nbin, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
		h = h.T
		plt.imshow(h, interpolation=interp, cmap=cmap, extent=[xcentre - threshold, xcentre + threshold, ycentre - threshold, ycentre + threshold])
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	plt.scatter(xcentre, ycentre, c='w', marker='*', s=ms)
	plt.plot([], [], c='none', label='Reyleigh')
	for i in range(len(probs)):
		plt.plot([], [], c=colors[i], ls=ls[i], label='%.2f' % (probs[i] * 100) + '%')
		r = stats.rayleigh.ppf(q=probs[i], loc=rayleigh_prms['loc'], scale=rayleigh_prms['scale'])
		#        q = rayleigh['scale'] * np.sqrt(-2 * np.log(probs[i]))
		#        r = stats.rayleigh.ppf(q=q, loc=rayleigh['loc'], scale=rayleigh['scale'])
		cir = Circle(xy=(float(xmean), float(ymean)), radius=r, color=colors[i], lw=lw, ls=ls[i])
		cir.set_facecolor('none')
		ax.add_artist(cir)

	if interp == None:
		cbar = plt.colorbar(h[3], ax=ax).set_label('counts', fontsize=fontsize)
	plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal') if ax_thresh != None else None
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.legend(ncol=3, fontsize=fontsize)

	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax


# COVARIANCE EIGENVALUES ---!
def eigsorted(cov):
	vals, vecs = np.linalg.eigh(cov)
	order = vals.argsort()[::-1]

	return vals[order], vecs[:, order]


# 2D HISTOGRAM WITH GAUSSIAN COVARIANCE CONFIDENCE INTERVAL ---!
def hist2d_gauss_CI(x, y, nbin=None, width=None, xcentre=0, ycentre=0, threshold=1, nstd=(1, 2, 3, 5), lw=4, ls=('-','--','-.',':'), cmap='gist_heat', colors=('k', 'r', 'orange', 'm'), ax_thresh=0.2, xlabel='x', ylabel='y', interp=None, ms=2e2, title='confidence intervals from theoretical distribution', fontsize=15, figsize=(10,8), rotation=0, filename='hist2d_CIgauss.png', usetex=False, sns_style=False, show=False):

	xmean = np.mean(x)
	ymean = np.mean(y)

	if width is None:
		width = threshold/nbin
	if nbin is None:
		nbin = int(threshold/width)
	if nbin is None and width is None:
		print('Error: set either nbin or width')

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	if interp == None:
		h = plt.hist2d(x, y, bins=nbin, cmap=cmap, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
	else:
		h, xedges, yedges = np.histogram2d(x, y, bins=nbin, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
		h = h.T
		plt.imshow(h, interpolation=interp, cmap=cmap, extent=[xcentre - threshold, xcentre + threshold, ycentre - threshold, ycentre + threshold])
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	plt.scatter(xcentre, ycentre, c='w', marker='*', s=ms)
	plt.plot([], [], c='none', label='gauss')
	for i in range(len(nstd)):
		plt.plot([], [], c=colors[i], ls=ls[i], label='%d sgm' % (nstd[i]))
		cov = np.cov(x, y)
		vals, vecs = eigsorted(cov)
		theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
		w, v = 2 * nstd[i] * np.sqrt(vals)
		ell = Ellipse(xy=(float(xmean), float(ymean)), width=w, height=v, angle=float(theta), color=colors[i], lw=lw, ls=ls[i])
		ell.set_facecolor('none')
		ax.add_artist(ell)

	if interp == None:
		cbar = plt.colorbar(h[3], ax=ax).set_label('counts', fontsize=fontsize)
	plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal') if ax_thresh != None else None
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.legend(ncol=3, fontsize=fontsize)

	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax

# 2D HISTOGRAM WITH GAUSSIAN COVARIANCE CONFIDENCE INTERVAL ---!
def contour_gauss_CI(x, y, nbin=None, width=None, xcentre=0, ycentre=0, threshold=1, nstd=(1, 2, 3, 5), colors=('k', 'r', 'orange', 'm'), ax_thresh=0.2, xlabel='x', ylabel='y', interp=None, title='confidence intervals from theoretical distribution', fontsize=15, figsize=(10, 8), rotation=0, filename='hist2d_CIgauss.png', usetex=False, sns_style=False, show=False):

	xmean = np.mean(x)
	ymean = np.mean(y)

	if width is None:
		width = threshold/nbin
	if nbin is None:
		nbin = int(threshold/width)
	if nbin is None and width is None:
		print('Error: set either nbin or width')

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	h = plt.hist2d(x, y, bins=nbin, cmap='jet', range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])

	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	plt.scatter(xcentre, ycentre, c='w', marker='*', s=1e2)
	plt.plot([], [], c='none', label='gauss')
	for i in range(len(nstd)):
		plt.plot([], [], c=colors[i], label='%d sigma' % (nstd[i]))
		cov = np.cov(x, y)
		vals, vecs = eigsorted(cov)
		theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
		w, v = 2 * nstd[i] * np.sqrt(vals)
		ell = Ellipse(xy=(xmean, ymean), width=w, height=v, angle=theta, color=colors[i], lw=2)
		ell.set_facecolor('none')
		ax.add_artist(ell)

	cbar = plt.colorbar(h[3], ax=ax).set_label('counts', fontsize=fontsize)
	plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal') if ax_thresh != None else None
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.legend(ncol=3, fontsize=fontsize)

	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax

# 2D HISTOGRAM MAP ---!
def hist2d_map(x, y, trials, nbin=None, width=None, xcentre=0, ycentre=0, threshold=1, ax_thresh=0.2, xlabel='x', ylabel='y', title='probability map', fontsize=15, figsize=(10,8), rotation=0, filename='hist2d_map.png', if_CI=None, rayleigh={'loc':0, 'scale':1}, nstd=(1, 2, 3, 5), colors=('k', 'r', 'orange', 'm'), probs=(0.6827, 0.9545, 0.9973, 0.99994), smooth=True, usetex=False, sns_style=False, show=False):

	if width is None:
		width = threshold/nbin
	if nbin is None:
		nbin = int(threshold/width)
	if nbin is None and width is None:
		print('Error: set either nbin or width')

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)
	h = ax.hist2d(x, y, bins=nbin, cmap='jet', vmin=0.0, vmax=trials,
	              range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
	if smooth:
		plt.clf()
		# hist2d with numpy (invert axis since imshow stumbles them) ---!
		# sigma=2
		# X = gaussian_filter(x, sigma)
		# Y = gaussian_filter(y, sigma)
		heatmap, xedges, yedges = np.histogram2d(x, y, bins=nbin, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
		ax.imshow(heatmap, extent=extent, cmap='jet', interpolation='gaussian', filterrad=1, filternorm=True, resample=False, origin='lower')

	plt.scatter(xcentre, ycentre, c='w', marker='*', s=1e2)

	if if_CI is None:
		pass

	elif if_CI.lower() == 'rayleigh':
		xmean = np.mean(x)
		ymean = np.mean(y)
		plt.plot([], [], c='none', label='Reyleigh')
		for i in range(len(probs)):
			plt.plot([], [], c=colors[i], label='%.2f' % (probs[i] * 100) + '%')
			r = stats.rayleigh.ppf(q=probs[i], loc=rayleigh['loc'], scale=rayleigh['scale'])
			cir = Circle(xy=(xmean, ymean), radius=r, color=colors[i], lw=2)
			cir.set_facecolor('none')
			ax.add_artist(cir)

	elif if_CI.lower() == 'gauss' or if_CI.lower == 'covariance' or if_CI.lower == 'cov':
		xmean = np.mean(x)
		ymean = np.mean(y)
		plt.plot([], [], c='none', label='gauss')
		for i in range(len(nstd)):
			plt.plot([], [], c=colors[i], label='%.2f' % (nstd[i] * 100) + '%')
			cov = np.cov(x, y)
			vals, vecs = eigsorted(cov)
			theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
			w, v = 2 * nstd[i] * np.sqrt(vals)
			ell = Ellipse(xy=(xmean, ymean), width=w, height=v, angle=theta, color=colors[i], lw=2)
			ell.set_facecolor('none')
			ax.add_artist(ell)
	else:
		print('Error: if_CI parameter value not understood')

	m = plt.cm.ScalarMappable(cmap='jet')
	m.set_clim(0., trials/100)
	cbar = plt.colorbar(m, boundaries=np.linspace(0, 100, 11)).set_label('cts \%', fontsize=fontsize)
	plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal')
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize) if title!=None else None
	plt.axis('equal')

	plt.tight_layout()
	fig.savefig(filename)
	# show fig ---!
	plt.show() if show == True else None
	plt.close()
	return fig, ax


# WILKS THEOREM DIST FOR EMPTY FIELDS ---!
def ts_wilks(x, trials, df=1, nbin=None, width=None, ylim=None, xlim=None, show=False, fontsize=15, figsize=(15,12), rotation=0, xlabel='TS', ylabel='normalised counts', title='TS distribution (empty fields)', filename='wilks_preTrials.png', usetex=False, sns_style=False):

	if width is None:
		width = (max(x)-min(x))/nbin
	if nbin is None:
		nbin = int((max(x)-min(x))/width)
	if nbin is None and width is None:
		print('Error: set either nbin or width')

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111, yscale='log')
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)

	h, edges = np.histogram(x, bins=int(nbin), density=False, range=(min(x), max(x)))
	yerr = np.sqrt(h)/trials
	h = h/trials
	cbin = (edges[1:] + edges[:-1]) / 2
	xerr = (edges[:-1] - edges[1:]) / 2

	x2 = np.arange(0, 30, 1)
	plt.errorbar(cbin, h, fmt='k+', yerr=yerr, xerr=xerr, markersize=5, label='')
	plt.plot(x2, stats.chi2.pdf(x2, df=df), c='orange', lw=1, ls='-.', label='chi2(dof=%d)' %df)
	plt.plot(x2, stats.chi2.pdf(x2, df=df)/2, c='b', lw=1, ls='--', label='chi2/2(dof=%d)' %df)

	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.legend(loc=0, fontsize=fontsize)
	plt.xlim(xlim) if xlim is not None else None
	plt.ylim(ylim) if ylim is not None else None

	plt.grid()
	plt.tight_layout()
	fig.savefig(filename)

	# show fig ---!
	plt.show() if show == True else None
	plt.close()

	return fig, ax

# WILKS THEOREM P-VALUES FOR EMPTY FIELDS ---!
def p_values(x, trials, df=1, nbin=None, width=None, ylim=None, xlim=None, show=False, fontsize=15, figsize=(15,12), rotation=0, xlabel='h', ylabel='p-values', title='p-value (empty fields)', filename='pvalue_preTrials.png', usetex=False, sns_style=False):

	if width is None:
		width = (max(x)-min(x))/nbin
	if nbin is None:
		nbin = int((max(x)-min(x))/width)
	if nbin is None and width is None:
		print('Error: set either nbin or width')

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111, yscale='log')
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)

	h = np.empty(len(np.arange(int(max(x)))))
	p = np.empty(len(np.arange(int(max(x)))))
	cbin, xerr = [], []
	for i in range(int(max(x))):
		cbin.append(i+1)
		xerr.append(0.5)
		for idx, val in enumerate(x):
			if val >= i+1:
				h[i] += 1
	p = h/trials
	yerr = np.sqrt(h)/trials

	x2 = np.arange(min(x), max(x)+5, 1)
	plt.errorbar(cbin[0], p[0], yerr=yerr[0], xerr=xerr[0], fmt='k+', markersize=5)
	plt.errorbar(cbin[1:], p[1:], yerr=yerr[1:], xerr=xerr[1:], fmt='k+', markersize=5, label='')
	plt.plot(x2, (1 - stats.chi2.cdf(x2, df=df)), lw=1, ls='--', c='green', label='chi2(dof=%d)' %df)
	plt.plot(x2, (1 - stats.chi2.cdf(x2, df=df))/2, lw=1, ls='-.', c='maroon', label='chi2/2(dof=%d)' %df)
	# plt.legend(('chi2/2(dof=%d)', 'chi2(dof=%d)', 'ts'), loc=0, fontsize=fontsize)
	plt.axhline(3e-7, c='gray', ls=':', alpha=1, lw=2)
	plt.text(23, 2e-7, '5sigma', fontsize=fontsize, alpha=1)

	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.legend(loc=0, fontsize=fontsize)
	plt.xlim(xlim) if xlim is not None else None
	plt.ylim(ylim) if ylim is not None else None

	plt.grid()
	plt.tight_layout()
	fig.savefig(filename)

	# show fig ---!
	plt.show() if show == True else None
	plt.close()

	return fig, ax

# WILKS THEOREM P-VALUES FOR EMPTY FIELDS ---!
def ts_wilks_cumulative(x, trials, df=1, nbin=None, width=None, ylim=None, xlim=None, show=False, fontsize=15, figsize=(15,12), rotation=0, xlabel='h', ylabel='cumulative probability', title='p-value (empty fields)', filename='cumulative_preTrials.png', usetex=False, sns_style=False):

	if width is None:
		width = (max(x)-min(x))/nbin
	if nbin is None:
		nbin = int((max(x)-min(x))/width)
	if nbin is None and width is None:
		print('Error: set either nbin or width')

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	sns.set() if sns_style else None

	ax = plt.subplot(111)
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)

	h = np.empty(len(np.arange(int(max(x)))))
	p = np.empty(len(np.arange(int(max(x)))))
	cbin, xerr = [], []
	for i in range(int(max(x))):
		cbin.append(i+1)
		xerr.append(0.5)
		for idx, val in enumerate(x):
			if val >= i+1:
				h[i] += 1
	p = 1 - h/trials
	yerr = np.sqrt(h)/trials

	x2 = np.arange(min(x), max(x)+5, 1)
	plt.errorbar(cbin[0], p[0], yerr=yerr[0], xerr=xerr[0], fmt='k+', markersize=5)
	plt.errorbar(cbin[1:], p[1:], yerr=yerr[1:], xerr=xerr[1:], fmt='k+', markersize=5, label='')
	plt.plot(x2, stats.chi2.cdf(x2, df=df), lw=1, ls='-.', c='maroon', label='P(dof=%d)' %df)
	# plt.legend(('chi2/2(dof=%d)', 'chi2(dof=%d)', 'ts'), loc=0, fontsize=fontsize)
	plt.axhline(1-3e-7, c='gray', ls=':', lw=2, alpha=1)
	plt.text(1, 0.95, '5sigma', fontsize=fontsize, alpha=1)

	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.legend(loc=0, fontsize=fontsize)
	plt.xlim(xlim) if xlim is not None else None
	plt.ylim(ylim) if ylim is not None else None

	plt.grid()
	plt.tight_layout()
	fig.savefig(filename)

	# show fig ---!
	plt.show() if show == True else None
	plt.close()

	return fig, ax

def chi2_reduced(x, trials, df=1, nbin=None, width=None, var=True):
	np.seterr(divide='ignore', invalid='ignore')
	if width is None:
		width = (max(x)-min(x))/nbin
	if nbin is None:
		nbin = int((max(x)-min(x))/width)
	if nbin is None and width is None:
		print('Error: set either nbin or width')

	h, edges = np.histogram(x, bins=int(nbin), density=False, range=(0., max(x)))
	yerr = np.sqrt(h)/trials
	h = h/trials
	cbin = (edges[1:] + edges[:-1])/2
	p = (1 - stats.chi2.pdf(cbin, df=df))/2
	#p = stats.chi2.pdf(cbin, df=df)/2
	#err = yerr/h

	#print('values', h, '\nerrors', yerr, '\nerror perc', err)

	with np.errstate(invalid='raise'):
		if var:
			chi2 = 2*np.sum((h[1:] - p[1:])**2/h[1:])
			#chi2 = np.sum((h[1:] - p[1:])**2/err[1:])
		else:
			chi2 = 2*np.sum((h[1:] - p[1:])**2/h[1:])
			#chi2 = np.sum((h[1:] - p[1:])**2/err[1:])
		h[1:] = np.array(h[1:])
		N = np.count_nonzero(h[1:])
		chi2r = chi2 / (N - 1)
	return chi2, chi2r


# MANUAL NORMALISED HISTOGRAM ---!
def normedHist(x, trials=None, step=None, nbin=None, ylim=None, xlim=None, show=False, normed=True, xscale='linear', yscale='log', fontsize=15, figsize=(15,12), rotation=0, xlabel='x', ylabel='normalised counts', leglabel='legend', title='normed histogram', usetex=False, sns_style=False, usesns=False, filename='normed_histogram.png'):

	x = np.sort(x)
	if step is None:
		step = (max(x)-min(x))/nbin
	if nbin is None:
		nbin = int((max(x)-min(x))/step)
	if nbin is None and step is None:
		print('Error: set either nbin or step')
	if trials is None:
		trials = len(x)

	fig = plt.figure(figsize=figsize)
	if usetex:
		plt.rc('text', usetex=usetex)
	if sns_style:
		sns.set()
	else:
		plt.grid()

	ax = plt.subplot(111, yscale=yscale, xscale=xscale)
	plt.xticks(fontsize=fontsize, rotation=rotation)
	plt.yticks(fontsize=fontsize, rotation=rotation)

	h = np.empty(len(np.arange(nbin)))
	cbin, xerr = [], []
	for i in range(nbin):
		cbin.append(step*i + step/2)
		xerr.append(step/2)
		for idx, val in enumerate(x):
			if val <= cbin[i]:
				h[i] += 1
		x=x[(x>=cbin[i])]
	h_norm = h/trials

	if normed:
		yerr = np.sqrt(h) / trials
		plt.errorbar(cbin, h_norm, yerr=yerr, xerr=xerr, fmt='k+', markersize=5, label=leglabel)
	else:
		yerr = h / trials
		plt.errorbar(cbin, h, yerr=yerr, xerr=xerr, fmt='k+', markersize=5, label=leglabel)
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.legend(loc=0, fontsize=fontsize)
	plt.xlim(xlim) if xlim is not None else None
	plt.ylim(ylim) if ylim is not None else None

	plt.tight_layout()
	fig.savefig(filename)

	# show fig ---!
	plt.show() if show else None
	plt.close()

	return fig, ax