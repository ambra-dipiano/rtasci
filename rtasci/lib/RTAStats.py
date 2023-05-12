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
from matplotlib.patches import Rectangle
from scipy import stats
from scipy.stats import rayleigh, norm, chi2
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Circle
from scipy.interpolate import interp1d

extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
extra2 = Line2D([0], [0], ls='-.', color='k', lw='1')

def hist1d(x, mean, true=None, nbin=20, hist=True, fit=True, fontsize=20, color='b', xscale='linear', figsize=(15,12), rotation=0, alpha=0.5, lw=3, ls=('--', '-.', ':'), title='gaussian fit', yscale='linear', ax_thresh=None, xlabel='x', ylabel='y', leglabel='data', filename='hist1d_gauss.png', usetex=False, sns_style=False, show=True):
    '''Generate multiple 1d histogram in a single plot. Optionally, only the histogram or the fit can be visualised.'''

    fig = plt.figure(figsize=figsize)
    if usetex:
        plt.rc('text', usetex=usetex)
    sns.set() if sns_style else None

    ax = plt.subplot(111, xscale=xscale, yscale=yscale)
    plt.xticks(fontsize=fontsize, rotation=rotation)
    plt.yticks(fontsize=fontsize, rotation=rotation)
    for index, el in enumerate(x):
        if el[0] is list():
            el=el[0]
        if fit:
            sns.distplot(el, bins=nbin, kde=False, hist=hist, fit=norm, norm_hist=True, fit_kws={"color": color[index], "ls": ls[index]}, color=color[index], hist_kws={'alpha':alpha})
            if mean != None:
                plt.axvline(mean[index], c=color[index], ls=ls[index], lw=lw, label=leglabel[index])
        else:
            sns.distplot(el, bins=nbin, kde=False, hist=hist, fit=None, norm_hist=True, color=color[index], hist_kws={'alpha':alpha}, label=leglabel[index])
        if true != None and type(true) == list:
            plt.axvline(true[index], c='k', ls=ls[index], lw=lw, label=f'{leglabel[index]} (expected)')
    if true != None and type(true) == float:
        plt.axvline(true, c='k', ls='-', lw=lw, label=f'expected')
    plt.title(title, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.legend(fontsize=fontsize)

    plt.grid() if not sns_style else None
    plt.tight_layout()
    fig.savefig(filename)
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax


# HIST 1D GAUSSIAN DISTRIBUTION
def hist1d_gauss(x, mean, true=None, loc=0, threshold=1, nbin=20, width=None, hist=True, fontsize=20, figsize=(15,12), color='b', alpha=0.5, lw=3, ls=('--', '-.', ':'), title='gaussian fit', ax_thresh=0.2, xlabel='x', ylabel='y', leglabel='data', rotation=0, filename='hist1d_gauss.png', usetex=False, sns_style=False, show=True):
    '''Generate multiple 1d histogram in a single plot. Optionally, only the histogram or the fit can be visualised.'''

    if nbin == None:
        if width == None:
            raise ValueError('Either nbin or width must be not None')
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
        sns.distplot(el, bins=nbin, kde=False, hist=hist, fit=norm, norm_hist=True, fit_kws={"color": color[index], "ls": ls[index], "lw:": lw}, color=color[index], hist_kws={'alpha':alpha, 'range':[loc-threshold, loc+threshold]}, label=leglabel[index])
        plt.axvline(mean[index], c=color[index], ls=ls[index], lw=lw) if mean != None else None
    if true != None:
        plt.axvline(true, c='k', ls='-', lw=lw, label='true')
    plt.axvline(loc, c='k', ls='-', lw=lw, label='true ~ %.3fdeg' %loc) if loc != None else None
    plt.title(title, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.xlim([loc-ax_thresh, loc+ax_thresh])

    plt.grid() if not sns_style else None
    plt.tight_layout()
    fig.savefig(filename)
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax


# HIST 1D RAYLEIGH DISTRIBUTION
def hist1d_rayleigh(x, mean, rayleigh_prms={'loc':0, 'scale':[1]}, threshold=1, nbin=None, width=None, hist=True, fontsize=15, figsize=(15,12), rotation=0, color='b', alpha=0.5, lw=3, ls=('-', '--', '-.', ':'),title='rayleigh fit', ax_thresh=0.2, xlabel='x', ylabel='y', leglabel='data', filename='hist1d_rayleigh.png', usetex=False, sns_style=False, show=True):
    '''Generate multiple 1d histogram in a single plot. Optionally, only the histogram or the fit can be visualised.'''

    if width == None:
        width = threshold/nbin
    if nbin == None:
        nbin = int(threshold/width)
    if nbin == None and width == None:
        raise ValueError('Either nbin or width must be not None')

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
            plt.axvline(rayleigh_prms['scale'][index], c=color[index], ls='-', lw=lw, label='mode ~ %.3fdeg' %rayleigh_prms['scale'][index])
    plt.title(title, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.xlim([rayleigh_prms['loc'], rayleigh_prms['loc']+ax_thresh]) if rayleigh_prms['loc'] != None else None

    plt.grid() if not sns_style else None
    plt.tight_layout()
    fig.savefig(filename)
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax


# RAYLEIGH CDF WITH CONFIDENCE INTERVAL
def rayleigh_cdf(x, loc=0, scale=1, if_CI=True, probs=(0.6827, 0.9545, 0.9973, 0.99994), xlabel='x', title='x ~ RA(gamma) CDF', colors=('k', 'r', 'orange', 'm'), fontsize=15, figsize=(15,12), rotation=0, filename='theo_rayleigh_cdf.png', usetex=False, sns_style=False, show=False):
    '''Plots the cumulative of Rayleigh distributed data, with given confidence interval.'''

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
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax


# RAYLEIGH PDF WITH CONFIDENCE INTERVAL
def rayleigh_pdf(x, loc=0, scale=1, if_CI=True, probs=(0.6827, 0.9545, 0.9973, 0.99994), xlabel='x', title='x ~ RA(gamma) CDF', colors=('k', 'r', 'orange', 'm'), fontsize=15, figsize=(15,12), rotation=0, filename='theo_rayleigh_cdf.png', usetex=False, sns_style=False, show=False):
    '''Plots the probability distribution of Rayleigh distributed data, with given confidence interval.'''

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
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax


# 2D HISTOGRAM WITH RAYLEIGH CONFIDENCE INTERVAL
def hist2d_rayleigh_CI(x, y, nbin=None, width=None, rayleigh_prms={'loc':0, 'scale':1}, xcentre=0, ycentre=0, interp=None, threshold=1, probs=(0.6827, 0.9545, 0.9973, 0.99994), colors=('k', 'r', 'orange', 'm'), ls=('-','--','-.',':'), cmap='gist_heat', lw=4, ms=2e2, ax_thresh=0.2, xlabel='x', ylabel='y', title='confidence intervals from theoretical distribution', fontsize=20 , figsize=(12,10), rotation=0, filename='hist2d_CIrayleigh.png', usetex=False, sns_style=False, show=False):
    '''Plots a 2d histrogram with Rayleigh confidence regions.'''

    xmean = np.mean(x)
    ymean = np.mean(y)

    if width is None:
        width = threshold/nbin
    if nbin is None:
        nbin = int(threshold/width)
    if nbin is None and width is None:
        raise ValueError('Either nbin or width must be not None')

    fig = plt.figure(figsize=figsize)
    if usetex:
        plt.rc('text', usetex=usetex)
    sns.set() if sns_style else None

    ax = plt.subplot(111)
    if interp == None:
        h = plt.hist2d(x, y, bins=nbin, cmap=cmap, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
        plt.colorbar(h[3], ax=ax).set_label('counts', fontsize=fontsize)
    else:
        h, edges = np.histogram2d(x, y, bins=nbin, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
        h = h.T
        im = plt.imshow(h, interpolation=interp, cmap=cmap, extent=[xcentre - threshold, xcentre + threshold, ycentre - threshold, ycentre + threshold])
        plt.colorbar(im, ax=ax).set_label('counts', fontsize=fontsize)
    plt.xticks(fontsize=fontsize, rotation=rotation)
    plt.yticks(fontsize=fontsize, rotation=rotation)
    plt.scatter(xcentre, ycentre, c='w', marker='*', s=ms)
    plt.plot([], [], c='none', label='Rayleigh')
    for i in range(len(probs)):
        plt.plot([], [], c=colors[i], ls=ls[i], label='%.2f' % (probs[i] * 100) + '%')
        r = stats.rayleigh.ppf(q=probs[i], loc=rayleigh_prms['loc'], scale=rayleigh_prms['scale'])
        #        q = rayleigh['scale'] * np.sqrt(-2 * np.log(probs[i]))
        #        r = stats.rayleigh.ppf(q=q, loc=rayleigh['loc'], scale=rayleigh['scale'])
        cir = Circle(xy=(float(xmean), float(ymean)), radius=r, color=colors[i], lw=lw, ls=ls[i])
        cir.set_facecolor('none')
        ax.add_artist(cir)

    plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal') if ax_thresh != None else None
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(ncol=3, fontsize=fontsize)

    plt.grid()
    plt.tight_layout()
    fig.savefig(filename)
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax


# COVARIANCE EIGENVALUES
def eigsorted(cov):
    '''Returns covariance eigenvalues.'''
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:, order]


# 2D HISTOGRAM WITH GAUSSIAN COVARIANCE CONFIDENCE INTERVAL
def hist2d_gauss_CI(x, y, nbin=None, width=None, xcentre=0, ycentre=0, threshold=1, nstd=(1, 2, 3, 5), lw=4, ls=('-','--','-.',':'), cmap='gist_heat', colors=('k', 'r', 'orange', 'm'), ax_thresh=0.2, xlabel='x', ylabel='y', interp=None, ms=2e2, title='confidence intervals from theoretical distribution', fontsize=20, figsize=(12,10), rotation=0, filename='hist2d_CIgauss.png', usetex=False, sns_style=False, show=False):
    '''Plots a 2d histrogram with Gaussian confidence regions.'''

    xmean = np.mean(x)
    ymean = np.mean(y)

    if width is None:
        width = threshold/nbin
    if nbin is None:
        nbin = int(threshold/width)
    if nbin is None and width is None:
        raise ValueError('Either nbin or width must be not None')

    fig = plt.figure(figsize=figsize)
    if usetex:
        plt.rc('text', usetex=usetex)
    sns.set() if sns_style else None

    ax = plt.subplot(111)
    if interp == None:
        h = plt.hist2d(x, y, bins=nbin, cmap=cmap, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
        cbar = plt.colorbar(h[3], ax=ax).set_label('counts', fontsize=fontsize)
    else:
        h, xedges, yedges = np.histogram2d(x, y, bins=nbin, range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
        h = h.T
        img = plt.imshow(h, interpolation=interp, cmap=cmap, extent=[xcentre - threshold, xcentre + threshold, ycentre - threshold, ycentre + threshold])
        cbar = plt.colorbar(img, ax=ax).set_label('counts', fontsize=fontsize)
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

    plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal') if ax_thresh != None else None
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(ncol=3, fontsize=fontsize)

    plt.grid()
    plt.tight_layout()
    fig.savefig(filename)
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax

# 2D HISTOGRAM WITH GAUSSIAN COVARIANCE CONFIDENCE INTERVAL
def contour_gauss_CI(x, y, nbin=None, width=None, xcentre=0, ycentre=0, threshold=1, nstd=(1, 2, 3, 5), colors=('k', 'r', 'orange', 'm'), ax_thresh=0.2, xlabel='x', ylabel='y', interp=None, title='confidence intervals from theoretical distribution', fontsize=20, figsize=(10, 8), rotation=0, filename='hist2d_CIgauss.png', usetex=False, sns_style=False, show=False):
    '''Plots Gaussian contour map.'''

    xmean = np.mean(x)
    ymean = np.mean(y)

    if width is None:
        width = threshold/nbin
    if nbin is None:
        nbin = int(threshold/width)
    if nbin is None and width is None:
        raise ValueError('Either nbin or width must be not None')

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
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax

# 2D HISTOGRAM MAP
def hist2d_map(x, y, nbin=None, width=None, xcentre=0, ycentre=0, threshold=1, ax_thresh=0.2, xlabel='x', ylabel='y', title='probability map', fontsize=20, figsize=(12,10), rotation=0, filename='hist2d_map.png', if_CI=None, rayleigh={'loc':0, 'scale':1}, nstd=(1, 2, 3, 5), colors=('k', 'r', 'orange', 'm'), probs=(0.6827, 0.9545, 0.9973, 0.99994), smooth=True, usetex=False, sns_style=False, show=False):
    '''Produces a smoothed count map.'''

    trials = len(x)
    if width is None:
        width = threshold/nbin
    if nbin is None:
        nbin = int(threshold/width)
    if nbin is None and width is None:
        raise ValueError('Either nbin or width must be not None')

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
        raise ValueError('Ivalid if_CI parameter value')

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
    # show fig
    plt.show() if show == True else None
    plt.close()
    return fig, ax


# WILKS THEOREM DIST FOR EMPTY FIELDS
def ts_wilks(x, df=1, nbin=None, width=None, trials=None, xrange=None, ylim=None, xlim=None, show=False, fontsize=15, figsize=(15,12), rotation=0, xlabel='TS', ylabel='normalised counts', title='TS distribution (empty fields)', filename='wilks_preTrials.png', usetex=False, sns_style=False, overlay=['chi2'], write_data=False, dpi=400, fmt='+', ecolor='red', markersize=1, elinewidth=1, alpha=0.8):
    '''Plots a TS distribution comparison with chi2 and chi2/2.'''
    assert width == None or nbin == None, 'Define either nbin or width but bot both.'

    x = np.array(x)
    filename = str(filename)
    fig = plt.figure(figsize=figsize)
    if usetex:
        plt.rc('text', usetex=usetex)
    sns.set() if sns_style else None

    ax = plt.subplot(111, yscale='log')
    plt.xticks(fontsize=fontsize, rotation=rotation)
    plt.yticks(fontsize=fontsize, rotation=rotation)

    for n, el in enumerate(x):
        # checks
        if xrange == None:
            xrange = (min(el), max(el))
        el = el[(el >= xrange[0]) & (el <= xrange[1])]
        if trials == None:
            trials = len(el)
        if width is None:
            width = (xrange[1]-xrange[0])/nbin
        if nbin is None:
            nbin = int((xrange[1]-xrange[0])/width)
        if nbin is None and width is None:
            raise ValueError('Either nbin or width must be not None')
        if type(overlay) != list:
            overlay = [overlay]

        # compute the histogram
        h, edges = np.histogram(el, bins=int(nbin), density=False, range=(xrange[0], xrange[1]))
        yerr = np.sqrt(h)/trials
        h = h/trials
        cbin = (edges[1:] + edges[:-1]) / 2
        xerr = (edges[:-1] - edges[1:]) / 2

        # plot the histogram
        plt.errorbar(cbin, h, yerr=yerr, xerr=xerr, fmt=fmt, ecolor=ecolor, markersize=markersize, elinewidth=elinewidth, alpha=alpha, label=f'hist ({n})')
        if 'mplt' in overlay:   
            plt.hist(el, bins=nbin, density=False, histtype='step', align='mid', range=(xrange[0], xrange[1]), label=f'mplt_{n}')

        # save the istogram
        if write_data:
            save_hist_on_file(x=cbin, y=h, xerr=xerr, yerr=yerr, filename=filename.replace('.png', f'_{n}.txt'))

    # overlay theoretical dist
    if 'chi2' in overlay:
        x2 = np.linspace(xrange[0], xrange[1], nbin)
        plt.plot(x2, chi2.pdf(x2, df=df), c='orange', lw=1, ls='-.', label=f'chi2(dof={df})')
        plt.plot(x2, chi2.pdf(x2, df=df)/2, c='b', lw=1, ls='--', label=f'chi2/2(dof={df})')

    # decorations
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(loc=0, fontsize=fontsize)
    plt.xlim(xlim) if xlim is not None else None
    plt.ylim(ylim) if ylim is not None else None

    plt.grid()
    plt.tight_layout()
    fig.savefig(filename, dpi=dpi)

    # show fig
    plt.show() if show == True else None
    plt.close()

    return fig, ax

# WILKS THEOREM P-VALUES FOR EMPTY FIELDS
def p_values(x, df=1, nbin=None, width=None, trials=None, xrange=None, ylim=None, xlim=None, show=False, fontsize=15, figsize=(15,12), rotation=0, xlabel='h', ylabel='p-values', slabel=None, title='p-value (empty fields)', filename='pvalue_preTrials.png', usetex=False, sns_style=False, overlay=['chi2'], sigma5=True, write_data=False, dpi=400, fmt='+', ecolor='red', markersize=1, elinewidth=1, alpha=0.8, legendprop={"size": 6}, legendloc=0):
    '''Plots a p-values distribution comparison with chi2 and chi2/2.'''
    assert width == None or nbin == None, 'Define either nbin or width but bot both.'
    markers_colors = ['#2ca25f', '#8856a7', '#e34a33', '#4a0a41']
    error_lines_colors = ['#99d8c9', '#9ebcda', '#fdbb84', '#f571e3']
    if len(x) > 4:
        print("Warning: colors are optimized for a maximum of three series.")
    x = np.array(x)
    filename = str(filename)
    fig = plt.figure(figsize=figsize)
    if usetex:
        plt.rc('text', usetex=usetex)
    sns.set() if sns_style else None

    ax = plt.subplot(111, yscale='log')
    plt.xticks(fontsize=fontsize, rotation=rotation)
    plt.yticks(fontsize=fontsize, rotation=rotation)

    if slabel is None:
        slabel = list(range(len(x)))

    for n, el in enumerate(x):
        # checks
        if xrange == None:
            xrange = (min(el), max(el))
        el = el[(el >= xrange[0]) & (el <= xrange[1])]
        if trials == None:
            trials = len(el)
        if width is None:
            width = (xrange[1]-xrange[0])/nbin
        if nbin is None:
            nbin = int((xrange[1]-xrange[0])/width)
        if nbin is None and width is None:
            raise ValueError('Either nbin or width must be not None')
        if type(overlay) != list:
            overlay = [overlay]

        # compute pvalues
        xerr = np.full((len(range(nbin))), width/2)
        h, edges = np.histogram(el, bins=nbin, range=xrange)
        cumul = np.cumsum(h[::-1])[::-1] 

        p = cumul/trials
        yerr = np.sqrt(cumul)/trials
        cbin = (edges[1:] + edges[:-1]) / 2

        # plot the pvalues
        plt.errorbar(cbin, p, yerr=yerr, xerr=xerr, fmt=fmt, markerfacecolor=markers_colors[n], markeredgecolor=markers_colors[n], ecolor=error_lines_colors[n], markersize=markersize, elinewidth=elinewidth,  alpha=alpha, label=f'p-values ({slabel[n]})')
        
        if 'mplt' in overlay:
            plt.hist(el, bins=nbin, density=True, histtype='step', align='mid', range=(xrange[0], xrange[1]), cumulative=-1, label=f'mplt_{n}')

        # save the pvalues
        if write_data:
            save_hist_on_file(x=cbin, y=p, xerr=xerr, yerr=yerr, filename=filename.replace('.png', f'_{n}.numpy.txt'))

    # overlay theoretical dist
    if 'chi2' in overlay:
        x2 = np.linspace(xrange[0], xrange[1], nbin)
        plt.plot(x2, (1 - chi2.cdf(x2, df=df)), lw=1, ls='--', c='green', label='chi2(dof=%d)' %df)
        plt.plot(x2, (1 - chi2.cdf(x2, df=df))/2, lw=1, ls='-.', c='maroon', label='chi2/2(dof=%d)' %df)
        plt.legend(('chi2/2(dof=%d)', 'chi2(dof=%d)', 'ts'), loc=0, fontsize=fontsize)

    # 5 sigma threshold
    if sigma5:
        plt.axhline(3e-7, c='gray', ls=':', alpha=1, lw=2)
        plt.text(xrange[0]*1.2, 3e-7, '5sigma', fontsize=fontsize, alpha=1)
        f = interp1d(cbin, p, fill_value="extrapolate", kind='linear')
        ts = f(3e-7)
        print(f'Significance (5sgm) == {ts}')

    # decorations
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(loc=legendloc, fontsize=fontsize, prop=legendprop, markerscale=10)
    plt.xlim(xlim) if xlim is not None else None
    plt.ylim(ylim) if ylim is not None else None

    plt.grid()
    plt.tight_layout()
    fig.savefig(filename, dpi=dpi)

    # show fig 
    plt.show() if show == True else None
    plt.close()

    return fig, ax

# WILKS THEOREM P-VALUES FOR EMPTY FIELDS
def ts_wilks_cumulative(x, df=1, nbin=None, width=None, trials=None, xrange=None, ylim=None, xlim=None, show=False, fontsize=15, figsize=(15,12), rotation=0, xlabel='h', ylabel='cumulative probability', title='p-value (empty fields)', filename='cumulative_preTrials.png', usetex=False, sns_style=False, overlay=['chi2'], write_data=False, sigma5=True, dpi=400, fmt='+', ecolor='red', markersize=1, elinewidth=1, alpha=0.8):
    '''Plots a TS cumulative distribution.'''
    assert width == None or nbin == None, 'Define either nbin or width but bot both.'

    x = np.array(x)
    filename = str(filename)
    fig = plt.figure(figsize=figsize)
    if usetex:
        plt.rc('text', usetex=usetex)
    sns.set() if sns_style else None

    ax = plt.subplot(111, yscale='log')
    plt.xticks(fontsize=fontsize, rotation=rotation)
    plt.yticks(fontsize=fontsize, rotation=rotation)

    for n, el in enumerate(x):
        # checks
        if xrange == None:
            xrange = (min(el), max(el))
        el = el[(el >= xrange[0]) & (el <= xrange[1])]
        if trials == None:
            trials = len(el)
        if width is None:
            width = (xrange[1]-xrange[0])/nbin
        if nbin is None:
            nbin = int((xrange[1]-xrange[0])/width)
        if nbin is None and width is None:
            raise ValueError('Either nbin or width must be not None')
        if type(overlay) != list:
            overlay = [overlay]

        # compute cumulative
        xerr = np.full((len(range(nbin))), width/2)
        h, edges = np.histogram(el, bins=nbin, range=xrange)
        cumul = np.cumsum(h[::-1])[::-1] 

        p = 1 - cumul/trials
        yerr = np.sqrt(cumul)/trials
        cbin = (edges[1:] + edges[:-1]) / 2

        # plot the cumulative
        plt.errorbar(cbin, p, yerr=yerr, xerr=xerr, fmt=fmt, markersize=markersize, ecolor=ecolor, elinewidth=elinewidth, alpha=alpha, label='cumulative')
        if 'mplt' in overlay:
            plt.hist(el, bins=nbin, density=True, histtype='step', align='mid', range=(xrange[0], xrange[1]), cumulative=-1, label=f'mplt_{n}')

        # save the pvalues
        if write_data:
            save_hist_on_file(x=cbin, y=p, xerr=xerr, yerr=yerr, filename=filename.replace('.png', f'_{n}.txt'))

   # overlay theoretical dist
    if 'chi2' in overlay:
        x2 = np.linspace(xrange[0], xrange[1], nbin)
        plt.plot(x2, (1 - chi2.cdf(x2, df=df)), lw=1, ls='--', c='green', label='chi2(dof=%d)' %df)
        plt.plot(x2, (1 - chi2.cdf(x2, df=df))/2, lw=1, ls='-.', c='maroon', label='chi2/2(dof=%d)' %df)
        plt.legend(('chi2/2(dof=%d)', 'chi2(dof=%d)', 'ts'), loc=0, fontsize=fontsize)

    # 5 sigma threshold
    if sigma5:
        plt.axhline(1-3e-7, c='gray', ls=':', alpha=1, lw=2)
        plt.text(xrange[0]*1.2, 1-3e-7, '5sigma', fontsize=fontsize, alpha=1)
        f = interp1d(cbin, p, fill_value="extrapolate", kind='linear')
        ts = f(1-3e-7)
        print(f'Significance (5sgm) == {ts}')

    # decorations
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(loc=0, fontsize=fontsize)
    plt.xlim(xlim) if xlim is not None else None
    plt.ylim(ylim) if ylim is not None else None

    plt.grid()
    plt.tight_layout()
    fig.savefig(filename, dpi=dpi)

    # show fig 
    plt.show() if show == True else None
    plt.close()

    return fig, ax

def chi2_reduced(x, df=1, nbin=None, width=None, var=True, xrange=None):
    '''Returns the chi2 and chi2 reduced.'''
    trials = len(x)
    np.seterr(divide='ignore', invalid='ignore')

    if xrange == None:
        xrange = (min(x), max(x))
    if trials == None:
        trials = len(x)
    if width is None:
        width = (xrange[1]-xrange[0])/nbin
    if nbin is None:
        nbin = int((xrange[1]-xrange[0])/width)
    if nbin is None and width is None:
        raise ValueError('Either nbin or width must be not None')

    h, edges = np.histogram(x, bins=int(nbin), density=False, range=(0., xrange[1]))
    yerr = np.sqrt(h)/trials
    h = h/trials
    cbin = (edges[1:] + edges[:-1])/2
    p = (1 - stats.chi2.pdf(cbin, df=df))/2
    #p = chi2.pdf(cbin, df=df)/2
    #err = yerr/h

    #print('values', h, '\nerrors', yerr, '\nerror perc', err)

    with np.errstate(invalid='raise'):
        if var:
            chi2n = 2*np.sum((h[1:] - p[1:])**2/h[1:])
            #chi2n = np.sum((h[1:] - p[1:])**2/err[1:])
        else:
            chi2n = 2*np.sum((h[1:] - p[1:])**2/h[1:])
            #chi2n = np.sum((h[1:] - p[1:])**2/err[1:])
        h[1:] = np.array(h[1:])
        N = np.count_nonzero(h[1:])
        chi2r = chi2n / (N - 1)
    return chi2n, chi2r


# MANUAL NORMALISED HISTOGRAM
def make_hist(x, step=None, nbin=None, width=None, normed=True, write_data=False, xrange=None, trials=None, filename='histogram.txt'):
    
    x = np.sort(x)
    if xrange == None:
        xrange = (min(x), max(x))
    if trials == None:
        trials = len(x)
    if width is None:
        width = (xrange[1]-xrange[0])/nbin
    if nbin is None:
        nbin = int((xrange[1]-xrange[0])/width)
    if nbin is None and width is None:
        raise ValueError('Either nbin or width must be not None')

    h = np.zeros_like(np.empty(len(np.arange(nbin))))
    cbin, xerr = [], []
    for i in range(nbin):
        cbin.append(step*i + step/2)
        xerr.append(step/2)
        for idx, val in enumerate(x):
            if val <= cbin[i]:
                h[i] += 1
        x=x[(x>=cbin[i])]

    if normed:
        yerr = np.sqrt(h) / trials
        h = h/trials
    else:
        yerr = np.sqrt(h)

    if write_data:
        save_hist_on_file(x=cbin, y=h, xerr=xerr, yerr=yerr, filename=filename.replace('.png', '.txt'))
    return cbin, h, xerr, yerr

def normed_hist_plot(x, step=None, nbin=None, width=None, ylim=None, xlim=None, trials=None, xrange=None, show=False, normed=True, xscale='linear', yscale='log', fontsize=15, figsize=(15,12), rotation=0, xlabel='x', ylabel='normalised counts', leglabel='legend', title='normed histogram', usetex=False, sns_style=False, usesns=False, filename='normed_histogram.png', write_data=False):
    '''Generates a manually normalised histogram.'''

    x = np.sort(x)
    if xrange == None:
        xrange = (min(x), max(x))
    if trials == None:
        trials = len(x)
    if width is None:
        width = (xrange[1]-xrange[0])/nbin
    if nbin is None:
        nbin = int((xrange[1]-xrange[0])/width)
    if nbin is None and width is None:
        raise ValueError('Either nbin or width must be not None')

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

    h = np.zeros_like(np.empty(len(np.arange(nbin))))
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

    if write_data:
        save_hist_on_file(x=cbin, y=h, xerr=xerr, yerr=yerr, filename=filename.replace('.png', '.txt'))

    # show fig
    plt.show() if show else None
    plt.close()

    return fig, ax


def save_hist_on_file(x, y, xerr, yerr, filename='histogram.txt', hdr='x xerr y yerr'):
    log = open(filename, 'w+')
    if hdr != None:
        if '\n' not in hdr:
            hdr += '\n'
        log.write(hdr)
    for i in range(len(x)):
        log.write(f'{x[i]} {xerr[i]} {y[i]} {yerr[i]}\n')
    log.close()
    return

def save_data_on_file(x, filename='histogram.txt', hdr='x xerr y yerr'):
    log = open(filename, 'w+')
    if hdr != None:
        if '\n' not in hdr:
            hdr += '\n'
        log.write(hdr)
    for el in x:
        log.write(f'{el}\n')
    log.close()
    return

def get_sigma_from_pvalue(pval, decimals=3):
    return np.abs(np.round(norm.ppf(pval), decimals))

def get_prob_from_sigma(sigma, decimals=8):
    return np.round(1-(norm.sf(sigma)*2), decimals)

def get_prob_from_pvalue(pval, decimals=8):
    return np.round(1-pval*2, decimals)

def get_pvalue_from_sigma(sigma, decimals=8):
    p = get_prob_from_sigma(sigma, decimals=decimals)
    return np.round((1-p)/2, decimals)