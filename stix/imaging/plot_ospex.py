#!/usr/bin/env python
"""
    Plot ospex pipeline results
    Author: Hualin Xiao (hualin.xiao@fhnw.ch)
    Date: Aug 16, 2022
"""
from astropy.io import fits
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter


def as_si(x,ndp, error=None):
    """
    format legend entries
    Parameters:
        x: value
        y: number of decimal points
        error: error of x
    Returns:
     a latex string

    """
    if np.abs(x) <10**ndp:
        return r'{m:.{ndp:d}f}'.format(m=x, ndp=ndp) if error is None else r'{m:.{ndp:d}f}$\pm${error:.{ndp:d}f}'.format(
        m=x,error=error, ndp=ndp)
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    if error is not None:
        se = error/(10**int(e))
        print(se)
        res= '{x:0.{ndp:d}f}'.format(x=se, ndp=ndp)
        if int(e)!=0:
            return r'({m:s} $\pm $ {se: 0.{ndp:d}f})$\times$ 10$^{{{e:d}}}$'.format(m=m, 
                                                                                      se=se, ndp=ndp, e=int(e))
        
        return r'{m:s} $\pm $ {se: 0.{ndp:d}f}'.format(m=m, se=se, ndp=ndp)
    
    return r'{m:s}$\times$ 10$^{{{e:d}}}$ '.format(m=m, e=int(e)) if int(e)!=0 else s
def create_param_legend(fitfun, params, sigmas):
    """
    params[0] - Emission measure in units of 10^49
  params[1] - KT, plasma temperature in keV
  params[2] - Abundance relative to coronal for Fe, Ni, Si, and Ca. S as well at half the deviation from 1 as Fe.
  params[3] - Total integrated electron flux, in units of 10^35 electrons sec^-1.
  params[4] - Power-law index of the electron flux distribution function below the break energy.
  params[5] - Break energy in the electron flux distribution function in keV
  params[6] - Power-law index of the electron flux distribution function above the break energy.
  params[7] - Low energy cutoff in the electron flux distribution function in keV.
  params[8] - High energy cutoff in the electron flux distribution function in keV.
  """

    entries = [
        f'EM: {as_si(params[0]*1e49 ,2, sigmas[0]*1e49)} cm$^{{-3}}$ ',
        f'KT: {as_si(params[1] ,2, sigmas[1])} keV ']
    if 'thick2' in fitfun:
        entries.extend([
        r'$a_0$: {txt:s}  electrons s$^{{-1}}$'.format(txt=as_si(params[3]*1e35 ,2, sigmas[3]*1e35)),
        f'$\gamma$: {as_si(params[4] ,2, sigmas[4])}',
        r'E$_{{c}}$: {txt:s} keV'.format(txt=as_si(params[7] ,2, sigmas[7])),
    ])
    return entries
    

def plot_ospex_results(fname):
    """
     plot ospex results and print the created plot to png or pdf
    """
    meta={}
    hdu=fits.open(fname)
    ebands=hdu['ENEBAND'].data
    ebins=[i[1] for i in ebands]
    ebins_keV=[i[2]-i[1] for i in ebands]
    ebins_mid=[(i[2]+i[1])*0.5 for i in ebands]
    ebins.append(ebands[-1][2])
    ibins_low=range(len(ebins))
    ibins_mid=[i+0.5 for i in ibins_low]

    det_area=hdu['RATE'].header['GEOAREA'] # detector area in units of cm^2
    to_model_units=lambda x:  (x/ebins_keV)/det_area
    #convert counts/sec to cnts/sec*keV*cm2
    xray_flux=to_model_units(hdu['RATE'].data['RATE'][0])
    #xray_flux
    xray_flux_err=to_model_units(hdu['RATE'].data['STAT_ERR'][0])
    bkg_flux=to_model_units(hdu['RATE'].data['BK_RATE'][0])
    bkg_flux_err=to_model_units(hdu['RATE'].data['BK_ERROR'][0])

    #create a canvas 
    fig = plt.figure()
    gs = gridspec.GridSpec(nrows=2, ncols=1, figure=fig, height_ratios=[3, 1])
    ax = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0],sharex=ax)
    plt.subplots_adjust(hspace=.0)

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    #plot signal and bkg
    flux_lines= ax.stairs(xray_flux,ebins, label='STIX spectrum', color=colors[0])
    ax.errorbar(ebins_mid, xray_flux, yerr=xray_flux_err, fmt='none', color=colors[0])
    ax.stairs(bkg_flux,ebins, label='BKG spectrum', color=colors[1])
    ax.errorbar(ebins_mid, bkg_flux, yerr=bkg_flux_err, fmt='none', color=colors[1])

    #fit components
    fitfun=hdu['RATE'].header['FITFUNC']
    components={'MODEL_VTH':'vth'}
    if 'thick2' in fitfun:
        components['MODEL_THICK']='thick2'
        components['MODEL_TOTAL']=fitfun

    icolor=2
    goodness_of_fit=hdu['Fit Components'].data['GOODNESS_OF_FIT']
    for key,name in components.items():
        try:
            data=hdu['Fit Components'].data[key][0]
            ax.stairs(data, ebins, label=name, color=colors[icolor])
            icolor+=1
        except (KeyError, IndexError):
            continue

    #customize plot1
    max_y=np.max(xray_flux)
    min_y=np.min(xray_flux)/10.
    ax.set_ylim(min_y, 2*max_y)
    ax.set_xlim(ebins[0], ebins[-1])
    title=f"{hdu[0].header['DATE-OBS']} – {hdu[0].header['DATE-END']}"
    ax.grid('on')
    ax.legend(loc='upper right', fontsize='small')#, bbox_to_anchor=(1., 1.05))
    ax.set_title(title)

    ax.set_ylabel(r"Flux [cnts s$^{-1}$ keV$^{-1}$ cm$^{-2}$]")
    ax.set_yscale('log')

    # display fit results
    params=hdu['RATE'].data['PARAMS'][0]
    sigmas=hdu['RATE'].data['SIGMAS'][0]

    meta={
            'fitfun':fitfun,
            'params':params,
            'sigmas':sigmas,
            'goodness_of_fit':goodness_of_fit

            }
    res_text=create_param_legend(fitfun, params, sigmas)

    legend_entries=[goodness_of_fit]+res_text
    ax.text(ebins[0]*1.1,min_y*1.2,'\n'.join(legend_entries), alpha=0.7, fontweight='bold',fontsize='small', 
            color='navy', fontstyle='normal', fontfamily='sans-serif')

    #residual plots

    residuals=to_model_units(hdu['RATE'].data['RESIDUAL'][0])
    ax2.stairs(residuals, ebins)
    ax2.set_xlabel('Energy (keV)')
    ax2.set_xscale('log')
    ax2.set_ylabel('Residuals')
    plt.tick_params(axis='x', which='minor')
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))

    return  (fig, gs, ax, ax2, meta)
