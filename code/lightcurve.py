import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy.coordinates import SkyCoord
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats


PATH = '/Users/erinkimbro/Projects/merian_variable'


#define arguments
parser = argparse.ArgumentParser()
parser.add_argument("ra", type=float, help="RA of object in decimal degrees")
parser.add_argument("dec", type=float, help="dec of object in decmial degrees")
args = parser.parse_args()

pos = SkyCoord(ra=args.ra, dec=args.dec, unit='deg')
name = pos.to_string('hmsdms')

hdr = fits.open(PATH + '/tables/measurements/'+name+'_results.fits')
t = hdr[1].data

hdr = fits.open(PATH + '/tables/models/'+name+'_model.fits')
model_t = hdr[1].data

#need to read in fit table -> search for object
fit = ascii.read(PATH+'/tables/fit.csv')
mask = fit['ra'] == args.ra
fit = fit[mask]

#need to read in tpl table -> search for object
tpl = ascii.read(PATH + '/tables/tpl_mags.csv')
mask = tpl['ra'] == args.ra
tpl = tpl[mask]

mag_mean, mag_median, mag_std = sigma_clipped_stats(t['mag'], sigma=3.0)

mask = (t['mag']<mag_mean+3*mag_std)&(t['mag']>mag_mean-3*mag_std)

mag = t['mag'][mask]
date = t['date'][mask]
mag_err = t['mag_err'][mask]


fig, ax = plt.subplots()
ax.errorbar(date, mag, yerr=mag_err, fmt='o', color='black', ecolor='grey', markersize=4, label='Data')
ax.hlines(tpl['mag'], np.min(date)-100, np.max(date)+100, color = 'black', label='Template Mag')
ax.fill_between(model_t['date'], model_t['mag']-model_t['mag_err'], model_t['mag']+model_t['mag_err'], 
                 alpha=0.3, color='teal')
ax.plot(model_t['date'], model_t['mag'], color='teal', label='Model')
ax.text(.05, .95, r'$\sigma_{vary}$='+str(np.round(fit['signif_vary'], 3)), transform=ax.transAxes)
ax.text(.05, .90, r'$\sigma_{QSO}$='+str(np.round(fit['signif_qso'], 3)), transform=ax.transAxes)
ax.text(.05, .85, r'$\sigma_{notQSO}$='+str(np.round(fit['signif_not_qso'], 3)), transform=ax.transAxes)

ax.set_xlabel('MJD')
ax.set_ylabel(r'$m_g$')
ax.set_xlim(np.min(date)-50, np.max(date)+50)
ax.set_ylim(np.max(mag)+5*np.std(mag), np.min(mag)-5*np.std(mag))
ax.legend()
#ax.gca().invert_yaxis()
plt.savefig(PATH+'/plot/light_curves/'+name+'_lightcurve.png')