import glob as glob
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats, CircularAperture, CircularAnnulus
import numpy as np
from astropy.table import Table
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from photutils.aperture import SkyCircularAperture, SkyCircularAnnulus
import csv
from astropy.wcs import WCS

#get files
PATH = '/Users/erinkimbro/Projects/merian_variable'
files = glob.glob(PATH + '/WorkingDir/1_1/s_rc_ztf*.fits')

tpl_file = glob.glob(PATH + '/WorkingDir/1_1/tpl.fits')

#RA=232.9219
#DEC=20.32498

#define arguments
parser = argparse.ArgumentParser()
parser.add_argument("ra", type=float, help="RA of object in decimal degrees")
parser.add_argument("dec", type=float, help="dec of object in decmial degrees")
args = parser.parse_args()

#photometry
mean = np.zeros(len(files))
median = np.zeros(len(files))
std = np.zeros(len(files))
peaks = np.zeros(len(files))
peaks_err = np.zeros(len(files))

seeing = np.zeros(len(files))
date = np.zeros(len(files))
flux = np.zeros(len(files))
flux_err = np.zeros(len(files))
mag = np.zeros(len(files))
mag_err = np.zeros(len(files))

#get template image
pos = SkyCoord(ra=args.ra, dec=args.dec, unit='deg')
name = pos.to_string('hmsdms')
hdr = fits.open(tpl_file[0])
wcs = WCS(hdr[0].header)
data = hdr[0].data
tpl_bkg_aper = SkyCircularAnnulus(pos, 5*u.arcsec, 10*u.arcsec)
tpl_bkg_aperstats = ApertureStats(data, tpl_bkg_aper, wcs=wcs)
aper = SkyCircularAperture(pos, 2*u.arcsec)
aperstats = ApertureStats(data, aper, local_bkg=tpl_bkg_aperstats.median, wcs=wcs)
tpl_peak = aperstats.sum
tpl_flux = hdr[0].header['GAIN']*tpl_peak/hdr[0].header['EXPOSURE'] #e-/s
tpl_flux_err = hdr[0].header['GAIN']*aperstats.std/hdr[0].header['EXPOSURE']
tpl_mag = -2.5*np.log10(tpl_flux) + hdr[0].header['MAGZP'] + hdr[0].header['CLRCOEFF']*hdr[0].header['CLRMED'] + hdr[0].header['APCOR1']
tpl_mag_err =  -2.5/(tpl_flux*np.log(10))*tpl_flux_err

#get diff flux
for i in range(len(files)):
    hdr = fits.open(files[i])
    data = hdr[0].data
    wcs = WCS(hdr[0].header)
    date[i] = hdr[0].header['OBSMJD']
    mean[i], median[i], std[i] = sigma_clipped_stats(data, sigma=3.0)
    bkg_aper = SkyCircularAnnulus(pos, 5*u.arcsec, 10*u.arcsec)
    bkg_aperstats = ApertureStats(data, bkg_aper, wcs=wcs)
    aper = SkyCircularAperture(pos, 2*u.arcsec)
    aperstats = ApertureStats(data, aper, local_bkg=bkg_aperstats.median, wcs=wcs) 
    peaks[i] = aperstats.sum
    peaks_err[i] = aperstats.std
    flux[i] = hdr[0].header['GAIN']*peaks[i]/hdr[0].header['EXPOSURE']  + tpl_flux# e-/s 
    flux_err[i] = np.sqrt((hdr[0].header['GAIN']*peaks_err[i]/hdr[0].header['EXPOSURE'])**2 + tpl_flux_err**2)

    mag[i] = -2.5*np.log10(flux[i]) + hdr[0].header['MAGZP'] + hdr[0].header['CLRCOEFF']*hdr[0].header['CLRMED'] + hdr[0].header['APCOR1']

    mag_err[i] = 2.5/(flux[i]*np.log(10))*flux_err[i]  


#build table w/ output to run thru qso_fit
#save all template photometry in one file 

t = Table([date, mag, mag_err], names=('date', 'mag', 'mag_err'))
t.sort(keys='date')


#file naming convention based on coordinates?
#for now just make it something to test pipeline 
#file save location?
#NAMING CONVENTION
#DIRECTORY
t.write(PATH + '/tables/measurements/'+name+'_results.fits', overwrite=True) 

tpl_dict = {
    "name": name,
    "ra": args.ra,
    "dec": args.dec,
    "mag": tpl_mag,
    "mag_err": tpl_mag_err
}

with open(PATH + "/tables/tpl_mags.csv", "a", newline="") as fp:
    # Create a writer object
    writer = csv.DictWriter(fp, fieldnames=tpl_dict.keys())

    # Write the header row
    #writer.writeheader()

    # Write the data rows
    writer.writerow(tpl_dict)
    print('Done writing dict to a csv file')