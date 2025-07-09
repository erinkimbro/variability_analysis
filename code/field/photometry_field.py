import glob as glob
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats, CircularAperture, CircularAnnulus, aperture_photometry
import numpy as np
from astropy.table import Table, vstack
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from photutils.aperture import SkyCircularAperture, SkyCircularAnnulus
import csv
from astropy.wcs import WCS
from photutils.utils import calc_total_error
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources
from photutils.segmentation import detect_sources
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.utils import calc_total_error, CutoutImage
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D, MedianBackground


parser = argparse.ArgumentParser()
parser.add_argument("field", type=int, help="field number")
parser.add_argument("ccdid", type=int, help="ccdid")
parser.add_argument("qid", type=int, help="qid")
parser.add_argument("filter", type=str, help="filter")
args = parser.parse_args()

#get files
PATH = '../../..'
PATH2 = '../../../results_'+args.filter+'/'+ str(int(args.field)) + '_' + str(int(args.ccdid)) + '_' + str(int(args.qid))
dir = ['/diapl2/WorkingDir/1_1/','/diapl2/WorkingDir/1_2/','/diapl2/WorkingDir/2_1/','/diapl2/WorkingDir/2_2/']

for z in range(len(dir)):
    
    files = ascii.read(PATH + dir[z]+'simages.txt', names=['file'])

    tpl_file = glob.glob(PATH + dir[z]+'tpl.fits')



    hdr = fits.open(tpl_file[0])
    wcs = WCS(hdr[0].header)
    tpl_data = hdr[0].data
    tpl_hdr = hdr[0].header

    #photometry
    mean = np.zeros(len(files))
    median = np.zeros(len(files))
    std = np.zeros(len(files))
    peaks = np.zeros(len(files))
    peaks_err = np.zeros(len(files))

    seeing = np.zeros(len(files))
    date = np.zeros(len(files))
    flux = []
    flux_err = []
    mag = []
    mag_err = []

    bkg_estimator = MedianBackground()
    bkg = Background2D(tpl_data, (160, 160), filter_size=(3, 3),
                    bkg_estimator=bkg_estimator)
    tpl_data -= bkg.background  # subtract the background

    threshold = 1.5 * bkg.background_rms

 
    segment_map = detect_sources(tpl_data, threshold, npixels=10)

    segm_deblend = deblend_sources(tpl_data, segment_map,
                                npixels=10, nlevels=32, contrast=0.001,
                                progress_bar=False)

    wcs = WCS(hdr[0].header)
    effective_gain = 6.2
    error = calc_total_error(tpl_data, bkg.background_rms, effective_gain)
    cat = SourceCatalog(tpl_data, segm_deblend, background=bkg.background, wcs=wcs, error=error)


    #labels = [1, 2]
    #cat_subset = cat.get_labels(labels)
    columns = ['label', 'xcentroid', 'ycentroid', 'sky_centroid', 'equivalent_radius','segment_flux','segment_fluxerr', 'segment_fluxerr', 'background_centroid', 'background_mean',
            'background_sum']
    tbl4 = cat.to_table(columns=columns)

    pos = tbl4['sky_centroid']

    #hdul = fits.open('/Users/erinkimbro/Projects/merian_variable/data/nsa_v0_1_2.fits')
    hdul = fits.open(PATH + '/workspace/data/merian_dwarfs.fits')
    merian_data = hdul[1].data
    ra = merian_data['coord_ra_Merian']*u.deg
    dec = merian_data['coord_dec_Merian']*u.deg
    catalog = SkyCoord(ra=ra, dec=dec)

    max_sep = 1.0 * u.arcsec #how much separation? base it on astrometry error?
    idx, d2d, d3d = pos.match_to_catalog_3d(catalog)
    sep_constraint = d2d < max_sep
    pos = pos[sep_constraint]
    catalog_matches_merian = merian_data[idx[sep_constraint]]
    tbl4_merian = tbl4[sep_constraint]

    hdul = fits.open(PATH + '/workspace/data/nsa_v0_1_2.fits')
    nsa_data = hdul[1].data
    ra = nsa_data['RA']*u.deg
    dec = nsa_data['DEC']*u.deg
    catalog = SkyCoord(ra=ra, dec=dec)
    pos = tbl4['sky_centroid']

    max_sep = 1.0 * u.arcsec #how much separation? base it on astrometry error?
    idx, d2d, d3d = pos.match_to_catalog_3d(catalog)
    sep_constraint = d2d < max_sep
    pos = pos[sep_constraint]
    catalog_matches_nsa = nsa_data[idx[sep_constraint]]
    tbl4_nsa = tbl4[sep_constraint]

    gh7_data = ascii.read(PATH + '/workspace/data/GH7.csv')
    ra = gh7_data['ra']*u.deg
    dec = gh7_data['dec']*u.deg
    catalog = SkyCoord(ra=ra, dec=dec)
    pos = tbl4['sky_centroid']

    max_sep = 1.0 * u.arcsec #how much separation? base it on astrometry error?
    idx, d2d, d3d = pos.match_to_catalog_3d(catalog)
    sep_constraint = d2d < max_sep
    pos = pos[sep_constraint]
    catalog_matches_gh7 = gh7_data[idx[sep_constraint]]
    tbl4_gh7 = tbl4[sep_constraint]

    if (len(tbl4_nsa)>0) & (len(tbl4_merian)>0):
        tbl4_nsa.add_column(np.repeat('nsa', len(tbl4_nsa)), name='catalog')
        tbl4_merian.add_column(np.repeat('merian', len(tbl4_merian)), name='catalog')
        tbl4 = vstack([tbl4_merian, tbl4_nsa])

    if (len(tbl4_nsa)>0) & (len(tbl4_merian)==0): 
        tbl4_nsa.add_column(np.repeat('nsa', len(tbl4_nsa)), name='catalog')
        tbl4 = tbl4_nsa


    if (len(tbl4_nsa)==0) & (len(tbl4_merian)>0):
        tbl4_merian.add_column(np.repeat('merian', len(tbl4_merian)), name='catalog')
        tbl4 = tbl4_merian

    if (len(tbl4_gh7) > 0) & ((len(tbl4_nsa)>0 | len(tbl4_merian)>0)) :
        tbl4_gh7.add_column(np.repeat('misc', len(tbl4_gh7)), name='catalog')
        tbl4 = vstack([tbl4, tbl4_gh7])

    if (len(tbl4_gh7) > 0) & ((len(tbl4_nsa)==0 & len(tbl4_merian)==0)) :
        tbl4_gh7.add_column(np.repeat('misc', len(tbl4_gh7)), name='catalog')
        tbl4 =  tbl4_gh7

    if (len(tbl4_gh7) == 0) & ((len(tbl4_nsa)==0 & len(tbl4_merian)==0)) :
        tbl4 = []


    
    if len(tbl4)>0:


        segment_mag = -2.5*np.log10(tbl4['segment_flux']) + tpl_hdr['MAGZP']
        tbl4.add_column(segment_mag, name='segment_mag')
        tbl4.write(PATH2 + '/tables/init_tbl_'+str(args.field) + str(args.ccdid) + str(args.qid)+ str(int(z))+'.csv', overwrite=True)
        pos = tbl4['sky_centroid']


        for i in range(len(tbl4)): 
            cutout = CutoutImage(tpl_data, (tbl4['ycentroid'][i], tbl4['xcentroid'][i]), (20, 20))
            norm = ImageNormalize(stretch=LogStretch())
            fig, ax = plt.subplots()
            ax.imshow(cutout, norm=norm)
            plt.savefig(PATH2+'/plots/cutouts/cutout'+tbl4['catalog'][i]+str(int(args.field)) + str(int(args.ccdid)) + str(int(args.qid)) + str(z) + str(i)+'.png')
            plt.close()



        tpl_bkg_aper = SkyCircularAnnulus(tbl4['sky_centroid'], 2.5*u.arcsec, 5*u.arcsec)
        tpl_bkg_aperstats = ApertureStats(tpl_data, tpl_bkg_aper, wcs=wcs, error=error)
        aper = SkyCircularAperture(tbl4['sky_centroid'], 2.25*u.arcsec)
        aperstats = ApertureStats(tpl_data, aper,  local_bkg=tpl_bkg_aperstats.median, wcs=wcs, error = error)
        tpl_peak = aperstats.sum
        tpl_flux_uncor = hdr[0].header['GAIN']*tpl_peak/hdr[0].header['EXPOSURE'] #e-/s

        tpl_flux = tpl_flux_uncor * (10**(-hdr[0].header['MAGZP']/2.5)) * (10**(-hdr[0].header['APCOR3']/2.5))

        exp = hdr[0].header['EXPOSURE'] #s
        dc = hdr[0].header['DARKCUR'] #e-/s
        rn = hdr[0].header['READNOI'] #e-, time independent
        n = aperstats.sum_aper_area
        meas = tpl_flux_uncor * exp #e-
        sky = bkg.background_median * hdr[0].header['GAIN'] #e-

        tpl_flux_err_ccd = np.sqrt(meas + n.value * (sky + dc*exp + rn**2)) / exp

        a = 10**(-hdr[0].header['MAGZP']/2.5)
        b = 10**(-hdr[0].header['APCOR3']/2.5)
        a_err = np.sqrt(((- np.log(10) * a / 2.5 ) * hdr[0].header['MAGZPUNC'])**2)
        b_err = np.sqrt(((- np.log(10) * b / 2.5 ) * hdr[0].header['APCORUN3'])**2)
        tpl_flux_err = np.sqrt((a*b*tpl_flux_err_ccd)**2 + (a*tpl_flux_uncor * b_err)**2 + (b*tpl_flux_uncor*a_err)**2)

        tpl_mag = -2.5*np.log10(tpl_flux) 
        tpl_mag_err =  np.sqrt(((2.5*tpl_flux_err)/(tpl_flux*np.log(10)))**2 )
            #tpl_mag = -2.5*np.log10(tpl_flux) + hdr[0].header['MAGZP'] + hdr[0].header['CLRCOEFF']*hdr[0].header['CLRMED'] + hdr[0].header['APCOR1']
            #tpl_mag_err =  2.5/(tpl_flux*np.log(10))*tpl_flux_err

        tbl4.add_columns([tpl_flux_uncor, tpl_flux_err_ccd, tpl_flux, tpl_flux_err], names=['nuc_flux', 'nuc_flux_err','nuc_flux_cor', 'nuc_flux_err_cor' ])


        '''
        aper_array = np.array([1, 1.5, 2, 3, 5, 7])
        aper_diff = np.abs(np.subtract(aper_array.reshape(1, 6), tbl4['equivalent_radius'].value.reshape(len(tbl4), 1)))
        ind = np.argmin(aper_diff, axis=1)
        aper_corr = np.array(['APCOR1', 'APCOR2', 'APCOR3', 'APCOR4', 'APCOR5', 'APCOR6'])[ind]

        for i in range(len(tbl4)):
            tpl_seg_mag[i] = -2.5*np.log10(tbl4['segment_flux'][i]) + hdr[0].header['MAGZP'] + hdr[0].header['CLRCOEFF']*hdr[0].header['CLRMED'] + hdr[0].header[aper_corr[i]]
            
        tpl_seg_mag_err =  2.5/(tbl4['segment_flux']*np.log(10))*tbl4['segment_fluxerr']
        '''


        seeing = np.zeros(len(files))
        date = np.zeros(len(files))
        aper_flux = np.zeros((len(tbl4), len(files)))
        aper_flux_err = np.zeros((len(tbl4), len(files)))
        aper_mag = np.zeros((len(tbl4), len(files)))
        aper_mag_err = np.zeros((len(tbl4), len(files)))


        for i in range(len(files)):
            hdr = fits.open(PATH + dir[z] + 's_' + files['file'][i])
            data = hdr[0].data
            if np.sum(data) > 0:

                wcs = WCS(hdr[0].header)
                date[i] = hdr[0].header['OBSMJD']

                bhdr = fits.open(PATH + dir[z] + files['file'][i])
                bdata = bhdr[0].data

                bkg_estimator = MedianBackground()
                bbkg = Background2D(data, (160, 160), filter_size=(3, 3),
                            bkg_estimator=bkg_estimator)
                #bdata -= bbkg.background  # subtract the background


                #bkg_estimator = MedianBackground()
                #bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                #               bkg_estimator=bkg_estimator)
                #data -= bkg.background  # subtract the background
                bkg_aper = SkyCircularAnnulus(tbl4['sky_centroid'], r_in=2.5*u.arcsecond, r_out=5*u.arcsecond)
                bkg_aperstats = ApertureStats(data, bkg_aper, wcs=wcs)
                aper = SkyCircularAperture(tbl4['sky_centroid'], r=2.25*u.arcsecond)
                aperstats = ApertureStats(data, aper, wcs=wcs)
                peaks = aperstats.sum - (bkg_aperstats.mean * aperstats.sum_aper_area.value)
                flux1_uncor = hdr[0].header['GAIN']*peaks/hdr[0].header['EXPOSURE']
                flux1 = flux1_uncor * (10**(-hdr[0].header['MAGZP']/2.5)) * (10**(-hdr[0].header['APCOR3']/2.5))
            

                exp = hdr[0].header['EXPOSURE'] #s
                dc = hdr[0].header['DARKCUR'] #e-/s
                rn = hdr[0].header['READNOI'] #e-, time independent
                n = 1 * np.pi * (hdr[0].header['PIXSCALE'] ** 2) # number of pixels
                meas = flux1_uncor * exp #e-
                mean, median, std = sigma_clipped_stats(bdata, sigma=3.0)
                bkg = bbkg.background_median * hdr[0].header['GAIN']
                meas = (tbl4['nuc_flux'] + flux1_uncor) * exp
                diff_flux_err = np.sqrt(tbl4['nuc_flux_err']**2 + (np.sqrt((meas + n*(bkg + (dc * exp) + rn**2)))/exp)**2)

                a = 10**(-hdr[0].header['MAGZP']/2.5)
                b = 10**(-hdr[0].header['APCOR3']/2.5)
                a_err = np.sqrt(((- np.log(10) * a / 2.5 ) * hdr[0].header['MAGZPUNC'])**2)
                b_err = np.sqrt(((- np.log(10) * b / 2.5 ) * hdr[0].header['APCORUN3'])**2)
                diff_flux_err_corr = np.sqrt((a*b*diff_flux_err)**2 + (a*flux1_uncor * b_err)**2 + (b*flux1_uncor*a_err)**2)

                tot_flux = tbl4['nuc_flux'] + flux1_uncor
                tot_flux_err = np.sqrt(diff_flux_err**2 + tbl4['nuc_flux_err']**2)

                tot_flux_cor =tbl4['nuc_flux_cor']+flux1

                flux_err_cor = np.sqrt(tbl4['nuc_flux_err_cor']**2 + diff_flux_err_corr**2)
        

                #print(flux1)
                for j in range(len(tbl4)):
                    if tot_flux[j] > 0:

                        tot_mag = -2.5*np.log10(tot_flux_cor[j])
                        tot_mag_err = (2.5*flux_err_cor[j])/(tot_flux_cor[j]*np.log(10)) 
                        aper_flux[j, i] = tot_flux_cor[j]# e-/s 
                        aper_flux_err[j,i] = flux_err_cor[j]
                        aper_mag[j, i] = tot_mag
                        aper_mag_err[j,i] = tot_mag_err
                #build table w/ output to run thru qso_fit
                #save all template photometry in one file 



        nsa_mask = tbl4['catalog']=='nsa'
        merian_mask = tbl4['catalog']=='merian'
        misc_mask = tbl4['catalog']=='misc'

        name = []
        for i in range(len(tbl4)):

            name.append(str(args.field) + str(args.ccdid) + str(args.qid) + str(int(z)) + str(int(i)))
            mask = aper_mag[i] > 0
            t = Table([date[mask], aper_mag[i][mask], aper_mag_err[i][mask]], names=('date', 'mag', 'mag_err'))
            t.sort(keys='date')
            t.write(PATH2 + '/tables/measurements/results_'+tbl4['catalog'][i]+'_'+name[i]+'.fits', overwrite=True) 


        avg_mag = np.nanmean(np.where(aper_mag > 0, aper_mag, np.nan), axis=1)
        avg_std = np.nanstd(np.where(aper_mag > 0, aper_mag, np.nan), axis=1)
        avg_mag_err = np.sqrt(np.nanmean(np.where(aper_mag_err > 0, aper_mag_err, np.nan)**2, axis=1))
        if len(np.array(name)[merian_mask])>0:
            tpl_t_m = Table([np.array(name)[merian_mask], catalog_matches_merian['objectId_Merian'], pos.ra.value[merian_mask], pos.dec.value[merian_mask], catalog_matches_merian['logmass_gaap'], 
                    catalog_matches_merian['logsfr_gaap'], catalog_matches_merian['z500'], avg_mag[merian_mask], avg_mag_err[merian_mask], avg_std[merian_mask],
                    np.repeat(args.field, len(avg_mag[merian_mask])), np.repeat(args.ccdid, len(avg_mag[merian_mask])), np.repeat(args.qid, len(avg_mag[merian_mask]))],
                    names=('name', 'objectId_Merian', 'ra', 'dec', 'logmass', 'logsfr', 'z' ,'mag', 'mag_err', 'mag_std', 'field', 'ccdid', 'qid'))
            tpl_t_m.write(PATH2 + "/tables/tpl_mags_merian_"+str(args.field) + str(args.ccdid) + str(args.qid)+ str(int(z))+".csv", format='csv', overwrite=True)

        if len(np.array(name)[nsa_mask])>0:
            tpl_t_n = Table([np.array(name)[nsa_mask], catalog_matches_nsa['NSAID'], pos.ra.value[nsa_mask], pos.dec.value[nsa_mask], catalog_matches_nsa['MASS'], 
                avg_mag[nsa_mask], avg_mag_err[nsa_mask],avg_std[nsa_mask], np.repeat(args.field, len(avg_mag[nsa_mask])), 
                np.repeat(args.ccdid, len(avg_mag[nsa_mask])), np.repeat(args.qid, len(avg_mag[nsa_mask]))], 
                names=('name', 'objectId_Merian', 'ra', 'dec', 'logmass', 'mag', 'mag_err', 'mag_std', 'field', 'ccdid', 'qid'))
            tpl_t_n.write(PATH2 + "/tables/tpl_mags_nsa_"+str(args.field) + str(args.ccdid) + str(args.qid)+ str(int(z))+".csv", format='csv', overwrite=True)

        
        if len(np.array(name)[misc_mask])>0:
            tpl_t_misc = Table([np.array(name)[misc_mask], pos.ra.value[misc_mask], pos.dec.value[misc_mask], 
                avg_mag[misc_mask], avg_mag_err[misc_mask], np.repeat(args.field, len(avg_mag[misc_mask])), 
                np.repeat(args.ccdid, len(avg_mag[misc_mask])), np.repeat(args.qid, len(avg_mag[misc_mask]))], 
                names=('name', 'ra', 'dec', 'mag', 'mag_err', 'field', 'ccdid', 'qid'))

            tpl_t_misc.write(PATH2 + "/tables/tpl_mags_misc_"+str(args.field) + str(args.ccdid) + str(args.qid)+ str(int(z))+".csv", format='csv', overwrite=True)


        #master_tpl = ascii.r ead(PATH + "/tables/tpl_mags.csv")

        
        #all_tpl = vstack([master_tpl, tpl_t])
        #all_tpl.write(PATH + "/tables/tpl_mags.csv", format='csv', overwrite=True)

        '''
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
        '''

    else:
        with open('../../../results/field_log.txt','a') as f:
            f.write(str(int(args.field)) + ' ' + str(int(args.ccdid)) + ' ' + str(int(args.qid)) + ' ' + str(z) + ' No Objects in Field')
