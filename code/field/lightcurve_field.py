import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy.coordinates import SkyCoord
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
import glob as glob

parser = argparse.ArgumentParser()
parser.add_argument("field", type=int, help="field number")
parser.add_argument("ccdid", type=int, help="ccdid")
parser.add_argument("qid", type=int, help="qid")
parser.add_argument("filter", type=str, help="filter")
args = parser.parse_args()


PATH = '../../..'
PATH2 = '../../../results_'+args.filter+'/'+ str(int(args.field)) + '_' + str(int(args.ccdid)) + '_' + str(int(args.qid))


#tpl_table = ascii.read(PATH + '/tables/tpl_mags.csv')
#name = tpl_table['name']
#ra = tpl_table['ra']
#dec = tpl_table['dec']

prop_files = glob.glob(PATH2 + '/tables/prop_*'+str(args.field) + str(args.ccdid) + str(args.qid) +'*.csv')


for j in range(len(prop_files)):
 
    prop = ascii.read(prop_files[j])
    catalog = prop_files[j].split('_')[-2]

    for i in range(len(prop)):

        name = str(prop['name'][i])

        hdr = fits.open(PATH2 + '/tables/measurements/results_'+catalog+'_'+name+'.fits')
        t = hdr[1].data

        hdr = fits.open(PATH2 + '/tables/models/model_'+catalog+'_'+name+'.fits')
        model_t = hdr[1].data

        if (len(t)>0)&(len(model_t)>0):

            #need to read in fit table -> search for object
            #mask = tp['name'] == str(args.field) + str(args.ccdid) + str(args.qid) + str(i)
            #fit = fit[mask]

            #need to read in tpl table -> search for object
            #tpl = ascii.read(PATH + '/tables/tpl_mags'+str(args.field) + str(args.ccdid) + str(args.qid) +'.csv')
            #mask = tpl['name'] == name[i]
            #tpl = tpl[mask]

            mag_mean, mag_median, mag_std = sigma_clipped_stats(model_t['mag'], sigma=3.0)

            #mask = (model_t['mag']<mag_mean+5*mag_std)&(model_t['mag']>mag_mean-5*mag_std)

            mag = model_t['mag']
            date = model_t['date']
            mag_err = model_t['mag_err']


            fig, ax = plt.subplots()
            ax.errorbar(model_t['date'], model_t['mag'], yerr=model_t['mag_err'], fmt='o', color='black', ecolor='grey', markersize=4, label='Data')
            ax.hlines(np.mean(model_t['mag']), np.min(date)-100, np.max(date)+100, color = 'black', label='Template Mag')
            ax.fill_between(model_t['date'], model_t['model_mag']-model_t['model_mag_err'], model_t['model_mag']+model_t['model_mag_err'], 
                            alpha=0.3, color='teal')
            ax.plot(model_t['date'], model_t['mag'], color='teal', label='Model')
            ax.text(.05, .95, r'$\sigma_{vary}$='+str(np.round(prop['signif_vary'][i], 3)), transform=ax.transAxes)
            ax.text(.05, .90, r'$\sigma_{QSO}$='+str(np.round(prop['signif_qso'][i], 3)), transform=ax.transAxes)
            ax.text(.05, .85, r'$\sigma_{notQSO}$='+str(np.round(prop['signif_not_qso'][i], 3)), transform=ax.transAxes)

            ax.set_xlabel('MJD')
            ax.set_ylabel(r'$m_r$')
            ax.set_xlim(np.min(model_t['date'])-50, np.max(model_t['date'])+50)
            ax.set_ylim(np.mean(model_t['mag'])+10*np.std(model_t['mag'])+0.5, np.mean(model_t['mag'])-10*np.std(model_t['mag'])-0.5)
            ax.legend()
            #ax.gca().invert_yaxis()
            plt.savefig(PATH2+'/plots/light_curves/'+catalog+'_'+name+'_lightcurve.png') 
            plt.close()
        else:
            pass