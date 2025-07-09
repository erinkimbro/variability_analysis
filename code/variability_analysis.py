'''
self contained functions to run variability analysis after photometry is performed
specifically rerun qso_fit
construct light curves 
'''
import glob as glob
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
import numpy as np
import qso_fit as qf 
from astropy.table import Table, vstack, hstack, join
import matplotlib.pyplot as plt
from taufit import simulate_drw
from taufit import fit_drw, fit_carma
import astropy.units as u
import subprocess
import os
import pandas as pd



def get_metatable(RA, DEC):
#get metatable
    cmd = "curl -o out.tbl 'https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?POS="+str(RA)+","+str(DEC)+"&ct=ipac_table'"

    #subprocess.run(cmd, shell=True)
    return cmd


def filter_table(filter, seeing, RA, DEC):
    '''
    filter metatable
    choose filter: zg, zr, zi
    choose seeing limit: CREATE DEFAULT VALUE?
    '''
    data = ascii.read("out.tbl")  
    data = data[data['filtercode']==filter]
    data = data[data['seeing']<seeing]
    ra_mask = (np.abs(RA - data['ra1']) >= 150/3600) & (np.abs(RA - data['ra2']) >= 150/3600) & (np.abs(RA - data['ra3']) >= 150/3600) &(np.abs(RA - data['ra4']) >= 150/3600)
    data = data[ra_mask]
    dec_mask = (np.abs(DEC - data['dec1']) >= 150/3600) & (np.abs(DEC - data['dec2']) >= 150/3600) & (np.abs(DEC - data['dec3']) >= 150/3600) &(np.abs(DEC - data['dec4']) >= 150/3600)
    data = data[dec_mask]
    return data

#get infor from meta table to retrieve image
def download_images(data, RA, DEC):
    for i in range(len(data)):
        filefracday = str(data['filefracday'][i])
        year = str(data['filefracday'][i])[0:4]
        monthday = str(data['filefracday'][i])[4:8]
        fracday = str(data['filefracday'][i])[8:]
        field = str(data['field'][i])
        filtercode = str(data['filtercode'][i])
        ccdid = str(data['ccdid'][i])
        imgtypecode = str(data['imgtypecode'][i])
        qid = str(data['qid'][i])

        paddedfield = ('000000'+field)[-6:]

        paddedccdid = ('00' + ccdid)[-2:]

        url = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/'+year+'/'+monthday+'/'+fracday+'/ztf_'+filefracday+'_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_'+imgtypecode+'_q'+qid+'_sciimg.fits?center='+str(RA)+','+str(DEC)+'&size=300arcsec&gzip=false'

        cmd = "wget --output-document=ztf_" +filefracday +'.fits "' + url + '"'

        subprocess.run(cmd, shell=True)

        cmd = "mv ztf_"+filefracday+".fits " + "../CurrentOriginalData"

        subprocess.run(cmd, shell=True)

def get_image_list():
    files = glob.glob('../../diapl2/WorkingDir/*.fits')

    with open('../../diapl2/WorkingDir/images.txt', 'w') as f:
        for line in files:
            f.write(f"{os.path.basename(line)}\n")

def diapl_setup():
    #SET SATURATION LEVEL
    files = glob.glob('../../../diapl2/WorkingDir/ztf*.fits')

    sat = np.zeros(len(files))
    bg_mean = np.zeros(len(files))
    bg_med = np.zeros(len(files))
    bg_std = np.zeros(len(files))

    for i in range(len(files)):
        hdr = fits.open(files[i])
        data = hdr[0].data
        header = hdr[0].header

        bg_mean[i], bg_med[i], bg_std[i] = sigma_clipped_stats(data, sigma=3.0)

        sat[i] = header['SATURATE']

        #print(np.mean(header['SATURATE']), np.std(header['SATURATE']))

    #edit instrument.par
    with open('../../../diapl2/WorkingDir/instrument.par', 'r') as f:
        lines = f.readlines()

    lines[12] = 'SAT_LEVEL = '+str(3*np.median(sat))+'\n'

    with open('../../../diapl2/WorkingDir/instrument.par', 'w') as f:
        f.writelines(lines)
        f.close()


    #edit fwhmm.par
    with open('../../../diapl2/WorkingDir/fwhmm.par', 'r') as f:
        lines = f.readlines()

    lines[17] = 'MIN_PEAK = '+str(np.median(3*bg_med))+'\n'
    lines[20] = 'MAX_SKY = '+str(np.median(1.5*bg_med))+'\n'



    with open('../../../diapl2/WorkingDir/fwhmm.par', 'w') as f:
        f.writelines(lines)
        f.close()

def set_ref_im():  
    data = pd.read_csv('../../diapl2/WorkingDir/fwhms.txt', sep="\t", names=['file', 'mean fwhm', 'background', 'med fwhm', 'numstars'], index_col=False)
    ind = np.argmin(data['background'])
    ref_name = data['file'][ind]

    #edit diapl_setup.par
    with open('../../diapl2/WorkingDir/diapl_setup.par', 'r') as f:
        lines = f.readlines()

    lines[22] = 'REFIM="'+ref_name+'"\n'

    with open('../../diapl2/WorkingDir/diapl_setup.par', 'w') as f:
        f.writelines(lines)
        f.close()

def run_diapl():

    cmd = "bash fwhms.bash"
    subprocess.run(cmd, shell=True)

    set_ref_im()

    cmd = "bash mktpllist.bash"
    subprocess.run(cmd, shell=True)

    cmd = "bash shifts.bash" 
    subprocess.run(cmd, shell=True)
 
    cmd = "bash template.bash"
    subprocess.run(cmd, shell=True)

    cmd = "bash subtract1.bash"
    subprocess.run(cmd, shell=True)




def load_lc(filter, ztf_id, merian_id):
    '''
    load light curve data from computer, chooses combined field over, original lc
    cuts points with large err or signficant deviation from mean 
    '''
    file_name1 = '/Users/erinkimbro/Projects/merian_variable_NEW/results_%s/*/tables/measurements/results_merian_%s.fits'%(filter, ztf_id)
    file_name2 = '/Users/erinkimbro/Projects/merian_variable_NEW/results_%s/combined/tables/combined_measurements/results_merian_%s.fits'%(filter, merian_id)

    file1 = glob.glob(file_name1)
    file2 = glob.glob(file_name2)
    print(file1, file2)

    if len(file2)>0:
        file = file2
        lc = fits.getdata(file[0])
        x = (lc['date'] - np.min(lc['date']))
        mags = lc['mag']
        mag_err = lc['magerr']

        mean_err, med_err, std_err = sigma_clipped_stats(mag_err, sigma=3, maxiters=5)
        mean, med, std = sigma_clipped_stats(mags, sigma=3, maxiters=5)

        mask = (np.abs(mags-med) < 3*std) & (np.abs(mag_err-med_err)<2*med_err)
        mags = mags[mask]
        mag_err = mag_err[mask]
        x = x[mask]

    else:
        file = file1
        lc = fits.getdata(file[0])
        x = (lc['date'] - np.min(lc['date']))
        mags = lc['mag']
        mag_err = lc['mag_err']


        mean_err, med_err, std_err = sigma_clipped_stats(mag_err, sigma=5, maxiters=5)
        mean, med, std = sigma_clipped_stats(mags, sigma=5, maxiters=5)

        mask = (np.abs(mags-mean) < 5*std) & (np.abs(mag_err-mean_err)<2*med_err)
        mags = mags[mask]
        mag_err = mag_err[mask]
        x = x[mask]
    
    return x, mags, mag_err


def run_qso_fit(x, mags, mag_err, filter):
    
    '''
    runs qso fit
    x: date
    returns table with orignal input data and fitted model 
    returns dict of qso fit statistics 
    '''

    fit_model = qf.qso_fit(x, mags, mag_err, filter=filter, return_model=True )

    fit = qf.qso_fit(x, mags, mag_err, filter=filter, return_model=False )

    model_t = Table([x, mags, mag_err, fit_model['model'], fit_model['dmodel']], 
                    names=('date', 'mag', 'mag_err','model_mag', 'model_mag_err'))
    
    return fit, model_t



def plot_lc(model_t, fit, id, PATH='', savefig=False, ):

    mag = model_t['mag']
    date = model_t['date']
    mag_err = model_t['mag_err']


    fig, ax = plt.subplots()
    ax.errorbar(model_t['date'], model_t['mag'], yerr=model_t['mag_err'], fmt='o', mec='black', color='none', ecolor='grey', markersize=4, label='Data')
    ax.hlines(np.mean(model_t['mag']), np.min(date)-100, np.max(date)+100, color = 'black', label='Mean Mag')
    ax.fill_between(model_t['date'], model_t['model_mag']-model_t['model_mag_err'], model_t['model_mag']+model_t['model_mag_err'], 
                    alpha=0.3, color='orange')
    ax.plot(model_t['date'], model_t['mag'], color='black', label='Model')
    ax.text(.05, .95, r'$\sigma_{vary}$='+str(np.round(fit['signif_vary'], 3)), transform=ax.transAxes)
    ax.text(.05, .90, r'$\sigma_{QSO}$='+str(np.round(fit['signif_qso'], 3)), transform=ax.transAxes)
    ax.text(.05, .85, r'$\sigma_{notQSO}$='+str(np.round(fit['signif_not_qso'], 3)), transform=ax.transAxes)
    #ax.text(.05, .85, r'$M_{*}$='+str(np.round(fit['signif_not_qso'], 3)), transform=ax.transAxes)


    ymax = np.max(model_t['mag']+model_t['mag_err'])+5*np.std(model_t['mag'])
    ymin = np.min(model_t['mag']+model_t['mag_err'])-5*np.std(model_t['mag'])

    ax.set_xlabel('MJD')
    ax.set_ylabel(r'$m_r$')
    ax.set_xlim(np.min(model_t['date'])-50, np.max(model_t['date'])+50)
    ax.set_ylim(ymax, ymin)
    ax.set_title(str(id))
    ax.legend()
    if savefig==True:
        plt.savefig(PATH + str(id)+'_lc.png')
    else:
        pass
    #ax.gca().invert_yaxis()
    
def fit_lc(x, y, yerr, plot=False, verbose=False):
    # Fit DRW
    # Note the inputs to fit_drw must have astropy units
    gp, samples, fig = fit_drw(x*u.day, y*u.mag, yerr*u.mag, plot=plot, verbose=verbose)
    log_tau_drw_recovered = np.log10(1/np.exp(np.median(samples[:,1])))
    return x, y, yerr, gp, samples

