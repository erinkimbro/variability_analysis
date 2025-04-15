import subprocess
import argparse
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
import glob as glob
import numpy as np
from astropy.io import fits


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



def remove_bad_images():
#move good files to Current Working directory
#make a cut on seeing, do this from out table -> DONE
#use lowest background frames

    #PATH = '/Users/erinkimbro/Projects/merian_variable'
    files = glob.glob('../CurrentOriginalData/ztf*.fits')

    mean = np.zeros(len(files))
    median = np.zeros(len(files))
    std = np.zeros(len(files))
    seeing = np.zeros(len(files))
    x = np.zeros(len(files))
    y = np.zeros(len(files))


    for i in range(len(files)):
        try:
            hdr = fits.open(files[i])
            data = hdr[0].data
            seeing[i] = hdr[0].header['SEEING']
            mean[i], median[i], std[i] = sigma_clipped_stats(data, sigma=3.0)
            x[i] = hdr[0]['NAXIS1']*hdr[0]['PIXSCALE']
            y[i] = hdr[0]['NAXIS2']*hdr[0]['PIXSCALE']
        except:
            np.delete(seeing, i)
            np.delete(mean, i)
            np.delete(median, i)
            np.delete(std, i)


    bg_mean, bg_median, bg_std = sigma_clipped_stats(median, sigma=3.0, maxiters=None)
    mask = (median<bg_median+3*bg_std) & (median>bg_median-3*bg_std) 

    good_files = np.array(files)[mask]

    return good_files, bg_median, bg_std

def move_files(good_files):
    for i in range(len(good_files)):

        cmd = "mv " + good_files[i] + " " + "../WorkingDir"

        subprocess.run(cmd, shell=True)

def clean_directory():
    cmd = "rm ../CurrentOriginalData/*.fits"
    subprocess.run(cmd, shell=True)


#create actual process for these functions
#maybe move the functions to a new file just for organization sake? 
parser = argparse.ArgumentParser()
parser.add_argument("ra", type=float, help="RA of object in decimal degrees")
parser.add_argument("dec", type=float, help="dec of object in decmial degrees")
args = parser.parse_args()

#get metatable
cmd = get_metatable(args.ra, args.dec)
print(cmd)
subprocess.run(cmd, shell=True)

#filter table, add filter/seeing as positional argument?
data = filter_table('zr', 2, args.ra, args.dec)

if len(data)>=20:

    #download images
    download_images(data, args.ra, args.dec)

    #remove high bg images
    good_files, bg_median, bg_std = remove_bad_images()
    min_thresh = bg_median + 3*bg_std

    move_files(good_files)

    clean_directory()

else:
    if len(data)>0:
        clean_directory()

PATH = '/Users/erinkimbro/Projects/merian_variable'
#files = glob.glob(PATH + '/CurrentWD/*.fits')

#edit fwhm.par
'''with open(PATH + '/WorkingDir/fwhmm.par', 'r') as f:
    lines = f.readlines()

lines[17] = 'MIN_PEAK=    '+str(min_thresh)+'\n'

with open(PATH + '/WorkingDir/fwhmm.par', 'w') as f:
    f.writelines(lines)
    f.close()'''

#edit aga.par
#with open(PATH + '/WorkingDir/fwhm.par', 'r') as f:
#    lines = f.readlines()
#
#lines[19] = 'REFIM="'+str(min_thresh)+'"\n'
#
#with open(PATH + '/WorkingDir/fwhm.par', 'w') as f:
#    f.writelines(lines)
#    f.close()

