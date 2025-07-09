import subprocess
import argparse
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
import glob as glob
import numpy as np
from astropy.io import fits


#get infor from meta table to retrieve image
def download_images(data):
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

        url = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/'+year+'/'+monthday+'/'+fracday+'/ztf_'+filefracday+'_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_'+imgtypecode+'_q'+qid+'_sciimg.fits'

        cmd = "wget --output-document=../../../diapl2/CurrentData/ztf_" +filefracday +'.fits "' + url + '"'

        subprocess.run(cmd, shell=True)

        #cmd = "mv ztf_"+filefracday+".fits " + "../../CurrentOriginalData"

        #subprocess.run(cmd, shell=True)

def move_files(good_files):
    for i in range(len(good_files)):

        cmd = "mv " + good_files[i] + " " + "../../../diapl2/WorkingDir"

        subprocess.run(cmd, shell=True)

def clean_directory():
    cmd = "rm ../../../diapl2/CurrentData/*.fits"
    subprocess.run(cmd, shell=True)

data = ascii.read("masked_out.tbl")



#download images
download_images(data)

    
files = glob.glob('../../../diapl2/CurrentData/ztf*.fits')

move_files(files)

clean_directory()    
