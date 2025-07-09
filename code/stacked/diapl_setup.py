#DIAPL SETUP

import glob as glob
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats


#SET SATURATION LEVEL
files = glob.glob('../../../diapl2/WorkingDir/*.fits')

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

