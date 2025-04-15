import glob as glob
import os
import pandas as pd
import numpy as np

PATH = '/Users/erinkimbro/Projects/merian_variable/WorkingDir'
#files = glob.glob(PATH + '/CurrentWD/*.fits')

data = pd.read_csv(PATH+'/fwhms.txt', sep="\t", names=['file', 'mean fwhm', 'background', 'med fwhm', 'numstars'], index_col=False)
ind = np.argmin(data['background'])
ref_name = data['file'][ind]

#edit diapl_setup.par
with open(PATH + '/diapl_setup.par', 'r') as f:
    lines = f.readlines()

lines[22] = 'REFIM="'+ref_name+'"\n'

with open(PATH + '/diapl_setup.par', 'w') as f:
    f.writelines(lines)
    f.close()

