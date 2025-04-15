import glob as glob
import os
import pandas as pd
import numpy as np

#files = glob.glob(PATH + '/CurrentWD/*.fits')

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

