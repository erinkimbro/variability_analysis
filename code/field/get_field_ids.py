

import numpy as np
import os
from astropy.io import ascii

def get_unique_ids(c):
    field = []
    ccdid = []
    qid = []
    funique, counts = np.unique(c['field'], return_counts=True)

    for i in range(len(funique)):
        mask = c['field']==funique[i]
        cunique, counts = np.unique(c['ccdid'][mask], return_counts=True)
        #print(funique[i], cunique)

        for j in range(len(cunique)):
            mask = (c['field']==funique[i]) & (c['ccdid']==cunique[j])
            qunique, counts = np.unique(c['qid'][mask], return_counts=True)
            for k in range(len(qunique)):
                print(funique[i], cunique[j], qunique[k])
                field.append(funique[i])
                ccdid.append(cunique[j])
                qid.append(qunique[k])
    return field, ccdid, qid

file = '/Users/erinkimbro/Projects/merian_variable_NEW/workspace/data/merian_dwarf_ztf_results.tbl'

c = ascii.read(file)

ra_mask = (c['ra_01']>(150.11916667-2))&(c['ra_01']<(150.11916667+2))
dec_mask = (c['dec_01']>(2.20583333-2))&(c['dec_01']<(2.20583333+2))

c = c[ra_mask&dec_mask]

field, ccdid, qid = get_unique_ids(c)
filter = 'zr'

with open('/Users/erinkimbro/Projects/merian_variable_NEW/workspace/data/fields.dat', 'w') as f:
    for i in range(len(field)):
        line = str(field[i]) + ' ' + str(ccdid[i]) + ' ' + str(qid[i]) + ' ' + filter + '\n'
        f.write(line)