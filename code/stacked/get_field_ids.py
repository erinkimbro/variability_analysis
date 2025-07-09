

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

c = ascii.read('/Users/erinkimbro/Projects/merian_variable_NEW/GH_obj/field/ztf_gh7_catalog_search_results.tbl')

field, ccdid, qid = get_unique_ids(c)
filter = 'zr'

with open('/Users/erinkimbro/Projects/merian_variable_NEW/workspace/data/fields.dat', 'w') as f:
    for i in range(len(field)):
        line = str(field[i]) + ' ' + str(ccdid[i]) + ' ' + str(qid[i]) + ' ' + filter + '\n'
        f.write(line)