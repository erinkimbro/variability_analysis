from astropy.table import Table, vstack, hstack, join
import glob as glob
import numpy as np
import os
import shutil
from astropy.io import fits, ascii
import qso_fit as qf
from astropy.io import fits, ascii
import csv
from astropy.table import Table, vstack, hstack, join
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clipped_stats
import glob as glob



def combine_prop_tables(catalog, filter):
    '''
    combine property tables from all fields
    catalog: string of desired catalog -- merian, nsa, misc
    '''
    PATH = '/Users/erinkimbro/Projects/merian_variable_NEW/results_%s'%(filter)
    files = glob.glob(PATH+'/*/tables/prop_'+catalog+'*')

    tab = []
    for i in range(len(files)):
        tab1 = ascii.read(files[i])
        tab.append(tab1)
    tab = vstack(tab)
    return tab


def duplicates(tab, catalog, filter):
    PATH = '/Users/erinkimbro/Projects/merian_variable_NEW/results_%s'%(filter)


    funique, idx, counts = np.unique(tab['objectId_Merian'],return_index=True, return_counts=True)
    idx_mult = idx[counts>1]
    idx_single = idx[counts==1]
    tab_new = tab[idx_mult]
    tab_new.write(PATH + '/combined/tables/mult_original_prop_'+catalog+'.csv', overwrite=True)
    tab_old = tab[idx_single]
    tab_old.write(PATH + '/combined/tables/single_prop_'+catalog+'.csv', overwrite=True)

    for i in range(len(funique[counts>1])):

        name = tab['name'][tab['objectId_Merian']==funique[counts>1][i]]
        objid = tab['objectId_Merian'][tab['objectId_Merian']==funique[counts>1][i]][0]

        meas_file=[]
        for j in range(len(name)): 
            meas1_file = glob.glob(PATH + '/*/tables/measurements/results_'+catalog+'_'+str(name[j])+'*')
            print(meas1_file)
            meas_file.append(meas1_file)
        
        meas_tab = []
        meas_file = np.concatenate(meas_file)
        if len(meas_file)>1:
            for j in range(len(meas_file)):
                meas1 = fits.open(meas_file[j])
                meas1 = Table(meas1[1].data, names=('date', 'mag', 'magerr'))
                meas_tab.append(meas1)
            
            meas = vstack(meas_tab)

            unique1, idx1, counts1 = np.unique(meas['date'],return_index=True, return_counts=True)
            mag_new = np.zeros(len(unique1))
            mag_err_new = np.zeros(len(unique1))
            date_new = np.zeros(len(unique1))
            for j in range(len(unique1)):
                mag = meas['mag'][meas['date']==unique1[j]]
                mag_err = meas['magerr'][meas['date']==unique1[j]]

                mag_new[j] = np.sum(mag)/len(mag)
                mag_err_new[j] = np.sum(mag_err**2)/len(mag_err)
                date_new[j] = unique1[j]

            #print(mag_new, mag_err_new, date_new)
            new_meas = Table([date_new, mag_new, mag_err_new], names=('date', 'mag', 'magerr'))

            #delete original file or just row from master table
            #create new folder for repeat measurements?
            #rerun qso fit
            #rerun light curve
            #should i do this in a bash script? 
            for j in range(len(name)):
                meas1_file = glob.glob(PATH + '/*/tables/measurements/results_'+catalog+'_'+str(name[j])+'*')
                print(meas1_file)

                try:
                    shutil.copyfile(meas1_file[0], PATH + '/combined/tables/measurements/results_'+catalog+'_'+str(name[j])+'.fits')
                except shutil.SameFileError:
                    pass

            new_meas.write(PATH + '/combined/tables/combined_measurements/results_'+catalog+'_'+str(objid)+'.fits', overwrite=True)



def qsofit(catalog, filter):
    PATH = '/Users/erinkimbro/Projects/merian_variable_NEW/'
    PATH2 = '/Users/erinkimbro/Projects/merian_variable_NEW/results_%s/combined/'%(filter)

    
    all_fits = []

    files = glob.glob(PATH2+'tables/combined_measurements/results_'+catalog+'_*.fits')

    name_all = []

    for i in range(len(files)):

        try:

            #load in table
            #create arguments for file to select files
            diff_table = fits.open(files[i])

            name = files[i].split('_')[-1][:-5]
            name_all.append(name)

            data = diff_table[1].data

            mag_mean, mag_median, mag_std = sigma_clipped_stats(data['mag'], sigma=3.0)
            mask = (data['mag']<mag_mean+3*mag_std)&(data['mag']>mag_mean-3*mag_std)
            data = data[mask]

            fit_model = qf.qso_fit(data['date'], data['mag'], data['magerr'], filter='r', return_model=True )

            fit = qf.qso_fit(data['date'], data['mag'], data['magerr'], filter='r', return_model=False )


            model_t = Table([data['date'], data['mag'], data['magerr'], fit_model['model'], fit_model['dmodel']], 
                            names=('date', 'mag', 'mag_err','model_mag', 'model_mag_err'))

            model_t.write(PATH2 + '/tables/models/model_'+catalog+'_'+name+'.fits', overwrite=True)

            # Open a csv file for writing
            # COME UP WITH NAMING CONVENTION
            # ALSO FIX DIRECTORY
            # SAVE MODEL AND FIT STATISTICS SEPARELATELY
            # APPEND FIT STATISTICS FOR ALL OBJECTS IN ONE FILE
            fit.update({'objectId_Merian':int(name)})
            #fit.update({"ra":ra[i]})
            #fit.update({"dec":dec[i]})

            all_fits.append(fit)
            
            
        except FileNotFoundError:
            pass

        except ValueError:
            pass

    #write file by field

    t_n = Table(rows=all_fits)

    tpl_t_n = ascii.read(PATH2 + "/tables/mult_original_prop_"+catalog+".csv", format='csv')


    prop_t_n = join(t_n, tpl_t_n, keys='objectId_Merian', table_names=['', 'SINGLE_'],

           uniq_col_name='{table_name}{col_name}') 

    prop_t_n.write(PATH2 + "/tables/mult_prop_combined_"+catalog+".csv", format='csv', overwrite=True)

    t_n.write(PATH2 + "/tables/fit_combined_"+catalog+".csv", format='csv', overwrite=True)
    #print(t)


def combine_tables(catalog, filter):
    PATH2 = '/Users/erinkimbro/Projects/merian_variable_NEW/results_%s/combined/'%(filter)
    single_tab = ascii.read(PATH2 + "/tables/single_prop_"+catalog+".csv")
    mult_tab = ascii.read(PATH2 + "/tables/mult_prop_combined_"+catalog+".csv")

    tab = vstack([single_tab, mult_tab], join_type='inner')
    tab.write(PATH2 + "/tables/all_prop_"+catalog+".csv", format='csv', overwrite=True)




filter = 'zg'

#deal with duplicate files
nsa = combine_prop_tables('nsa', filter)
print(len(nsa))
merian = combine_prop_tables('merian', filter)
print(len(merian))

#combine duplicates                
duplicates(merian, 'merian', filter)
duplicates(nsa, 'nsa', filter)

#run qsofit
qsofit('merian', filter)
qsofit('nsa', filter)


#combine properties from single and mult tables 
combine_tables('merian', filter)
combine_tables('nsa', filter)