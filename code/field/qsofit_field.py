import qso_fit as qf
from astropy.io import fits, ascii
import csv
from astropy.table import Table, vstack, hstack, join
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clipped_stats
import glob as glob


parser = argparse.ArgumentParser()
parser.add_argument("field", type=int, help="field number")
parser.add_argument("ccdid", type=int, help="ccdid")
parser.add_argument("qid", type=int, help="qid")
args = parser.parse_args()


PATH = '../../..'
PATH2 = '../../../results/'+ str(int(args.field)) + '_' + str(int(args.ccdid)) + '_' + str(int(args.qid))


tpl_files = glob.glob(PATH2 + '/tables/tpl_mags_*.csv')

print(tpl_files)



for j in range(len(tpl_files)):

    

    tpl_table = ascii.read(tpl_files[j])

    catalog = tpl_files[j].split('_')[-2]

    all_fits = []

    files = glob.glob(PATH2+'/tables/measurements/results_'+catalog+'_'+tpl_files[j].split('_')[-1][:-4] + '*.fits')



    for i in range(len(files)):

        try:

            #load in table
            #create arguments for file to select files
            diff_table = fits.open(files[i])

            name = files[i].split('_')[-1][:-5]

            data = diff_table[1].data

            mag_mean, mag_median, mag_std = sigma_clipped_stats(data['mag'], sigma=3.0)
            mask = (data['mag']<mag_mean+3*mag_std)&(data['mag']>mag_mean-3*mag_std)
            data = data[mask]

            fit_model = qf.qso_fit(data['date'], data['mag'], data['mag_err'], filter='r', return_model=True )

            fit = qf.qso_fit(data['date'], data['mag'], data['mag_err'], filter='r', return_model=False )


            model_t = Table([data['date'], data['mag'], data['mag_err'], fit_model['model'], fit_model['dmodel']], 
                            names=('date', 'mag', 'mag_err','model_mag', 'model_mag_err'))

            model_t.write(PATH2 + '/tables/models/model_'+catalog+'_'+name+'.fits', overwrite=True)

            # Open a csv file for writing
            # COME UP WITH NAMING CONVENTION
            # ALSO FIX DIRECTORY
            # SAVE MODEL AND FIT STATISTICS SEPARELATELY
            # APPEND FIT STATISTICS FOR ALL OBJECTS IN ONE FILE
            fit.update({"name":int(name)})
            #fit.update({"ra":ra[i]})
            #fit.update({"dec":dec[i]})

            all_fits.append(fit)
            
            
        except FileNotFoundError:
            pass

        except ValueError:
            pass

    #write file by field

    t_n = Table(rows=all_fits)

    tpl_t_n = ascii.read(PATH2 + "/tables/tpl_mags_"+catalog+'_'+tpl_files[j].split('_')[-1][:-4] +".csv", format='csv')


    prop_t_n = join(t_n, tpl_t_n, keys='name')

    prop_t_n.write(PATH2 + "/tables/prop_"+catalog+"_"+tpl_files[j].split('_')[-1][:-4] +".csv", format='csv', overwrite=True)

    t_n.write(PATH2 + "/tables/fit_"+catalog+"_"+tpl_files[j].split('_')[-1][:-4] +".csv", format='csv', overwrite=True)
    #print(t)

 

    #master_fit = ascii.read(PATH + "/tables/fit_all.csv")

    #all_fit_t = vstack([master_fit, t])
    #all_fit_t.write(PATH + "/tables/fit_all.csv", format='csv', overwrite=True)


    '''

    with open(PATH + "/tables/fit_all.csv", "w", newline="") as fp:
        # Create a writer object
        writer = csv.DictWriter(fp, fieldnames=all_fits[0].keys())

        # Write the header row
        writer.writeheader()

        # Write the data rows
        for i in range(len(all_fits)):
            writer.writerow(all_fits[i])
        print('Done writing dict to a csv file')
    '''
