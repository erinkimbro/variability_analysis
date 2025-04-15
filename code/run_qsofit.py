import qso_fit as qf
from astropy.io import fits
import csv
from astropy.table import Table
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clipped_stats


PATH = '/Users/erinkimbro/Projects/merian_variable'


#define arguments
parser = argparse.ArgumentParser()
parser.add_argument("ra", type=float, help="RA of object in decimal degrees")
parser.add_argument("dec", type=float, help="dec of object in decmial degrees")
args = parser.parse_args()

pos = SkyCoord(ra=args.ra, dec=args.dec, unit='deg')
name = pos.to_string('hmsdms')

#load in table
#create arguments for file to select files
diff_table = fits.open('../tables/measurements/'+name+'_results.fits')

data = diff_table[1].data

mag_mean, mag_median, mag_std = sigma_clipped_stats(data['mag'], sigma=3.0)
mask = (data['mag']<mag_mean+3*mag_std)&(data['mag']>mag_mean-3*mag_std)
data = data[mask]

fit_model = qf.qso_fit(data['date'], data['mag'], data['mag_err'], filter='g', return_model=True )

fit = qf.qso_fit(data['date'], data['mag'], data['mag_err'], filter='g', return_model=False )


model_t = Table([data['date'], fit_model['model'], fit_model['dmodel']], names=('date', 'mag', 'mag_err'))

model_t.write(PATH + '/tables/models/'+name+'_model.fits', overwrite=True)

# Open a csv file for writing
# COME UP WITH NAMING CONVENTION
# ALSO FIX DIRECTORY
# SAVE MODEL AND FIT STATISTICS SEPARELATELY
# APPEND FIT STATISTICS FOR ALL OBJECTS IN ONE FILE
fit.update({"name":name})
fit.update({"ra":args.ra})
fit.update({"dec":args.dec})

with open(PATH + "/tables/fit.csv", "a", newline="") as fp:
    # Create a writer object
    writer = csv.DictWriter(fp, fieldnames=fit.keys())

    # Write the header row
    #writer.writeheader()

    # Write the data rows
    writer.writerow(fit)
    print('Done writing dict to a csv file')