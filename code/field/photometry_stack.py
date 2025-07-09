import glob as glob
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats, CircularAperture, CircularAnnulus, aperture_photometry
import numpy as np
from astropy.table import Table, vstack
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from photutils.aperture import SkyCircularAperture, SkyCircularAnnulus
import csv
from astropy.wcs import WCS
from photutils.utils import calc_total_error
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources
from photutils.segmentation import detect_sources
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.utils import calc_total_error, CutoutImage
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D, MedianBackground

parser = argparse.ArgumentParser()
parser.add_argument("field", type=int, help="field number")
parser.add_argument("ccdid", type=int, help="ccdid")
parser.add_argument("qid", type=int, help="qid")
parser.add_argument("filter", type=str, help="filter")
args = parser.parse_args()

def bin_data(t, data, bin_size):
    ##bin data background subtraction?
    step = np.arange(0, len(t), bin_size)
    data_bin = []
    t_bin = []
    #np.sum(data_sim-np.median(data_sim), axis=0)/len(data_sim)
    for i in range(len(step)-1):
        t_bin.append(np.mean(t[step[i]:step[i+1]]))
        bkg = np.median(data[step[i]:step[i+1]])
        data_bin.append(np.sum(data[step[i]:step[i+1]] - bkg, axis=0)/len(data[step[i]:step[i+1]]))

    return t_bin, data_bin

