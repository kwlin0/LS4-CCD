#!/usr/bin/env python3

import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
from glob import glob
import os
import csv


# Image region excludes 3 edges by 20 pixels to avoid the glowing edges
image_region = {"x0" : 20,
                "x1" : 1030,
                "y0" : 20,
                "y1" : 4076}
overscan_region = {"x0" : 1031,
                   "x1" : 1051,
                   "y0" : 0,
                   "y1" : 4096}

def get_exptime(filename):
    return fits.getheader(filename)['EXPTIME']

def get_ccdinfo(filename):
    return fits.getheader(filename)['CCD_NAME'], fits.getheader(filename)['AMP_NAME'], fits.getheader(filename)['TAP_NAME']

def get_exptime_list(data):
    """
    Get a list of exposure times from a list of FITS image file names
    """
    exptime_list = list()
    for i in range(len(data.flatten())):
        exptime_list += [get_exptime(data.flatten()[i])]
    exptime_list = np.array(sorted(list(set(exptime_list)))) # literal raw exposure time list
    exptime_ints = np.round(exptime_list)[::2] # int exposure time
    return exptime_ints

def group_pairs(data):
    exptimes = get_exptime_list(data)
    files_by_exptime = np.empty((len(exptimes), 2), dtype='<U65')
    i = 0
    for file in range(len(data.flatten())):
        if i == 2:
            i = 0
        for exp in range(len(exptimes)):
            if int(get_exptime(data.flatten()[file])) == exptimes[exp]:
                files_by_exptime[exp][i] = data.flatten()[file]
                i = i + 1
    return files_by_exptime

def get_ptcdata(data, image_region, overscan_region):
    """
    image_region, overscan_region = dict definitions of pixel coordinates for the active (image)
    area and the overscan (bias) area
    """
    # paired_files = list of tuple filenames of the same exptime
    paired_files = group_pairs(data)
    ccd_name, amp_name, tap_name = get_ccdinfo(data[0])
    signal_mean = np.empty((len(paired_files)))
    signal_mean_diff = np.empty((len(paired_files)))
    signal_var = np.empty((len(paired_files)))
    
    
    for i, name in enumerate(paired_files):
        if amp_name == 'RIGHT':
            im1 = np.flip(fits.getdata(paired_files[i,0]), axis=1)
            im2 = np.flip(fits.getdata(paired_files[i,1]), axis=1)
        else:
            im1 = fits.getdata(paired_files[i,0])
            im2 = fits.getdata(paired_files[i,1])

        im1_active = (im1[image_region["y0"]:image_region["y1"], image_region["x0"]:image_region["x1"]])
        im1_overscan = im1[overscan_region["y0"]:overscan_region["y1"], overscan_region["x0"]:overscan_region["x1"]]
    
        im2_active = (im2[image_region["y0"]:image_region["y1"], image_region["x0"]:image_region["x1"]])
        im2_overscan = im2[overscan_region["y0"]:overscan_region["y1"], overscan_region["x0"]:overscan_region["x1"]]
    
        signal_mean[i] = (np.mean(im1_active) + np.mean(im2_active)) / 2
    
        im1_active = im1_active - np.median(im1_active)
        im2_active = im2_active - np.median(im2_active)
        diff = im1_active - im2_active
        signal_var[i] = np.var(diff.flatten())/2
        signal_mean_diff[i] = np.mean(diff)
    return signal_mean, signal_var

def get_ptcplot(data, image_region, overscan_region, plot=False):
    """
    Plots the PTC and returns the slope (gain factor)
    """
    signal_mean, signal_var = get_ptcdata(data, image_region, overscan_region)
    ccd_name, amp_name, tap_name = get_ccdinfo(data[0])
    df = pd.DataFrame(zip(signal_var[1:], signal_mean[1:]))
    df = df[(np.abs(stats.zscore(df)) < 2).all(axis=1)] # detects outliers with zscore > 3 which is then ignored from fitting
    
    coef = np.polyfit(signal_mean[1:11],signal_var[1:11],1)
    #coef = np.polyfit(df[1], df[0], 1) # fit to the 26th exposure
    poly1d_fn = np.poly1d(coef) 
    if plot:
        plt.figure(figsize=(10,7))
        plt.scatter(signal_mean[1:], signal_var[1:], s=15, c='k')
        plt.plot(signal_mean[1:25], poly1d_fn(signal_mean[1:25]), c='r', label='linear fit')
        plt.scatter(signal_mean[26], signal_var[26], c='b', s=15)
        #plt.axvline(signal_mean[26], ls='--', c='b', label= '%.0f' % signal_mean[26] + ' ADU')
        plt.xlabel('Signal mean [ADU]', size=23)
        plt.ylabel('Signal variance [ADU$^2$]', size=23)
        #plt.xscale('log')
        #plt.yscale('log')
        plt.ylim(0, 10000)
        plt.title(str(ccd_name) + '/' + str(amp_name) +' PTC' + str(poly1d_fn), size=23)
        plt.legend(prop = {'size': 15})
        plt.grid()
        plt.tick_params(direction='in',which='both',top=True,right=True, width=1, labelsize=16)
        plt.tight_layout()
        plt.show()
    return poly1d_fn[1]

def do_ptc(data, image_region, overscan_region, plot=False):
    ccd_name, amp_name, tap_name = get_ccdinfo(data[0])
    gain = get_ptcplot(data, image_region, overscan_region, plot)
    return ccd_name, amp_name, tap_name, gain

def get_std(data, overscan_region):
    """
    Routine should be run on overscan-subtracted images (Gaussian centered at 0)
    """
    overscan = fits.getdata(data)[overscan_region["y0"]:overscan_region["y1"], overscan_region["x0"]:overscan_region["x1"]]
    overscan_1d = overscan.flatten()
    mu, sigma = stats.norm.fit(overscan_1d[abs(overscan_1d) < 30])
    return sigma


def main():
    parser = argparse.ArgumentParser(description="Get PTCs for data ramp")
    parser.add_argument("path", nargs='+', help="FITS files input")
    parser.add_argument('--write', nargs=1, default='ptc_results.csv', type=str, help='Write CSV output file')
    parser.add_argument('-n', '--noise', action='store_true', help='Calculate overscan noise')
    parser.add_argument('-p', '--plot', action='store_true', help='Generate PTC plot (not implemented yet)')
    args = parser.parse_args()
    
    print('Going to write ' + str(args.write))
    #directory = '20240930/gain_sequence/*_???_*.fits'
    
    files = args.path

    # Sort by (1) exposure number, (2) controller number, (3) image number associated with controller
    sfiles = sorted(files,key=lambda x: (int(os.path.splitext(x)[0].split('_')[-2]), (os.path.splitext(x)[0].split('_')[-3]), int(os.path.splitext(x)[0].split('_')[-1]) ) )

    nout = 64 # number of FITS images that are output per exposure
    nexp = int(len(np.array(sfiles)) / nout) # total number of exposures

    # Reshape into (nexp, nout) array
    sfiles_shape = np.array(sfiles).reshape(nexp, nout) # sfiles_shape[:,AMPNUMBER] where AMPNUMBER is an integer out of 64 total amplifiers on the FPA
    
    ccd = list()
    amp = list()
    tap = list()
    gain = list()
    std = list()
    noise = list()
    
    for i in range(sfiles_shape.shape[1]):
        data = sfiles_shape[:,i][1:] # list of exposures from the same CCD
        ccd_name, amp_name, tap_name, measuredgain = do_ptc(data, image_region, overscan_region)
        ccd += [ccd_name]
        amp += [amp_name]
        tap += [tap_name]
        gain += [measuredgain]
        
        printout = [ccd_name, amp_name, tap_name, measuredgain]
        
        if args.noise:
            stdev = get_std(data[0], overscan_region) # sample the noise from the CCD using one of the exposures in the data list, e.g. data[0]
            std += [stdev]
            noise += [stdev/measuredgain]
            printout += [stdev, stdev/measuredgain]
        
        print(printout)
        
        if args.noise:
            df = pd.DataFrame({'CCD_NAME': ccd, 'AMP_NAME': amp, 'TAP_NAME': tap, 'GAIN': gain, 'STDEV': std, 'NOISE': noise})
            df.to_csv(args.write, index=False, header=True)
        else:
            df = pd.DataFrame({'CCD_NAME': ccd, 'AMP_NAME': amp, 'TAP_NAME': tap, 'GAIN': gain})
            df.to_csv(args.write, index=False, header=True)
    
if __name__ == '__main__':
    main()