#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A script to reduce DECam CCD image data for testing purposes.
May want to check the output destinations for each step just in case.

Want to add: get rid of destination args and make it all automatic and self-contained

"""

import os, glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import sep


# Overscan subtraction routine


def overscanSubtraction(filename, extension, destination):
    """ Overscan subtraction routine
    
    Keyword arguments:
    filename -- the path and file of image
    extension -- FITS image extension
    destination -- directory of output image (current working directory/destination directory)
    """
    im_dir = filename
    im_header = fits.getheader(im_dir, ext=extension)
    im_data   = fits.getdata(im_dir, ext=extension)     # Original dimensions: (4150, 1100)
    
    y_overscan = 4095      # pixel where top overscan begins
    x_overscan = 1030      # pixel where side overscan begins
    
    # Median of each row of side overscan pixels
    med = np.median(im_data[0:y_overscan, x_overscan:], axis=1)
    
    # Crop top overscan region
    im_topOverscanCrop = im_data[0:y_overscan, 0:]
    
    # Subtract median overscan value from data row-wise
    im_topOverscan_sub = np.empty((y_overscan, 1100))
    for i in range(len(med)):
        im_topOverscan_sub[i] = im_topOverscanCrop[i] - med[i]
    
    # Crop side overscan region
    im_overscan_sub = im_topOverscan_sub[:,:x_overscan]
    
    fits.writeto(os.getcwd()+ '/'+ destination + '/' + filename.rsplit('/',1)[1][:-3] + '_overscan_' + str(extension) +'.fits', im_overscan_sub, im_header, overwrite=True)
    print('Written in', os.getcwd()+ '/'+ destination + '/' + filename.rsplit('/',1)[1][:-3] + '_overscan_'+ str(extension) +'.fits')
    
    return None


# Make super-bias


def makeSuperbias(filelist, outputfile):
    """ Generate superbias routine
    
    filelist -- list of FITS files
    outputfile -- name of FITS file to be output, including destination directory
    """
    
    n = len(filelist)
    
    frame0 = fits.getdata(filelist[0])
    
    imsize_y, imsize_x = frame0.shape
    
    fits_stack = np.zeros((imsize_y, imsize_x , n))
    
    for ii in range(0, n):
        im = fits.getdata(filelist[ii])
        fits_stack[:,:,ii] = im
        
    med = np.median(fits_stack, axis=2) # take median pixel value of each pixel in stack

    fits.writeto(os.getcwd()+'/'+outputfile, med, fits.getheader(filelist[0]), overwrite=True)
    print('Written in', os.getcwd()+'/'+outputfile)


# Bias subtraction routine


def biasSubtraction(filelist, superbias, destination):
    """ Subtracts super bias from filelist images
    
    filelist -- list of images to be super bias subtracted
    superbias -- super bias FITS file (note extension)
    destination -- directory of output image (current working directory/destination directory)
    """    
    # Super bias data
    sb_data = fits.getdata(superbias)
    
    for im in filelist:
        im_data = fits.getdata(im)
        bias_sub = im_data - sb_data 
        
        fits.writeto(os.getcwd()+ '/'+ destination + '/bs_'+im.rsplit('/',1)[1], bias_sub, fits.getheader(im), overwrite=True)
        print('Written in:', os.getcwd()+'/'+ destination + '/bs_'+im.rsplit('/',1)[1])


# Merge 2 amps together by flipping x-axis


def mergeAmps(im_ex2, im_ex3):
    """ Combines two halves of CCD together from extensions 2 and 3
    
    im_ex2 -- list of image files with extension 2
    im_ex3 -- list of image files with extension 3, corresponding images must be in same sequential order as im_ex2
    
    """
    
    for im in range(len(im_ex3)):
        im_ex3_flipx = np.flip(fits.getdata(im_ex3[im]), axis=1)
        
        extcombined = np.concatenate((fits.getdata(im_ex2[im]), im_ex3_flipx), axis=1)
        
        fits.writeto(os.getcwd()+'/'+im_ex3[im].rsplit('/',1)[0]+'/merged_'+im_ex3[im][:-7].rsplit('/',1)[1]+'.fits', extcombined, fits.getheader(im_ex2[im]), overwrite=True)
        print('Written in:', os.getcwd()+'/'+im_ex3[im].rsplit('/',1)[0]+'/merged_'+im_ex3[im][:-7].rsplit('/',1)[1]+'.fits')

        

# Extract sources from image


def extract(filename, threshold):
    """
    Extract sources using SEP
    
    filename -- path to FITS file OR numpy.Ndarray
    threshold -- pixel value to extract above
    
    Returns: object with different params
    
    """
    if isinstance(filename, (np.ndarray)):
        data = filename.astype(np.float64)
    
    else:  
        data = fits.getdata(filename).astype(np.float64)
    
    sources = sep.extract(data, threshold)
    
    
    #xcoord = [tple[7] for tple in sources]
    #ycoord = [tple[8] for tple in sources]
    #a = [tple[15] for tple in sources]
    #b = [tple[16] for tple in sources]
    #flag = [tple[29] for tple in sources]
    
    return sources


# Generate a DS9 regions file corresponding to given x,y coordinates


def makeRegions(x, y, destination):
    """
    Make ds9 regions file given (x,y) coordinates
    
    x, y -- numpy arrays
    destination -- path to desired destination must include .reg
    
    """
    
    xcoord_reg = list(map(str, x))
    ycoord_reg = list(map(str, y))
    
    ds9reg = open(destination,"w+")
    ds9reg.write('# Region file format: DS9 version 7.3.2 \n')
    ds9reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
    ds9reg.write('IMAGE \n')
    for ix, iy in zip(xcoord_reg, ycoord_reg):
        ds9reg.write('circle '+ix+' '+iy+' 5 \n')
    ds9reg.close()
    
    print('Written in', destination)

# Executes SEP photometry routine with circles
    
def doPhotometry(data, extractobj, radius):
    """ SEP's photometry feature using circles
    
    data -- ndarray object
    extractobj -- output object from sep.extract
    radius -- radius of circle to extract
    """
    
    flux, fluxerr, flag = sep.sum_circle(data, extractobj['x'], extractobj['y'], radius)
    
    return flux, fluxerr, flag

# Computes a list of fluxes from doPhotometry for a list of files at one exposure time

def computeFlux(filelist, threshold, radius, makeRegionFile = True):
    """Computes fluxes of a list of image files by source extraction and photometry
       Uses doPhotometry()
       
       filelist -- a list of FITS images at a single exposure time
       threshold -- float for extraction threshold
       radius -- integer pixel radius of photometric circle for doPhotometry
       
       Output:
       A list of 1D numpy arrays of length len(filelist). Each 1D numpy array contains extracted flux values for each FITS image in filelist.
    """
    
    dat = []
    
    for im in filelist:
        filelist_obj = extract(im, threshold)
        
        if makeRegionFile:
            makeRegions(filelist_obj['x'], filelist_obj['y'], im.rsplit('/',1)[0]+'/'+im.rsplit('/',2)[2][:-5]+'.reg')
        
        data = fits.getdata(im).astype(np.float64)
        flux, fluxerr, flag = doPhotometry(data, filelist_obj, radius)
        
        dat += [flux]
    
    return dat
    
# Double Gaussian profile

def double_gaussian( x, c1, mu1, sigma1, c2, mu2, sigma2 ):
    result =   c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) ) \
          + c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )
    return result

# Compute gain of x-ray exposed image

def computeGain(fluxlist, nbins, fluxrange, plot=True):
    """
    fluxlist -- list of 1D arrays of extracted fluxes from computeFlux() at one exposure time
    nbins -- integer number of bins
    fluxrange -- tuple of range to restrict histogramming
    """
    
    if plot:
        
        nout = []
        binsout = []
        peak1 = []
        peak2 = []
        for i in range(len(fluxlist)):
            (n, bins, _) = plt.hist(fluxlist[i], bins=nbins, range=fluxrange, density=False, histtype='step')
            bin_centers = bins[:-1] + np.diff(bins) / 2
            
            # For log plot, turn 0s into 1
            n[n==0]=1
            popt, _ = curve_fit(double_gaussian, bin_centers, np.log(n), p0=[1,80000.0,100.0,1,88000.0,100.0])
            xfit = np.linspace(bins[0], bins[-1], 10000)
            plt.plot(xfit, np.exp(double_gaussian(xfit, *popt)), label='i='+str(i))
            peak1 += [popt[1]]
            peak2 += [popt[4]]
            
            # Store hist data
            nout += [n]
            binsout += [bins]
        plt.yscale('log')
        plt.xlim(np.mean([peak1])*0.8, 100000)
        plt.xlabel('Flux')
        plt.ylabel('Counts')
        plt.legend()
        plt.show()
        
        peak1 = np.array(peak1)
        peak2 = np.array(peak2)
        
        gain_from_peak1 = peak1 / ((1570.5 + 1573.48)/2)
        gain_from_peak2 = peak2 / 1730.50
        
        mean_gain = np.mean(np.concatenate((gain_from_peak1,gain_from_peak2)))
        
        return mean_gain, gain_from_peak1, gain_from_peak2
    
    else:
        
        peak1 = []
        peak2 = []
        for i in range(len(fluxlist)):
            (n, bins) = np.histogram(fluxlist[i], bins=nbins, range=fluxrange)
            bin_centers = bins[:-1] + np.diff(bins) / 2
            n[n==0]=1
            popt, _ = curve_fit(double_gaussian, bin_centers, np.log(n), p0=[1,80000.0,100.0,1,88000.0,100.0])
            peak1 += [popt[1]]
            peak2 += [popt[4]]
        peak1 = np.array(peak1)
        peak2 = np.array(peak2)
        
        gain_from_peak1 = peak1 / ((1570.5 + 1573.48)/2)
        gain_from_peak2 = peak2 / 1730.50
        
        mean_gain = np.mean(np.concatenate((gain_from_peak1,gain_from_peak2)))
        
        return mean_gain, gain_from_peak1, gain_from_peak2

# Compute read noise
    
def computeRN(filelist, gain, plot=True):
    """
    gain -- float value of measured gain 
    filelist -- list of raw images filenames (.fz)
    """
    
    readNoise_2 = [] # for ext. 2
    readNoise_3 = [] # for ext. 3
    
    for i in range(len(filelist)):
        
        im_dat_2 = fits.getdata(filelist[i], 2)
        im_dat_3 = fits.getdata(filelist[i], 3)
        
        # Overscan region data
        
        im_dat_2_os = im_dat_2[:,1030:]
        im_dat_3_os = im_dat_3[:,1030:]
        
        # Median subtract overscan region
        
        im_dat_2_os_medsub = np.zeros(np.shape(im_dat_2_os))
        im_dat_3_os_medsub = np.zeros(np.shape(im_dat_3_os))
        
        for j in range(len(im_dat_2_os)):
            im_dat_2_os_medsub[j] = im_dat_2_os[j] - np.median(im_dat_2_os, axis=1)[j]
        for k in range(len(im_dat_2_os)):
            im_dat_3_os_medsub[k] = im_dat_3_os[k] - np.median(im_dat_3_os, axis=1)[k]
        
        # Plots
        
        if plot:
        
            fig, (ax1,ax2) = plt.subplots(1,2,figsize=(15,5))
            (n10s_2, bins10s_2, patches10s_2) = ax1.hist(im_dat_2_os_medsub.flatten(), bins=np.arange(min(im_dat_2_os_medsub.flatten()), max(im_dat_2_os_medsub.flatten()) + 1, 1), range=(-800,800), density=True)
            im_dat_2_os_medsub_fitint = np.linspace(-800, 800, len(im_dat_2_os_medsub.flatten()[abs(im_dat_2_os_medsub.flatten()) < 800]))
            mu_2, sigma_2 = stats.norm.fit(im_dat_2_os_medsub.flatten()[abs(im_dat_2_os_medsub.flatten()) < 800])
            im_dat_2_os_medsub_fit = stats.norm.pdf(im_dat_2_os_medsub_fitint, mu_2, sigma_2)
            ax1.plot(im_dat_2_os_medsub_fitint, im_dat_2_os_medsub_fit, label='Gaussian Fit')
            ax1.set_xlim(-800,800)

            (n10s_3, bins10s_3, patches10s_3) = ax2.hist(im_dat_3_os_medsub.flatten(), bins=np.arange(min(im_dat_3_os_medsub.flatten()), max(im_dat_3_os_medsub.flatten()) + 1, 1), range=(-800,800), density=True)
            im_dat_3_os_medsub_fitint = np.linspace(-800, 800, len(im_dat_3_os_medsub.flatten()[abs(im_dat_3_os_medsub.flatten()) < 800]))
            mu_3, sigma_3 = stats.norm.fit(im_dat_3_os_medsub.flatten()[abs(im_dat_3_os_medsub.flatten()) < 800])
            im_dat_3_os_medsub_fit = stats.norm.pdf(im_dat_3_os_medsub_fitint, mu_3, sigma_3)
            ax2.plot(im_dat_3_os_medsub_fitint, im_dat_3_os_medsub_fit, label='Gaussian Fit')
            ax2.set_xlim(-800,800)
        
        
            ax1.legend(loc='upper right', frameon=False)
            ax2.legend(loc='upper right', frameon=False)
            ax1.set_title((str(filelist[i]))+' Ext 2 '+'(histogram normalized)')
            ax2.set_title((str(filelist[i]))+' Ext 3 '+'(histogram normalized)')
            ax1.set_xlabel('pixel value'); ax1.set_ylabel('counts')
            ax2.set_xlabel('pixel value'); ax2.set_ylabel('counts')
            plt.show()
            
            readNoise_2 += [sigma_2/gain]
            readNoise_3 += [sigma_3/gain]
            
        else:
            
            hist_2, bins_2 = np.histogram(im_dat_2_os_medsub.flatten(), bins=np.arange(min(im_dat_2_os_medsub.flatten()), max(im_dat_2_os_medsub.flatten()) + 1, 1), range=(-800,800))
            mu_2, sigma_2 = stats.norm.fit(im_dat_2_os_medsub.flatten()[abs(im_dat_2_os_medsub.flatten()) < 800])
            
            hist_3, bins_3 = np.histogram(im_dat_3_os_medsub.flatten(), bins=np.arange(min(im_dat_3_os_medsub.flatten()), max(im_dat_3_os_medsub.flatten()) + 1, 1), range=(-800,800))
            mu_3, sigma_3 = stats.norm.fit(im_dat_3_os_medsub.flatten()[abs(im_dat_3_os_medsub.flatten()) < 800])
            
            readNoise_2 += [sigma_2/gain]
            readNoise_3 += [sigma_3/gain]
            
    return readNoise_2, readNoise_3