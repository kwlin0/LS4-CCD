#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Identifying and creating a mask for bad columns and pixels.

Version 1.2

Kenneth Lin, 01/22

"""

# in-built
import os, glob

# standard

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import stats
from scipy.optimize import curve_fit

# specialty

from astropy.io import fits

##################################################


def gaussian(x, c1, mu1, sigma1):
    result =   c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )
    return result

def file_prop(filelist):
    """
    Arguments
    -----------
    filelist: list of FITS files
    
    Returns
    -----------
    nlist: length of a list of 2D arrays
    imsize_y, imsize_x: shape (y dim, x dim) of first element in the list
    
    """
    
    nlist = len(filelist)
    
    try:
        frame0 = fits.getdata(filelist[0])
    except:
        print('Invalid filelist input!')
    
    imsize_y, imsize_x = frame0.shape
    
    return nlist, imsize_y, imsize_x

def stackImages(filelist):
    """
    Create 3D array stack of images.
    (for the purposes of identifying bad pixels, they should be overscan subtracted)
    
    Arguments
    -----------
    filelist: list of FITS files
    
    Returns
    -----------
    fits_stack: 3D numpy array
    
    """
    
    nlist, imsize_y, imsize_x = file_prop(filelist)
    
    fits_stack = np.zeros((imsize_y, imsize_x , nlist))
    
    for ii in range(0, nlist):
        im = fits.getdata(filelist[ii])
        fits_stack[:,:,ii] = im
    
    return fits_stack

def medianImage(fits_stack_array):
    """
    Create median image from stacked images (3D array)
    
    Arguments
    -----------
    fits_stack_array: 3D numpy array
    
    Returns
    -----------
    array: 2D numpy array
    
    """
    return np.median(fits_stack_array, axis=2)

def sigImage(filelist, nlist, imsize_y, imsize_x):
    """
    Create sigma image.
    
    Arguments
    -----------
    filelist: list of FITS files
    nlist: number of elements in filelist (int)
    imsize_y: y-dimension of image 2D array (int)
    imsize_x: x-dimension of image 2D array (int)
    
    Returns
    -----------
    sig: 3D numpy array
    
    """
    print('Stacking images...')
    stackedImages = stackImages(filelist)
    
    print('Creating median image...')
    med = medianImage(stackedImages)
    
    print('Creating sig image...')
    sig = np.zeros((imsize_y, imsize_x , nlist))
    
    for ii in range(0, nlist):
        im = fits.getdata(filelist[ii])
        sig[:,:,ii] = im - med
    
    return sig

def varianceImage(filelist, write=False):
    """
    Create variance image. "Width image" also exists (what was it for?)
    
    Arguments
    -----------
    filelist: list of FITS files (which should be overscan subtracted)
    
    Returns
    -----------
    varim: variance image (2D numpy array)
    
    """
    
    nlist, imsize_y, imsize_x = file_prop(filelist)
    
    ms = int(nlist * 0.1585)
    
    varim = np.zeros((imsize_y, imsize_x))
    wdtim = np.zeros((imsize_y, imsize_x))
    
    sig = sigImage(filelist, nlist, imsize_y, imsize_x)
    
    print('Creating variance image...')
    for j in range(imsize_y):
        for i in range(imsize_x):

            sig_sort = np.sort(sig[j,i,:])
            
            var1 = abs(sig_sort[int(nlist/2)] - sig_sort[ms])
            var2 = abs(sig_sort[int(nlist/2)] - sig_sort[nlist-ms])
            
            var = max(var1,var2)
            # wdt = abs(var1-var2)
            
            varim[j,i] = var
            # wdtim[j,i] = wdt
    
    print('Done.')
    
    if write == True:
        outputname = 'varim.fits'
        fits.writeto(outputname, varim, overwrite=True)
        print('Written as', outputname)
    
    return varim

def displayImage(image):
    """
    A convenience function- display image with imshow
    
    Arguments
    -----------
    image: 2D numpy array
    
    Returns
    -----------
    None
    
    """
    
    fig, (ax1) = plt.subplots(1,1, figsize=(20,30))
    obj1 = ax1.imshow(image, vmin=np.mean(image)/4, vmax=np.mean(image)*2, origin='lower', cmap='gray')
    #obj2 = ax2.imshow(wdtim, vmin=0, vmax=10, origin='lower', cmap='gray')
    fig.colorbar(obj1, ax=ax1)
    #fig.colorbar(obj2, ax=ax2)
    #ax1.set_title(str(image))
    #ax2.set_title('wdt')
    fig.show()

    return None

def findClip(mu, sigma, signum):
    """ Given Gaussian fit mean and sigma, compute pixel value at which to clip above."""
    
    sigma = abs(sigma)
    
    sigmac = signum * sigma
    
    pxclipval = mu + sigmac
    
    return pxclipval

def findClipPx(varim, mu, sigma, signum=5.):
    """ Find indices of pixels with higher value than signum * sigma, where signum is by default 5.
    
    Arguments
    -----------
    varim: variance image array (2D numpy array)
    mu: Fitted variance image Gaussian mean (float)
    sigma: Fitted variance image Gaussian standard deviation (float)
    signum: Multiple of sigma (int or float) default = 5.
    
    Returns
    -----------
    pxclip indices, pxclip value (tuple)
    
    """
    
    pxclipval = findClip(mu, sigma, signum)
    pxclipidx = np.where(varim.flatten() > pxclipval)[0]
    
    return pxclipidx, pxclipval

def findPx(varim, pxclipval):
    """ Returns (x,y) of identified bad/hot pixels"""
    
    y = np.where(varim > pxclipval)[0]
    x = np.where(varim > pxclipval)[1]

    return x, y

def makeMask(varim, x, y, display=False):
    """ Returns mask image array
    
    Arguments
    -----------
    varim: 2D numpy array
        (nominally the variance image array)
    x, y: array_like
        pixel positions
        can be used as inputs to makeRegions for ds9 .reg generation
    display: boolean, optional
        imshow mask, default=False
    
    Returns
    -----------
    mask: 2D numpy array
    
    """
    
    imsize_y, imsize_x = varim.shape
    
    mask = np.zeros((imsize_y, imsize_x))
    
    mask[y,x] = 1.
    
    outputname = 'mask.fits'
    fits.writeto(outputname, mask, overwrite=True)
    print('Written in ', outputname)
    
    if display:
        displayImage(mask)
    
    return mask

# Pixel distribution fitting

def fitPx(varim, binnum=300, hrange=(0,50000), xlim=(-1000,30000)):
    """ Fits a Gaussian over normally-distributed pixels and uses findClipPx to identify unwanted pixels.
    
    Arguments
    -----------
    varim: 2D numpy array
    binnum: float
        Number of histogram bins
    hrange: tuple
        range parameter of histogram
    xlim: tuple
        plot x-axis limits
    
    Returns
    -----------
    pxclipval: float
        Pixel value cutoff
    mu_fit: float
        Gaussian fitting mu
    sigma_fit: float
        Gaussian fitting sigma
    
    """
    
    varim = varim.flatten()
    fitX = np.linspace(-2000, 10000, 800) # array range for curve fit
    
    plt.figure(figsize=(15,5))
    
    # Plot variance image distribution
    
    (n, bins, _) = plt.hist(varim, bins=binnum, range=hrange, histtype='stepfilled', 
                            density=False, facecolor='g', alpha=0.7, label='Variance Image')
    
    bin_centers = bins[:-1] + np.diff(bins) / 2
    
    # The main event
    
    mu_guess       = bins[np.argmax(n)] # Uses peak of histogram as initial guess (should be robust)
    sigma_guess    = 100
    constant_guess = 1
    
    popt, _ = curve_fit(gaussian, bin_centers, n, p0=[constant_guess, mu_guess, sigma_guess])
    
    mu_fit    = popt[1]
    sigma_fit = abs(popt[2])
    
    plt.plot(fitX, gaussian(fitX, *popt), label='Gaussian Fit')
    plt.axvline(mu_fit, ls=':', c='k', label='$x_0$ = '+str('%s' % float('%.2f' % mu_fit)))
    
    pxclip, pxclipval = findClipPx(varim, mu_fit, sigma_fit) # Finds pixel value at which to cutoff and its index
    
    plt.hist(varim[pxclip], bins=binnum, range=hrange, histtype='step', 
             edgecolor='r', hatch='//', label='Masked Px')
    
    #plt.hist(varim, bins=binnum, range=(0,50000), histtype='stepfilled', 
    #                     density=False, facecolor='k', alpha=0.1, label='Variance: Tail')
    plt.yscale('log')
    plt.xlabel('Pixel value', size=15)
    plt.ylabel('N', size=15)
    plt.legend(frameon=False, prop={'size': 15})
    plt.tight_layout()
    plt.xlim(xlim)
    plt.ylim(1e0,np.max(n)*5)
    plt.tick_params(direction='in', axis='both', which='both', labelsize=15)
    plt.grid()
    plt.show()
    
    return pxclipval, mu_fit, sigma_fit