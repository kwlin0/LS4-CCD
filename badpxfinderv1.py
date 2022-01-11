#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Identifying and creating a mask for bad columns and pixels.

Version 1.0

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
    n: length of a list of 2D arrays
    imsize_y, imsize_x: shape (y dim, x dim) of first element in the list
    
    """
    
    n = len(filelist)
    
    try:
        frame0 = fits.getdata(filelist[0])
    except:
        print('Invalid filelist input!')
    
    imsize_y, imsize_x = frame0.shape
    
    return n, imsize_y, imsize_x

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
    
    n, imsize_y, imsize_x = file_prop(filelist)
    
    fits_stack = np.zeros((imsize_y, imsize_x , n))
    
    for ii in range(0, n):
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

def sigImage(filelist, n, imsize_y, imsize_x):
    """
    Create sigma image.
    
    Arguments
    -----------
    filelist: list of FITS files
    n: number of elements in filelist (int)
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
    sig = np.zeros((imsize_y, imsize_x , n))
    
    for ii in range(0, n):
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
    
    n, imsize_y, imsize_x = file_prop(filelist)
    
    ms = int(n * 0.1585)
    
    varim = np.zeros((imsize_y, imsize_x))
    wdtim = np.zeros((imsize_y, imsize_x))
    
    sig = sigImage(filelist, n, imsize_y, imsize_x)
    
    print('Creating variance image...')
    for j in range(imsize_y):
        for i in range(imsize_x):

            sig_sort = np.sort(sig[j,i,:])
            
            var1 = abs(sig_sort[int(n/2)] - sig_sort[ms])
            var2 = abs(sig_sort[int(n/2)] - sig_sort[n-ms])
            
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

def findClipidx(varim, mu, sigma, signum=5.):
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

def makeMask(varim, pxclipval, write=True):
    """ Returns (x,y) of identified bad/hot pixels
    
    Arguments
    -----------
    varim: variance image array (2D numpy array)
    
    Returns
    -----------
    x, y: array_like
        pixel positions
        can be used as inputs to makeRegions for ds9 .reg generation
        
    """
    
    imsize_y, imsize_x = varim.shape
    
    mask = np.zeros((imsize_y, imsize_x))
    
    y = np.where(varim > pxclipval)[0]
    x = np.where(varim > pxclipval)[1]
    
    mask[y,x] = 1.
    
    if write == True:
        outputname = 'mask.fits'
        fits.writeto(outputname, mask, overwrite=True)
        print('Written in ', outputname)
    
    return x, y