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
import sep


# Overscan subtraction routine


def overscanSubtraction(filename, extension, destination):
    """ Overscan subtraction routine
    
    Keyword arguments:
    filename -- the path and file of image
    extension -- FITS image extension
    destination -- directory of output image (current working directory/destination.fits)
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


def biasSubtraction(filelist, superbias):
    """ Subtracts super bias from filelist images
    
    filelist -- list of images to be super bias subtracted
    superbias -- super bias FITS file (note extension)
    """    
    # Super bias data
    sb_data = fits.getdata(superbias)
    
    for im in filelist:
        im_data = fits.getdata(im)
        bias_sub = im_data - sb_data 
        
        fits.writeto(os.getcwd()+'/bs_'+im.rsplit('/',1)[1], bias_sub, fits.getheader(im), overwrite=True)
        print('Written in:', os.getcwd()+'/bs_'+im.rsplit('/',1)[1])


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
    
    filename -- path to FITS file
    threshold -- pixel value to extract above
    
    Returns: object with different params
    
    """
    
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