#!/usr/bin/env python3

import argparse
from astropy.io import fits
from astropy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
import os

#natroot = 'https://astroarchive.noirlab.edu'
#adsurl = f'{natroot}/api/adv_search'

def get_subtracted(FITS_file):
    hdulist = fits.open(FITS_file)
    overscan_subtracted_hdulist = fits.HDUList(fits.PrimaryHDU(header=hdulist[0].header)) # gain-UNcorrected FITS file
    gcoverscan_subtracted_hdulist = fits.HDUList(fits.PrimaryHDU(header=hdulist[0].header)) # gain-corrected FITS file
    left_arr = np.empty((62, 4095, 1023))
    right_arr = np.empty((62, 4095, 1023))
    
    gcleft_arr = np.empty((62, 4095, 1023))
    gcright_arr = np.empty((62, 4095, 1023))

    for hdun, hdu in tqdm(enumerate(hdulist[1:63]), colour='green'):
        # Header stuff
        gainA = hdu.header['GAINA']
        gainB = hdu.header['GAINB']
        N_or_S = hdu.header['DETPOS'][0]
    
        # Do the subtraction
        left_overscan = hdu.data[51:4146, 7:56]
        left_overscan_mean = np.median(left_overscan, axis=1) # row-by-row mean of left overscan
        left_subtracted = hdu.data[51:4146, 57:1080] - np.vstack(left_overscan_mean)
        left_arr[hdun] = left_subtracted
    
    
        right_overscan = hdu.data[51:4146, 2105:2154]
        right_overscan_mean = np.median(right_overscan, axis=1) # row-by-row mean of right overscan
        right_subtracted =hdu.data [51:4146, 1081:2104] - np.vstack(right_overscan_mean)
        right_arr[hdun] = right_subtracted
    
        overscan_subtracted = np.concatenate((left_arr[hdun], right_arr[hdun]), axis=1)    # stitch the two amps back together
        overscan_subtracted_hdulist.append(fits.CompImageHDU(overscan_subtracted, hdu.header))
    
        if N_or_S == 'N':
            gcleft_arr[hdun] = left_arr[hdun] * gainA
            gcright_arr[hdun] = right_arr[hdun] * gainB
        elif N_or_S == 'S':
            gcleft_arr[hdun] = left_arr[hdun] * gainB
            gcright_arr[hdun] = right_arr[hdun] * gainA
    
        gcoverscan_subtracted = np.concatenate((gcleft_arr[hdun], gcright_arr[hdun]), axis=1)    # stitch the two amps back together
        gcoverscan_subtracted_hdulist.append(fits.CompImageHDU(gcoverscan_subtracted, hdu.header))
    
    outname = 'osub_' + FITS_file.split('/')[-1]
    overscan_subtracted_hdulist.writeto(outname, overwrite = True)

    gcoutname = 'gc_osub_' + FITS_file.split('/')[-1]
    gcoverscan_subtracted_hdulist.writeto(gcoutname, overwrite = True)
    return outname, gcoutname

def get_masks(FITS_file):
    """
    FITS_file should be a FITS image of a DARK DECam exposure
    """
    dataA_arr = np.empty((70 - 1, 4095, 1023))
    maskA_arr = np.empty((70 - 1, 4095, 1023))
    dataB_arr = np.empty((70 - 1, 4095, 1023))
    maskB_arr = np.empty((70 - 1, 4095, 1023))
    new_mask = fits.HDUList()
    new_gc = fits.HDUList() # gain-corrected FITS file
    hdulist = fits.open(FITS_file)
    for hdun, hdu in tqdm(enumerate(hdulist[1:63]), colour='green'):
        gainA = hdu.header['GAINA']
        gainB = hdu.header['GAINB']
        dataB_arr[hdun] = hdu.data[51:4146, 1081:2104] * gainB
        dataA_arr[hdun] = hdu.data[51:4146, 57:1080] * gainA
        maskA_arr[hdun] = stats.sigma_clip(dataA_arr[hdun], sigma_lower=10, sigma_upper=10).mask
        maskB_arr[hdun] = stats.sigma_clip(dataB_arr[hdun], sigma_lower=10, sigma_upper=10).mask
        mask_arr = np.concatenate((maskA_arr[hdun], maskB_arr[hdun]), axis=1)    # stitch the two amps back together
        new_mask.append(fits.CompImageHDU(mask_arr, hdu.header))   # make image data compressed
        gc_arr = np.concatenate((dataA_arr[hdun], dataB_arr[hdun]), axis=1)    # stitch the two amps back together
        new_gc.append(fits.CompImageHDU(gc_arr, hdu.header))   # make image data compressed
    outname = 'crm_' + FITS_file.split('/')[-1]
    outname_gc = 'gc_' + FITS_file.split('/')[-1]
    new_mask.writeto(outname, overwrite = True)
    new_gc.writeto(outname_gc, overwrite = True)
    return outname

def get_trimmed(FITS_file):
    """
    Trim the stacked images to eliminate edge artifacts from processing
    
    """
    hdulist = fits.open(FITS_file)
    trimeddata = hdulist[1].data[220:2000,185:4200]
    outname = 'trimmed/trim_'+FITS_file.split('/')[-1]
    fits.writeto(outname, trimeddata, hdulist[1].header, overwrite=True)
    print('Wrote ' + outname)

def get_chucks(data, npix=500):
    """
    Get chucks from input data image array with NO overlaps between chucks
    i.e., input array data will be divided into as many npix x npix sized chucks as possible
    without overlap between any two chucks
    Parameters
    -----------
    data = input image array
    npix = number of pixels to divide into in a square chuck
    
    Returns
    --------
    data_chunked : ndarray of size nchucks of npix x npix size chucks
    """
    # Integer number of npix sized chunks that image can be divided into
    nchunks_y = data.shape[0] // npix
    nchunks_x = data.shape[1] // npix
    
    # Remainder when divided into npix sized chunks
    remainder_y = data.shape[0] % npix
    remainder_x = data.shape[1] % npix
    
    # Split data into nchunks_y chunks, change back to numpy array of shape (nchunks_y, npix, data.shape[1])
    data_vsplit = np.array(np.split(data[:-remainder_y,:], nchunks_y, axis=0))
    
    # Now split along the data.shape[1] axis, which corresponds to axis=2 of data_vsplit
    data_chunked = np.array(np.split(data_vsplit[:,:,:-remainder_x], nchunks_x, axis=2))
    
    return data_chunked

def get_rand_chuck(data, npix=500):
    """
    Get a single chuck from an image array
    
    Parameters
    -----------
    data = input image array
    npix = number of pixels to divide into in a square chuck
    
    Returns
    --------
    data_chunked : ndarray of size nchucks of npix x npix size chucks
    """
    xrand = np.random.randint(0, data.shape[1]-npix)
    yrand = np.random.randint(0, data.shape[0]-npix)
    data_chunked = data[yrand:yrand+npix, xrand:xrand+npix]
    return data_chunked, xrand, yrand

def get_rand_chucks(data, npix=500, nchucks=10):
    """
    Get nchuck number of subarrays from an image array by random selection
    
    Parameters
    -----------
    data = input image array
    npix = number of pixels to divide into in a square chuck
    
    Returns
    --------
    data_chunked : ndarray of size nchucks of npix x npix size chucks
    xrand : 1D array of randomly drawn x coordinate values from input data image array
    yrand : 1D array of randomly drawn y coordinate values from input data image array
    """
    data_chunked = np.zeros((nchucks, npix, npix))
    xrand = np.zeros((nchucks))
    yrand = np.zeros((nchucks))
    for i in range(nchucks):
        square, x, y = get_rand_chuck(data, npix)
        data_chunked[i] = square
        xrand[i] = x
        yrand[i] = y
    return data_chunked, xrand, yrand

def gaussian_noise(image, amount):
    """
    Generate simulated read noise.
    
    Parameters
    ----------
    
    image: numpy array
        Image whose shape the noise array should match.
    amount : float
        Amount of read noise, in electrons.
    gain : float, optional
        Gain of the camera, in units of electrons/ADU.
    """
    # Set up the random number generator, allowing a seed to be set from the environment
    seed = os.getenv('GUIDE_RANDOM_SEED', None)

    if seed is not None:
        seed = int(seed)
    
    # This is the generator to use for any image component which changes in each image, e.g. read noise
    # or Poisson error
    noise_rng = np.random.default_rng(seed)
    noise = noise_rng.normal(scale=amount, size=image.shape)
    
    return noise


def sky_background(image, sky_counts):
    """
    Generate sky background, optionally including a gradient across the image (because
    some times Moons happen).
    
    Parameters
    ----------
    
    image : numpy array
        Image whose shape the cosmic array should match.
    sky_counts : float
        The target value for the number of counts (as opposed to electrons or 
        photons) from the sky.
    gain : float, optional
        Gain of the camera, in units of electrons/ADU.
    """
    # Set up the random number generator, allowing a seed to be set from the environment
    seed = os.getenv('GUIDE_RANDOM_SEED', None)

    if seed is not None:
        seed = int(seed)
    
    # This is the generator to use for any image component which changes in each image, e.g. read noise
    # or Poisson error
    noise_rng = np.random.default_rng(seed)
    sky_im = noise_rng.poisson(sky_counts, size=image.shape)
    
    return sky_im

def add_noise(data, level, sigma=0):
    synthetic_image = np.zeros([data.shape[0], data.shape[1]])
    new_image = data + sky_background(synthetic_image, level) + gaussian_noise(synthetic_image, sigma)
    return new_image

def main():
    parser = argparse.ArgumentParser(description="Generate CR mask for a raw DECam image.")
    parser.add_argument("path", nargs='+', help="FITS file input")
    parser.add_argument('-m', '--mask', action='store_true', help='Apply mask on DECam dark FITS image')
    parser.add_argument('-s', '--subtraction', action='store_true', help='Apply overscan subtraction on DECam FITS image')
    args = parser.parse_args()

    inputdat = np.array(args.path)
    for file in inputdat:
        if args.mask:
            outname = get_masks(file)
            print('Wrote ' + str(outname))
        if args.subtraction:
            outname, gcoutname = get_subtracted(file)
            print('Wrote ' + str(outname) + ', ' + str(gcoutname))
            
        print("Done")

if __name__ == '__main__':
    main()