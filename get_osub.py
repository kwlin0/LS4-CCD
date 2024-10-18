#!/usr/bin/env python3

import argparse
import numpy as np
from astropy.io import fits

def get_ccdinfo(filename):
    return fits.getheader(filename)['CCD_NAME'], fits.getheader(filename)['AMP_NAME'], fits.getheader(filename)['TAP_NAME']

def overscan_subtraction(filename):
    """
    Row-by-row overscan subtraction
    """
    for imfileN, imfile in enumerate(filename):
        # Header stuff
        ccd_name, amp_name, tap_name = get_ccdinfo(imfile)
        print(imfile, amp_name)
        data = fits.getdata(imfile)
        flippedx_data = np.flip(data, axis=1) # copy a version of the data as a flipped image along x axis for the AMP_NAME = RIGHT so we can capture the overscan/image regions correctly
    
        if amp_name == 'LEFT':
            # Do the subtraction
            overscan_LEFT = data[:, 1031:]
            overscan_mean_LEFT = np.median(overscan_LEFT, axis=1) # row-by-row mean of left overscan
            subtracted_LEFT = data - np.vstack(overscan_mean_LEFT) # overscan subtract the whole image regardless of overscans/prescans
            osub_reconstructed_LEFT = subtracted_LEFT
            outname = 'osub_'+imfile.split('/')[-1]
            fits.writeto(outname, osub_reconstructed_LEFT, header=fits.getheader(imfile), overwrite=True)
        elif amp_name == 'RIGHT':
            overscan_RIGHT = flippedx_data[:, 1031:]
            overscan_mean_RIGHT = np.median(overscan_RIGHT, axis=1) # row-by-row mean of right overscan
            subtracted_RIGHT = flippedx_data - np.vstack(overscan_mean_RIGHT) # overscan subtract the whole image regardless of overscans/prescans
            subtracted_RIGHT_unflip = np.flip(subtracted_RIGHT, axis=1) # extra step here to unflip the image
            osub_reconstructed_RIGHT = subtracted_RIGHT_unflip
            outname = 'osub_'+imfile.split('/')[-1]
            fits.writeto(outname, osub_reconstructed_RIGHT, header=fits.getheader(imfile), overwrite=True)
            
def main():
    parser = argparse.ArgumentParser(description="Do overscan subtraction of images")
    parser.add_argument("path", nargs='+', help="FITS files input")
    args = parser.parse_args()
    files = args.path

    overscan_subtraction(files)
    print('Done!')
    
if __name__ == '__main__':
    main()