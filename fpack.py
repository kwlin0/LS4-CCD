#!/usr/bin/env python3

############################
# -*- coding: utf-8 -*-
#
# @Author: Kenneth Lin
# @Date: 2024-12-02
# @Filename: fpack.py
#
# use astropy.io utility to
# fpack the FITS files for a
# single exposure into a single
# compressed fits.fz file
#
###############################

import os
import argparse
from astropy.io import fits
from glob import glob

def combine_files(single_exposure_filenames, compress=True):
    primaryHDU = fits.PrimaryHDU()
    hdulist = [primaryHDU]
    for filenum in range(len(single_exposure_filenames)):
        data       = fits.open(single_exposure_filenames[filenum])[0].data
        header     = fits.open(single_exposure_filenames[filenum])[0].header
        imageHDU   = fits.CompImageHDU(data, header=header)
        if compress == False:
            imageHDU   = fits.ImageHDU(data, header=header)
        else:
            imageHDU   = fits.CompImageHDU(data, header=header)
        hdulist.append(imageHDU)
    allHDU = fits.HDUList(hdulist)
    if compress == False:
        allHDU.writeto(single_exposure_filenames[0].split('/')[-2]+".fits", overwrite=True)
    else:
        allHDU.writeto(single_exposure_filenames[0].split('/')[-2]+".fits.fz", overwrite=True)
    
    
def main():
    parser = argparse.ArgumentParser(description="Combine a single exposure focal plane image FITS files into a single fpacked file. Output files in the directory where this script is run.")
    parser.add_argument("path", nargs='+', help="Raw file input: (1) to compress single exposure, use path/to/exp_00000/*.fits, (2) to compress a number of exposures, use path/to/exp*")
    args = parser.parse_args()
    
    inputdat = sorted(args.path)
    
    if inputdat[0].lower().endswith(('.fits', '.fit', '.FITS')):
        print('FITS files given - compressing as single exposure')
        combine_files(inputdat)
    
    else:
        print('Compressing all exposures')
        for file in inputdat:
            # Verify that given path is an exposure directory
            if os.path.isdir(file):
                print(file)
                single_exposure_filenames = sorted(glob(file+'/test*.fits'))
                combine_files(single_exposure_filenames)
            # If not a directory, skip
            else:
                print('Skipping', file)

if __name__ == '__main__':
    main()