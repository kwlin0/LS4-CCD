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

def combine_files(single_exposure_filenames, compress=True, header=None):
    primaryHDU = fits.PrimaryHDU(header=header)
    hdulist = [primaryHDU]
    for filenum in range(len(single_exposure_filenames)):
        data       = fits.open(single_exposure_filenames[filenum])[0].data
        header     = fits.open(single_exposure_filenames[filenum])[0].header
        if compress == False:
            imageHDU   = fits.ImageHDU(data, header=header)
        else:
            imageHDU   = fits.CompImageHDU(data, header=header, compression_type='RICE_1')
        hdulist.append(imageHDU)
    allHDU = fits.HDUList(hdulist)
    if compress == False:
        outname = os.path.dirname(single_exposure_filenames[0]) + ".fits"
        allHDU.writeto(outname, overwrite=True)
        print('Wrote', outname)
    else:
        outname = os.path.dirname(single_exposure_filenames[0]) + ".fits.fz"
        allHDU.writeto(outname, overwrite=True)
        print('Wrote', outname)
    
    
def main():
    parser = argparse.ArgumentParser(description="Combine a single exposure focal plane image FITS files into a single fpacked file. Output files in the parent directory of raw image files.")
    parser.add_argument("path", nargs='+', help="Raw file input: (1) to compress single exposure, use path/to/exp_00000/*.fits, (2) to compress a number of exposures, use path/to/exp*")
    args = parser.parse_args()
    
    inputdat = sorted(args.path)
    
    if inputdat[0].lower().endswith(('.fits', '.fit', '.FITS')):
        print('FITS files given\ncompressing as single exposure')
        combine_files(inputdat)
    
    else:
        print('Compressing all exposures...')
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
