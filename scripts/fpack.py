#!/usr/bin/env python3

############################
# -*- coding: utf-8 -*-
#
# @Author: Kenneth Lin
# @Date: 2024-12-03
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

# sample PrimaryHDU header
hdu_0 = {'OBSMODE' : ['survey', 'Observing mode'],
         'FOCUS'   : [28.10000, 'Telescope focus position'],
         'INSTRUME': ['LS4Cam', 'La Silla Schmidt Southern Survey Camera'],
         'OBSERVAT': ['ESO La Silla', 'ESO La Silla Observatory'],
         'TELESCOP': ['ESO 1.0-m Schmidt'],
         'LATITUDE': [-29.25444, 'Telescope Latitude  (degrees north)'],
         'LONGITUD': [70.74167, 'Telescope Longitude (degrees west)'],
         'ELEVATIO': [2347, 'Telescope Elevation (meters)'],
         'CCDMODE' : ['dual', 'Dual or single amplifier readout mode'] }

def modhead(hdu, header_dict):
    """
    Modify FITS header given a FITS HDU and header data structured as a dictionary `header_dict`
    """
    for k,v in header_dict.items():
        hdu.header[k] = v[0]
        if len(v) == 2:
            hdu.header.comments[k] = v[1]

def combine_files(single_exposure_filenames, compress=True, primaryheader=None):
    primaryHDU = fits.PrimaryHDU()
    if primaryheader is not None:
        modhead(primaryHDU, primaryheader)
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
        size = os.stat(outname).st_size
        print("Wrote", outname, f"@ {round(size/(pow(1024,2)), 2)} MB")
    else:
        outname = os.path.dirname(single_exposure_filenames[0]) + ".fits.fz"
        allHDU.writeto(outname, overwrite=True)
        size = os.stat(outname).st_size
        print("Wrote", outname, f"@ {round(size/(pow(1024,2)), 2)} MB")
    
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
                combine_files(single_exposure_filenames, primaryheader=hdu_0)
            # If not a directory, skip
            else:
                print('Skipping', file)

if __name__ == '__main__':
    main()