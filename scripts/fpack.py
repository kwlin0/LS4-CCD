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

# Header entries
basic_info = {
                'INSTRUME': ['LS4Cam', 'La Silla Schmidt Southern Survey Camera'],
                'OBSERVAT': ['ESO La Silla', 'ESO La Silla Observatory'],
                'TELESCOP': ['ESO 1.0-m Schmidt'],
                'LATITUDE': [-29.25444, 'Telescope Latitude  (degrees north)'],
                'LONGITUD': [70.74167, 'Telescope Longitude (degrees west)'],
                'ELEVATIO': [2347, 'Telescope Elevation (meters)'],
                'CCDMODE' : ['dual', 'Dual or single amplifier readout mode'], 
                'DATASECL': ['[7:1030,1:4096]', 'Data section from amp L'],
                'BIASSECL' : ['[1031:1050,1:4096]', 'Overscan region L'],
                'PRESECL' : ['[1:6,1:4096]', 'Prescan region L'],
                'DATASECR' : ['[21:1044,1:4096]', 'Data section from amp R'],
                'BIASSECR' : ['[1:20,1:4096]', 'Overscan region R'],
                'PRESECR' : ['[1045:1050,1:4096]', 'Prescan region R']
             }

dynamic_info = {
                'PROJECT' : ['LEG', 'determines whether public or proprietary'],
                'TARGETID' : [0000, 'Target or field name'] 
               }

left_amp   = {
                'BIASSEC' : ['[1031:1050,1:4096]', 'Overscan region L']
             }

right_amp  = {
                'BIASSEC' : ['[1:20,1:4096]', 'Overscan region R']
             }

def modhead(hdu, header_dict, loc=None):
    """
    Modify FITS header given a FITS HDU and header data structured as a dictionary `header_dict`
    By default, the new header entries are appended at the end of the header file. Optionally, use `loc`
    argument to specify where the insert the header entry. Note FITS headers are zero-indexed from the top.
    """
    nloc = 0
    for k,v in header_dict.items():
        # First check if keyword already exists. If it does, update keyword in place.
        if k in hdu.header:
            hdu.header[k] = v[0]
            if len(v) == 2:
                hdu.header.comments[k] = v[1]
        # For new keywords, give option loc of where to put them in the header.
        else:
            if loc is not None:
                if len(v) == 1:
                    hdu.header.insert(loc+nloc, (k, v[0]))
                elif len(v) == 2:
                    hdu.header.insert(loc+nloc, (k, v[0], v[1]))
                nloc += 1
            # Otherwise append the new keywords at the end of the header.
            else:
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
        hdu_header = fits.open(single_exposure_filenames[filenum])[0].header 
        if compress == False:
            imageHDU   = fits.ImageHDU(data, header=hdu_header)
        else:
            imageHDU   = fits.CompImageHDU(data, header=hdu_header, compression_type='RICE_1')
        
        ## Header additions and updating ##
        modhead(imageHDU, basic_info, loc=8)
        modhead(imageHDU, dynamic_info, loc=21)

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
    """
    Script can be used in one of two ways:
    1. To compress one exposure, the inputs should be the 64 FITS files from the controller.
    2. To compress many exposures in one run, give the exposure directories with (or without) a wildcard.
    """
    parser = argparse.ArgumentParser(description="Combine a single exposure focal plane image FITS files into a single fpacked file. Output files in the parent directory of raw image files.")
    parser.add_argument("path", nargs='+', help="Raw file input: (1) to compress single exposure, use path/to/exp_00000/*.fits, (2) to compress a number of exposures, use path/to/exp*")
    args = parser.parse_args()
    
    inputdat = sorted(args.path)
    
    if inputdat[0].lower().endswith(('.fits', '.fit', '.FITS')):
        print('FITS files given\ncompressing as single exposure')
        combine_files(inputdat, primaryheader=basic_info)
    
    else:
        print('Compressing all exposures...')
        for file in inputdat:
            # Verify that given path is an exposure directory
            if os.path.isdir(file):
                print(file)
                single_exposure_filenames = sorted(glob(file+'/test*.fits'))
                combine_files(single_exposure_filenames, primaryheader=basic_info)
            # If not a directory, skip
            else:
                print('Skipping', file)

if __name__ == '__main__':
    main()
