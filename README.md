# LS4-CCD

An under-construction pipeline for DECam CCD testing for the LS4 Survey.

## Packages used

Numpy, Scipy.stats, Astropy, SEP

## Measuring gain

For exposed (X-ray) images, only the pixel counts of the X-ray signal is used in computing the gain. For each exposure time, this is done by extracting all signal from a filelist of the same exposure time of FITS images, doing photometry on the fluxes of these extracted signals, and then taking the difference (and sum) of the different combination pairs of images in the filelist. Mean counts comes from the mean from the summed image pairs and the standard deviation comes from the difference pairs, where the histogram of counts after taking the difference is fitted to a Gaussian distribution to find the standard deviation and then the FWHM. The variance is taking to be the FWHM^2 / 2.
