# LS4 CCD Testing Analysis Tools

An under-construction pipeline for DECam CCD testing for the LS4 Survey.

## Dependencies

* NumPy (1.20.0)
* SciPy.stats (1.6.2)
* Astropy (4.2.1)
* SEP (1.2.0)

## Measuring gain

There are two approaches used in the testing to obtain a gain measurement. 

### Photon transfer curve

A pair of flats at increasing exposure times can be used to produce a signal variance as a function of signal mean plot, the slope of which is the gain. This slope is in the shot noise regime. Two images are taken in order to subtract out the fixed pattern noise.

### X-ray

Using the x-ray approach involves manually placing a radioactive 55Fe source directly in front of the CCD. The idea is to extract the sources of the x-ray hits on the active surface and plot the flux distribution of these hits. The energies of the electron capture decays are well-known. The ratio of the peaks at 5.9 keV and 6.5 keV are used to find the gain.

## Measuring quantum efficiency

Images taken at a range of wavelengths (at fixed exposure time) can be used to generate a relative QE curve once dark-subtracted.
