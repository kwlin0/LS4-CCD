# LS4 CCD Testing Analysis Tools

An under-construction pipeline for DECam CCD testing for the LS4 Survey. The controller used for CCD testing was the 4-channel LTA developed by Fermilab. The controller to be used at the observatory will be the STA Archon controller.

## Requirements

* Numpy (>=1.20.0)
* Scipy (>=1.6.2)
* Astropy (>=4.2.1)
* sep (>=1.2.0)
* Optional: pylatex (>=1.4.1)

## CCD layout

The CCDs used here are identical to the DECam CCDs and are excess devices that remained after the completion of the DECam instrument. There are two readout amplifiers for each CCD and the ideal and typical readout mode is using the two amplifiers to clock out the charge in a split readout mode. The figure below demonstrates the silicon layout of these CCDs. The dark gray areas indicate "dead silicon" which are not photosensitive. The outer edges of this dead silicon demarcates the aluminum outline around the CCD die on the wafer. There is 100 µm between the aluminum lines of two die on the wafer where the saw cuts to separate the die. Typically, there would be an additional ~50 µm distance of silicon from this aluminum line to where the physical edge of the die depending on exactly where the dicing saw has cut.

![DECam_CCD_silicon](https://github.com/user-attachments/assets/6c940688-e0e6-4b76-abdd-af285ed2667b)

During CCD testing, the best performing CCDs were ranked and then installed in the central four positions of the focal plane. The criteria for best performance was primarily defined to be detectors that were measured to have the lowest readout noise levels and cosmetic issues with no evidence of inherent CCD problems such as charge transfer efficiency and linear response.

## Controller

Sensors were tested and characterized with the 4-channel LTA. With testing completed, LS4 will now exclusively use the Archon controller. The Archon controller communicates via a raw port and can send data at the full 1 GB network rate. It can also be configured and used with a GUI interface. We use four Archon controllers in a leader-follower configuration, such that all the clocking is synchronized. Each controller reads the video outputs of eight CCDs. The software for running the controllers synchronously is located in a separate repository (https://github.com/dlrabinowitz/ls4_control) and is built on top of the sdss-archon software.

### Measuring gain

There are two approaches used in the testing to obtain a gain measurement. The primary method was using the photon transfer transfer curve (PTCs) which relies on a reasonably flat illumination of the detector. The secondary method using X-rays was only used for a select number of detectors tested.

#### Photon transfer curve

A pair of flats at increasing exposure times can be used to produce a signal variance as a function of signal mean plot, the slope of which is the gain. This slope is in the shot noise regime. Two images are taken in order to subtract out the fixed pattern noise.

#### X-ray

Using the x-ray approach involves manually placing a radioactive 55Fe source directly in front of the CCD. The idea is to extract the sources of the x-ray hits on the active surface and plot the flux distribution of these hits. The energies of the electron capture decays are well-known. The ratio of the peaks at 5.9 keV and 6.5 keV are used to find the gain. A series of darks (we have been using 30) are taken at 60 second exposures. X-rays also may reveal evidence of charge transfer inefficiency, the telltale signs of which are pixel streaking on the x-ray hits.

### Measuring quantum efficiency

Images taken at a range of wavelengths (at fixed exposure time) can be used to generate a relative QE curve once dark-subtracted.

### Bad columns and pixels

Biases are used for identifying defects on the active surface. A large number of biases to work with is preferable.

### Charge transfer efficiency

Both X-rays and flat field data can be used to calculate the charge transfer efficiency (CTE) of the detector. For X-rays, the classic X-ray horizontal stacking plot can be used to determine if there is a slope in the X-ray event line between the leading and trailing edges of the CCD, which is indicative of charge transfer inefficiency. Using flat field data, the extended-pixel edge response technique can be used.

## File compression

The raw image output from the controllers for a given exposure is a single FITS image for each readout amplifier, yielding a total of 64 independent FITS image files for the entire focal plane of 32 CCDs in dual amplifier readout mode. Before files are transferred, the FITS files for each exposure should be combined and compressed into a single fpacked FITS file with 64 image HDU extensions. This can be done with the `fpack.py` script which uses the Rice compression algorithm by default.
