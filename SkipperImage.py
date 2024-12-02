#!/usr/bin/env python
"""
Load and process skipper images.
"""

from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
from pylab import *
import pylab as plt

#This version does not contain the Tracer Pixel variable

def getKeyword(header,key,dval):
    try:
        return header[key] #number of samples taken
    except:
        print('Keyword "'+key+'" not found, using default value '+str(dval))
    
    return dval

class SkipperImage:

    def __init__(self,filename=None,avgMargin=4,skipSamples=0,skipRows=0,deTrend=True,nPre=7,nCol=362,nSamp=1):
        self.filename=filename
        self.avgMargin=avgMargin
        self.skipSamples=skipSamples
        self.skipRows=skipRows
        self.nPre=nPre
        self.nCol=nCol
        self.nSamp=nSamp
        self.deTrend=deTrend
        self.imageProcessed=False

    @property
    def header(self,amp=1):
        return self.fits[amp+1].header
        
    def processImage(self,filename=None):
        if((not self.imageProcessed) or ((self.filename != filename) and (filename != None))):
            self.loadImage(filename)
            self.correctSamples()
            self.averageSamples()
        self.imageProcessed=True

    def loadImage(self,filename=None):
        """Load FITS image into memory and grab header values.

        Parameters
        ----------
        filename : str, None
            FITS filename to load. If `None`, loads self.filename

        Returns
        -------
        None
        """
        self.imageProcessed=False
        if(filename != None):
            self.filename=filename

        if(self.filename == None):
            raise ValueError('SkipperImage::loadImage -  No filename given')

        try:
            self.fits=fits.open(self.filename)
        except:
            raise IOError('Could not open FITS file '+str(filename))

        #Support single or extended format
        self.extended = (len(self.fits) > 1)
        self.hoffset=0
        if(self.extended):
            self.nSamp=int(float(getKeyword(self.fits[0].header,'NSAMP',self.nSamp)))
            self.nCol=int(int(getKeyword(self.fits[1].header,'CCDNCOL',self.nCol))/2)
            self.nPre=int(getKeyword(self.fits[1].header,'CCDNPRES',self.nPre))
            if(self.fits[0].header['NAXIS'] < 2):
                self.hoffset=1
        else:
            self.nSamp=int(float(getKeyword(self.fits[0].header,'NSAMP',self.nSamp)))
            self.nCol=int(int(getKeyword(self.fits[0].header,'CCDNCOL',self.nCol))/2)
            self.nPre=int(getKeyword(self.fits[0].header,'CCDNPRES',self.nPre))

        self.osStart=self.nCol+self.nPre
        self.ampNum=len(self.fits)-self.hoffset

        self.data=list()
        for i in range(0,self.ampNum):
            #reshape the image axis
            self.r,cs=self.fits[i+self.hoffset].data.shape
            self.c=int(cs/self.nSamp)
            data3d=self.fits[i+self.hoffset].data.reshape(self.r,self.c,self.nSamp)
            self.data.append(data3d)
    
    def correctSamples(self):
        for i in range(0,self.ampNum):
            data3d=self.data[i][self.skipRows:,:,:]
            if(self.deTrend):
                ##generate identity array for each row
                indArr=np.zeros_like(data3d)+np.arange(0,self.c).reshape(1,self.c,1)

            #separate over-scan region
            dataOs=data3d[:,self.osStart+self.avgMargin:-self.avgMargin,:] #get overscan minus average margin
    
            #compute correction from overscan
            #osMed=np.median(dataOs,axis=1,keepdims=True) #average overscan for each row and sample
            osMed=sigma_clip(dataOs,axis=1).mean(axis=1,keepdims=True).data #average overscan for each row and sample, rejecting outliers
            if(self.deTrend):
                shifted=np.roll(osMed,1,axis=0) #generate array which is shifted[i]=osMean[i-1]
                slope=(osMed-shifted)/(self.c-1) #compute slope between current and previous overscan   
                corrArr=indArr*slope+shifted #compute correction for each pixel and sample as slope*c+os[r-1]
                corrArr[0,:,:]=osMed[0,0,:] #first row just gets normal overscan correction
                dataCorr=data3d-corrArr #compute corrected image
            else:
                dataCorr=data3d-osMed
            self.data[i]=dataCorr
    
    def averageSamples(self):
        """Calculate the average of the samples."""
        self.imgs=list()
        for i in range(0,self.ampNum):
            #average along the sampling axis, ignoring the first SKIP samples
            self.imgs.append(np.mean(self.data[i][:,:,self.skipSamples:],axis=2))
        self.imageProcessed=True

    def getImage(self,amp=1,filename=None,full=False):
        return self.getImages(filename,full=full)[amp-1]
                
    def getImages(self,filename=None,full=False):
        self.processImage(filename)

        if(full):
            return self.imgs
        else:
            trimmed=list()
            for img in self.imgs:
                trimmed.append(img[:,self.nPre:self.osStart])
            return trimmed

    def getOverscan(self,amp=1,filename=None):
        return self.getOverscans(filename)[amp-1]

    def getOverscans(self,filename=None):
        self.processImage(filename)

        trimmed=list()
        for img in self.imgs:
            trimmed.append(img[:,self.osStart:])
        return trimmed
                
    def plotSamples(self,row,column,amp=None):
        self.processImage()

        if(amp == None):
            for i in range(0,self.ampNum):
                plot(self.data[i][row,column,:],label='Amp '+str(i))
        else:
            plot(self.data[amp][row,column,:],label=str(row)+', '+str(column))


    def getResolution(self,amp=None,valCut=100):
        """Get the resolution (standard deviation), as a function of
        number of samples.

        Parameters
        ----------
        amp : int
            amplifier number
        valCut : fload
            cut outliers greater than this value

        Returns
        -------
        samps, sigmas : 
            number of samples and standard deviations
        """
        self.processImage()

        samps=list()
        for i in range(1,int(sqrt(self.nSamp))+1):
            samps.append(i**2)

        if(amp != None):
            stds=list()
            for samp in samps:
                vals=np.mean(self.data[amp-1][:,self.osStart+self.avgMargin:-self.avgMargin,self.skipSamples:self.skipSamples+samp],axis=2).flatten()
                stds.append(np.std(vals[abs(vals) < valCut]))
            return samps,stds
        else:
            sigmas=list()
            for i in range(0,self.ampNum):
                stds=list()
                for samp in samps:
                    vals=np.mean(self.data[i][:,self.osStart+self.avgMargin:-self.avgMargin,self.skipSamples:self.skipSamples+samp],axis=2).flatten()
                    stds.append(np.std(vals[abs(vals) < valCut]))
                sigmas.append(stds)
            return samps,sigmas


    def plotImage(self,amp=1,vmin=None,vmax=None,full=False):
        self.processImage()

        if(full):
            plt.imshow(self.imgs[amp-1],vmin=vmin,vmax=vmax,origin='lower',aspect='auto')
        else:
            plt.imshow(self.imgs[amp-1][:,self.nPre:self.osStart],vmin=vmin,vmax=vmax,origin='lower',aspect='auto')
        plt.colorbar(shrink=0.8)
        plt.title('Amp '+str(amp))

    def plotOverscan(self,amp=1,vmin=None,vmax=None,full=False):
        self.processImage()
        plt.imshow(self.imgs[amp-1][:,self.osStart:],vmin=vmin,vmax=vmax,origin='lower',aspect='auto')
        plt.colorbar(shrink=0.8)
        plt.title('Amp '+str(amp))

    def plotImages(self,plotLabel=None,vmin=None,vmax=None,full=False):
        plt.figure(figsize=(14,7))
        for i in range(1,self.ampNum+1):
            plt.subplot(2,2,i)
            self.plotImage(i,vmin=vmin,vmax=vmax,full=full)
        plt.suptitle(plotLabel)
        plt.tight_layout()

    def histImage(self,amp=1,full=False,skipCols=0,skipRows=0,scale=1.0,**kwargs):
        """
        Parameters
        ----------
        ...
        kwargs : passed to plt.hist
        
        Returns
        -------
        None
        """
        self.processImage()
        if(full):
            plt.hist(self.imgs[amp-1][skipRows:,skipCols:].flatten()*scale,histtype='step',**kwargs)
        else:
            plt.hist(self.imgs[amp-1][skipRows:,self.nPre+skipCols:self.osStart].flatten()*scale,histtype='step',**kwargs)

    def histOverscan(self,amp=1,full=False,skipCols=0,skipRows=0,scale=1.0,**kwargs):
        self.processImage()
        plt.hist(self.imgs[amp-1][skipRows:,self.osStart+skipCols:].flatten()*scale,histtype='step',**kwargs)

    def histImages(self,full=False,skipCols=0,skipRows=0,scale=1.0,**kwargs):
        for i in xrange(0,self.ampNum):
            self.histImage(i+1,full=full,skipCols=skipCols,skipRows=skipRows,scale=scale,label='Amp '+str(i),**kwargs)

    def save(self,filename=None):
        self.processImage()

        if(filename == None):
            fname=self.filename[0:-5]+'_processed.fits'
        else:
            fname=filename

        print('Saving to '+fname)

        if(self.hoffset > 0):
            prim=fits.PrimaryHDU()
            prim.header=self.fits[0].header
            newFITS=fits.HDUList([prim])
        else:
            newFITS=fits.HDUList()
        for i in range(0,self.ampNum):
            newImg = fits.ImageHDU(self.imgs[i])
            existKeys=newImg.header.keys()
            for card in self.fits[i+self.hoffset].header.cards:
                k=card[0]
                if(k in existKeys):
                    continue
                newImg.header.append(card)
            newImg.header.append(('SKIPAVG',True,'Skipper averaging has been done'))
            newFITS.append(newImg)
        newFITS.writeto(fname)
        
