import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.style
matplotlib.style.use('classic')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.colors import LogNorm
from glob import glob
import sys,os,getopt

def main(argv):

    try:
        opts, args = getopt.getopt(argv,"hi:")
    except getopt.GetoptError:
        print('light_curves.py -i <inputfile>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('light_curves.py -i <inputfile>')
            sys.exit()
        elif opt == "-i":
            inputdir = arg

    print('Plotting light curves, reading from', inputdir)

    if os.path.isfile(inputdir+'/incidentLightCurve.fits'):
        incidentLightCurve = fits.open(inputdir+'/incidentLightCurve.fits')
        time = incidentLightCurve[1].data.field(1)
        flux = incidentLightCurve[1].data.field(3)
        plt.figure(figsize=(15,5))
        line = plt.plot(time,flux,linewidth=2.0)
        plt.setp(line, color='r')
        plt.xlabel('MJD', fontsize=14)
        plt.ylabel('Flux (electrons per exposure)', fontsize=14)
        plt.grid(1)
        ax = plt.subplot(111)
        ax.ticklabel_format(useOffset=False)
        plt.xlim(time[0],time[len(time)-1])    
        if not os.path.exists(os.getcwd()+'/png'):
            os.makedirs(os.getcwd()+'/png')
        plt.savefig('png/incidentLightCurve.png')

    if os.path.isfile(inputdir+'/extractedLightCurve.fits'):
        extractedLightCurve = fits.open(inputdir+'/extractedLightCurve.fits')
        time = extractedLightCurve[1].data.field(1)
        flux = extractedLightCurve[1].data.field(3)
        plt.figure(figsize=(15,5))
        plt.plot(time,flux,'o',ms=5)
        plt.xlabel('MJD', fontsize=14)
        plt.ylabel('Flux (ADU)', fontsize=14)
        ax = plt.subplot(1,1,1)
        ax.ticklabel_format(useOffset=False)
        plt.grid(1)
        plt.xlim(time[0],time[len(time)-1])    
        if not os.path.exists(os.getcwd()+'/png'):
            os.makedirs(os.getcwd()+'/png')
        plt.savefig('png/extractedLightCurve.png')

    if os.path.isfile(inputdir+'/idealLightCurve.fits'):
        idealLightCurve = fits.open(inputdir+'/idealLightCurve.fits')
        time = idealLightCurve[1].data.field(1)
        flux = idealLightCurve[1].data.field(3)
        plt.figure(figsize=(15,5))
        plt.plot(time,flux,'o',ms=5)
        plt.xlabel('MJD', fontsize=14)
        plt.ylabel('Flux (ADU)', fontsize=14)
        ax = plt.subplot(1,1,1)
        ax.ticklabel_format(useOffset=False)
        plt.grid(1)
        plt.xlim(time[0],time[len(time)-1])    
        if not os.path.exists(os.getcwd()+'/png'):
            os.makedirs(os.getcwd()+'/png')
        plt.savefig('png/idealLightCurve.png')
        
    if os.path.isfile(inputdir+'/incidentLightCurveWithPhotonNoise.fits'):
        if os.path.isfile(inputdir+'/incidentLightCurve.fits'):
            incidentLightCurveWithPhotonNoise = fits.open(inputdir+'/incidentLightCurveWithPhotonNoise.fits')
            incidentLightCurve = fits.open(inputdir+'/incidentLightCurve.fits')
            time = incidentLightCurveWithPhotonNoise[1].data.field(1)
            flux = incidentLightCurveWithPhotonNoise[1].data.field(3)
            inputTime = incidentLightCurve[1].data.field(1)
            inputFlux = incidentLightCurve[1].data.field(3)
            inputFlux *= inputFlux.size/flux.size
            plt.figure(figsize=(15,5))
            plt.plot(time,flux,'o',ms=5,label='Incident light curve with photon noise')
            plt.plot(inputTime,inputFlux,label='Incident light curve',linewidth=3.0,color='red')
            plt.xlabel('MJD', fontsize=14)
            plt.ylabel('Flux (electrons per stacked image)', fontsize=14)
            legend = plt.legend(loc='lower right',shadow=True,fontsize=14)
            ax = plt.subplot(1,1,1)
            ax.ticklabel_format(useOffset=False)
            plt.grid(1)
            plt.xlim(time[0],time[len(time)-1])    
            if not os.path.exists(os.getcwd()+'/png'):
                os.makedirs(os.getcwd()+'/png')
            plt.savefig('png/incidentLightCurveWithPhotonNoise.png')
        
    if len(glob(inputdir+'/*SCI_COR_Lightcurve-DEFAULT*.fits'))>0 and os.path.isfile(glob(inputdir+'/*SCI_COR_Lightcurve-DEFAULT*.fits')[0]):
        data_reduction_lc = fits.open(glob(inputdir+'/*SCI_COR_Lightcurve-DEFAULT*.fits')[0])
        time = data_reduction_lc[1].data.field(1)
        flux = data_reduction_lc[1].data.field(3)
        plt.figure(figsize=(15,5))
        plt.plot(time,flux,'o',ms=5,color='green')
        plt.xlabel('MJD', fontsize=14)
        plt.ylabel('Flux (ADU per stacked image)', fontsize=14)
        ax = plt.subplot(1,1,1)
        ax.ticklabel_format(useOffset=False)
        plt.grid(1)
        plt.xlim(time[0],time[len(time)-1])    
        if not os.path.exists(os.getcwd()+'/png'):
            os.makedirs(os.getcwd()+'/png')
        plt.savefig('png/dataReductionLightCurve.png')

    if os.path.isfile(inputdir+'/extractedNormalizedLightCurve.fits'):
        if os.path.isfile(inputdir+'/incidentNormalizedLightCurve.fits'):
            extractedNormalizedLightCurve = fits.open(inputdir+'/extractedNormalizedLightCurve.fits')
            incidentNormalizedLightCurve = fits.open(inputdir+'/incidentNormalizedLightCurve.fits')
            time = extractedNormalizedLightCurve[1].data.field(1)
            flux = extractedNormalizedLightCurve[1].data.field(3)
            inputTime = incidentNormalizedLightCurve[1].data.field(1)
            inputFlux = incidentNormalizedLightCurve[1].data.field(3)
            plt.figure(figsize=(15,5))
            plt.plot(time,flux,'o',ms=5,label='Extracted flux')
            plt.plot(inputTime,inputFlux,label='incident light curve',linewidth=3.0,color='red')
            plt.xlabel('MJD', fontsize=14)
            plt.ylabel('Flux normalized to mean outside transit', fontsize=14)
            plt.grid(1)
            legend = plt.legend(loc='lower right',shadow=True,fontsize=14)
            ax = plt.subplot(111)
            ax.ticklabel_format(useOffset=False)
            plt.xlim(time[0],time[len(time)-1])    
            if not os.path.exists(os.getcwd()+'/png'):
                os.makedirs(os.getcwd()+'/png')
            plt.savefig('png/extractedNormalizedLightCurve.png')

    if os.path.isfile(inputdir+'/idealNormalizedLightCurve.fits'):
        if os.path.isfile(inputdir+'/incidentNormalizedLightCurve.fits'):
            idealNormalizedLightCurve = fits.open(inputdir+'/idealNormalizedLightCurve.fits')
            incidentNormalizedLightCurve = fits.open(inputdir+'/incidentNormalizedLightCurve.fits')
            time = idealNormalizedLightCurve[1].data.field(1)
            flux = idealNormalizedLightCurve[1].data.field(3)
            inputTime = incidentNormalizedLightCurve[1].data.field(1)
            inputFlux = incidentNormalizedLightCurve[1].data.field(3)
            plt.figure(figsize=(15,5))
            plt.plot(time,flux,'o',ms=5,label='Extracted flux without detector effects')
            plt.plot(inputTime,inputFlux,label='incident light curve',linewidth=3.0,color='red')
            plt.xlabel('MJD', fontsize=14)
            plt.ylabel('Flux normalized to mean outside transit', fontsize=14)
            plt.grid(1)
            legend = plt.legend(loc='lower right',shadow=True,fontsize=14)
            ax = plt.subplot(111)
            ax.ticklabel_format(useOffset=False)
            plt.xlim(time[0],time[len(time)-1])    
            if not os.path.exists(os.getcwd()+'/png'):
                os.makedirs(os.getcwd()+'/png')
            plt.savefig('png/idealNormalizedLightCurve.png')

    if os.path.isfile(inputdir+'/extractedOMinusCLightCurve.fits'):
        extractedOMinusCLightCurve = fits.open(inputdir+'/extractedOMinusCLightCurve.fits')
        time = extractedOMinusCLightCurve[1].data.field(1)
        flux = extractedOMinusCLightCurve[1].data.field(3)
        plt.figure(figsize=(15,5))
        plt.plot(time,flux,'o',ms=5)
        plt.xlabel('MJD', fontsize=14)
        plt.ylabel('Observed - Calculated Flux', fontsize=14)
        plt.grid(1)
        ax = plt.subplot(111)
        ax.ticklabel_format(useOffset=False)
        plt.xlim(time[0],time[len(time)-1])    
        if not os.path.exists(os.getcwd()+'/png'):
            os.makedirs(os.getcwd()+'/png')
        plt.savefig('png/extractedOMinusCLightCurve.png')

    if os.path.isfile(inputdir+'/noiseCurve.fits'):
        noiseCurve = fits.open(inputdir+'/noiseCurve.fits')
        time = noiseCurve[1].data.field(0)
        noise = noiseCurve[1].data.field(1)
        plt.figure(figsize=(15,5))
        plt.plot(time,noise,linewidth=2.0)
        plt.xlabel('Integration time (hours)', fontsize=14)
        plt.ylabel('Noise standard deviation (ppm)', fontsize=14)
        plt.grid(1)
        xmin, xmax, ymin, ymax = plt.axis()
        plt.axis([0, 3, 0, ymax*1.05])
        ax = plt.subplot(111)
        ax.ticklabel_format(useOffset=False)
        plt.xlim(time[0],time[len(time)-1])    
        if not os.path.exists(os.getcwd()+'/png'):
            os.makedirs(os.getcwd()+'/png')
        plt.savefig('png/noiseCurve.png')
        
    if os.path.isfile(inputdir+'/idealNoiseCurve.fits'):
        noiseCurve = fits.open(inputdir+'/noiseCurve.fits')
        time = noiseCurve[1].data.field(0)
        noise = noiseCurve[1].data.field(1)
        idealNoiseCurve = fits.open(inputdir+'/idealNoiseCurve.fits')
        idealNoise = idealNoiseCurve[1].data.field(1)
        plt.figure(figsize=(15,5))
        plt.plot(time,noise,linewidth=2.0,label='With detector effects')
        plt.plot(time,idealNoise,linewidth=2.0,label='Without detector effects')
        plt.xlabel('Integration time (hours)', fontsize=14)
        plt.ylabel('Noise standard deviation (ppm)', fontsize=14)
        plt.grid(1)
        legend = plt.legend(loc='upper right',shadow=True,fontsize=14)
        xmin, xmax, ymin, ymax = plt.axis()
        plt.axis([0, 3, 0, ymax*1.05])
        ax = plt.subplot(111)
        ax.ticklabel_format(useOffset=False)
        plt.xlim(time[0],time[len(time)-1])    
        if not os.path.exists(os.getcwd()+'/png'):
            os.makedirs(os.getcwd()+'/png')
        plt.savefig('png/idealNoiseCurve.png')

if __name__ == "__main__":
   main(sys.argv[1:])
