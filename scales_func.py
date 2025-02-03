import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
#from photutils.background import Background2D, MedianBackground
from astropy.convolution import Gaussian2DKernel
import astropy.convolution.convolve as convolve
from astropy.visualization import (MinMaxInterval, SqrtStretch,LogStretch,ZScaleInterval,HistEqStretch,LinearStretch,
                                   PowerStretch,ImageNormalize,simple_norm)
from astropy.io import fits
from matplotlib.patches import Rectangle
from math import log, log10, ceil, floor, exp, sqrt
from scipy import integrate,interpolate
import numpy.random as random
import os
from scales_etc_lib import *

def scales_etc_snr(filter,itime,nframes,mag):
    # constants
    c_km = 2.9979E5      # km/s
    c = 2.9979E10       # cm/s
    h = 6.626068E-27    # cm^2*g/s
    k = 1.3806503E-16   # cm^2*g/(s^2*K)
    Ang = 1E-8          # cm
    mu = 1E-4           # cm
    darkcurrent = 0.002 #e- / s
    readnoise = 5 # e-
    area=76*10000 # cm2
    scale=0.02


    filters = {
    "K" : (1.95,2.45,2.185,200),
    "L" : (2.90,4.15,3.46,80),
    "M" : (4.5,5.2,5.83,140),
    "CH4" : (3.1,3.5,3.3,250),
    "ICE" : (2.0,4.0,2.82,50),
    "SED" : (2.0,5.0,3.16,35),   
    }

    wave_min=filters[filter][0]
    wave_max=filters[filter][1]
    wave_cen=filters[filter][2]
    R=filters[filter][3]

    R=float(R)
    dxspectrum = int(np.ceil( np.log10(wave_max/wave_min)/np.log10(1.0+1.0/(R*2.0)) ))
    wave=np.zeros(dxspectrum)
    wave[0]=wave_min
    for i in range(dxspectrum-1):
    	wave[i+1]=wave[i]*(1+(1/(2*R)))

    lambdac=wave_cen*10000

    ABmag = mag

    fnu = 10**(-0.4*(ABmag + 48.60))                 # erg/s/cm^2/Hz
    flambda = fnu*Ang/((lambdac*Ang)**2/c)          #erg/s/cm^2/Ang# 
    E_phot = (h*c)/(lambdac*Ang) # erg
    flux_phot=flambda/E_phot                       #photons/s/cm^2/Ang
    flux_phot=flux_phot*10000                       #photons/s/cm^2/micron

    psf_cube=fits.getdata('keck_psf_cube_2.0_5.2_n341.fits')
    wave_samp_orig=np.linspace(1.95,5.2,341)
    psf_new=resample_3d_cube(psf_cube,wave_samp_orig,wave)

    for i in range(len(psf_new)):
        psf_new[i,:,:]=psf_new[i,:,:]/psf_new[i,:,:].sum()
    filter_trans = np.ones(dxspectrum)
    intfilter = integrate.trapezoid(filter_trans,wave) 

    spec,specwave=get_model_spectra('vega')

    spec_temp=resample_data(specwave,spec,wave)
    spec_temp=spec_temp/spec_temp.sum()

    flux_phot_sec=flux_phot*area*intfilter

    spec_temp=spec_temp*flux_phot_sec

    for i in range(len(psf_new)):
        psf_new[i,:,:]=psf_new[i,:,:]*spec_temp[i]

    signal_cube=psf_new
    
    sky_data=np.loadtxt('mk_skybg_zm_16_15_ph.dat.txt')
    wave_sky=sky_data[:,0]*0.001 #microns
    flux_phot_sky=sky_data[:,1]*1000/10000 # photons/sec/cm2/microns

    flux_sky=resample_data(wave_sky, flux_phot_sky, wave) #photons/sec/cm2/microns/arcsec2

    flux_sky=flux_sky*(wave[1]-wave[0])*area*(scale**2)

    background_cube=np.zeros(signal_cube.shape)
    for i in range(108):
        for j in range(108):
            background_cube[:,i,j]=flux_sky

    readnoise = readnoise**2.0/itime
    noise_total=background_cube+readnoise
    signal = signal_cube*np.sqrt(itime*nframes)
    noise = np.sqrt(signal_cube+noise_total)
    snr=signal/noise
    return snr,wave
