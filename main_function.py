import numpy as np
#import scipy
from matplotlib import pyplot as plt
#%matplotlib inline
import random
from scipy.signal import hanning, tukey, boxcar
from scipy.spatial import distance
from scipy.ndimage import rotate
from math import sin, cos, radians, pi

from astropy.modeling import models
from astropy.modeling.models import Rotation2D, Gaussian2D

import pymc3 as pm

#!pip install iminuit # you can comment this out once it has run once
#import iminuit

#!pip install astropy
import astropy

from numpy.fft import ifft, fft, ifft2, fft2, fftshift, ifftshift, rfft, rfft2, fftfreq

from scipy.constants import au, c
z=au
k = 2*np.pi / 2.

def GaussHPBWVariationsLog(start,stop,step):
    '''
    Returns an array of values for hpbw in logspace
    '''      
    standardDev_range = np.logspace(start,stop, step)
    return (standardDev_range)

def get_gaussian_visibilities(theta_s, b, w=2): # theta_s = hpbw
    """
    produce 2D gaussian for
    - Gaussian FWHM theta_s (in arcseconds)
    - 2d array of baseline lengths b ( in m)
    - wavelength in m
    """
    A = np.radians(1)/3600.
    F = 2*np.sqrt(2*np.log(2))
    theta_r = A*theta_s
    b_hpbw = w/theta_r
    #print(b_hpbw)
    return np.exp(-0.5*(F*b/b_hpbw)**2)

def get_q_gaussian(wavelenght):
    '''
    Calculates a q for Gaussian case
    - wavelenght in (m)
    Return: q
    '''

    r_F = np.sqrt((wavelenght*z)/(2*pi)) # np.sqrt((z/k)=(wavelenght*z)/(2*pi))
    q_new_values = np.linspace(10/(r_F), 0.001/r_F, 50)

    arr = np.sqrt(q_new_values[:, None]**2 + q_new_values[None, :]**2)
    arr1 = np.empty((50,50))
    for x in range(50):
      for y in range(50):
        arr1[x,y] = arr[-x,y]
    arr2 = np.concatenate((arr, arr1))

    arr3 = np.empty((100,50))
    for x in range(100):
      for y in range(50):
        arr3[x,y] = arr2[x,-y]
    arr4 = np.concatenate((arr2, arr3), axis=1)

    arr2 = arr[::-1, :]
    arr3 = arr[:,::-1]
    arr4 = arr[::-1, ::-1]

    arr12 = np.concatenate((arr, arr2))
    arr34 = np.concatenate((arr3, arr4))
    arrF = np.concatenate((arr12, arr34), axis=1)
    return arrF

def get_q_doubles(wavelenght):
    '''
    Provides array of baselines for Double's case
    - wavelenght in (m)
    Return: arrF = q; qx = q with constant x for each column of y
    '''
    r_F = np.sqrt((wavelenght*z)/(2*pi))

    q_double = np.linspace(10/(r_F), 0.001/r_F, 50)

    q = np.hypot(q_double[None, :], q_double[:, None])                ## q - main for doubles

    qx = q_double[None, :] * np.ones([len(q_double), 1])

    arr = q

    arr2 = arr[::-1, :]

    arr2 = arr[::-1, :]
    arr3 = arr[:,::-1]
    arr4 = arr[::-1, ::-1]

    arr12 = np.concatenate((arr, arr2))


    arr34 = np.concatenate((arr3, arr4))

    arrF = np.concatenate((arr12, arr34), axis=1)
    return arrF, qx












def get_gaussian_ps(HalfPowerBeamWidth, Orientation):
    '''
    Provides a plot of power spectrum for Gaussian source
    - hpbw
    - angle of orientation

    '''

    arrF = get_q_gaussian(2)
    array_x = np.arange(-50,50, 1)

    standardDev = HalfPowerBeamWidth
    WAVELENGTH=2
    A = np.radians(1)/3600.
    F = 2*np.sqrt(2*np.log(2)) 
    C = 2*A*F*np.pi
    q_gauss = arrF                           ###           setting q
    b = q_gauss*z/k

    Vis = get_gaussian_visibilities(standardDev, b, 2)
    FF = FresnelFilterNew(2, q_gauss)
    #SI = ScintilationIndexNew(Vis,FF)                        ### Only need the positive values for the Power Spectrum
    powerSpectrum_p = np.sum(FF, axis=0)
    powerSpectrum_point = []
    for i in powerSpectrum_p: 
      if i > 0:
        powerSpectrum_point.append(i)
    plt.loglog(array_x, powerSpectrum_point)
    powerSpectrum_s = np.sum(Vis*FF, axis=0)
    powerSpectrum_source = []
    for i in powerSpectrum_s: 
      if i > 0:
        powerSpectrum_source.append(i)
    plt.loglog(array_x, powerSpectrum_source)
    plt.show()
    plt.plot(array_x, powerSpectrum_source)
    plt.show()


def visibilityFT(imageArray):
    '''
    Transforms image plane data into visibility plane
    Returns Visibility
    '''
    imageArrayFFT = np.fft.fft2(imageArray)
    imageArrayFFTabs = np.abs(imageArrayFFT)
    b = np.fft.fftshift(imageArrayFFTabs)
      # Shift peak to be at 1
    if b.max() == 0:
        b=b
    else:
        b= b/b.max()
    return b

def getQ(radius, wavelenght, x=5):
    '''
    Obtains q-array, x - n.pixels at first r(r-in pixels) / first zero
    ' 
    '''
    z = (1.496*(10)**11) # metres  =  1 AU
    r_F = np.sqrt((wavelenght*z)/(2*pi)) # np.sqrt((z/k)=(wavelenght*z)/(2*pi))
    q = x*radius/radius.max()/r_F

    return q