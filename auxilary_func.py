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
import iminuit

#!pip install astropy
import astropy

from numpy.fft import ifft, fft, ifft2, fft2, fftshift, ifftshift, rfft, rfft2, fftfreq

from scipy.constants import au, c


def point_pos(x1, y1, d, theta):
    """
    Calculates position of Source 2 for a Double,
    taking x and y coordinates of Source 1,
    separation and orientation angle
    - theta in radians
    - distance in arcseconds
    """ 
    theta_rad = radians(theta) + pi/2 #pi/2 - radians(theta)          Applicable for visibility plane 
    x2 = int(x1 + d*cos(theta_rad))
    y2 = int(y1 + d*sin(theta_rad))
    return x2,y2
    
def brightnessRatio(ratio):
    '''
     A fucntion that defines brightness for 2 sources, using random distribution and brightness ratio as input.
     - ratio input - float 
    '''    
    a = random.randint(0,2)
    if a > 0.8:
      LS1 = a*100
    else:
      LS1 = 100
    LS2 = int(ratio*LS1)
    return LS1, LS2

def createAGaussian(amplitude,standardDev, x, y, imageArray):
    '''
    Creates a source image using gaussian distribution and plugs it into an array that represents Image plane
    - amplitude
    - standard Deviation 
    - x coordinate of source
    - y coordinate
    - array in image plane
    '''
    gaussGridSize = 10                                                  ##      Guassian is circular
    gaussGridStep = 0.5
    coordXGauss,coordYGauss = np.mgrid[0:gaussGridSize:gaussGridStep,0:gaussGridSize:gaussGridStep]
    theta = 0
    theta_rad = radians(theta)
    g = Gaussian2D(amplitude,5,5,standardDev,standardDev,theta_rad)
  # Slicing a array and adding Objects into the picture
    slXstart = x1 - gaussGridSize
    slXstop = x1 + gaussGridSize
    slYstart = y1 - gaussGridSize
    slYstop = y1 + gaussGridSize
    array = np.nan_to_num(g(coordXGauss,coordYGauss))
    arr = np.nan_to_num(array)
    imageArray[slXstart:slXstop,slYstart:slYstop] += array
    return imageArray 

def circular_gaussian(stddev):        ###  Creates circular gaussian on a small grid 

  #theta = 30
  #theta_rad = radians(theta)
  x,y = np.mgrid[-30:30:0.5,-30:30:0.5]
  g = Gaussian2D(10,0,0,stddev,stddev)#,theta_rad)
  return (g(x,y))

def createCircular_gaussian(standardDev,imageArray):        ### Creates a circular gaussian in image plane using given SD, the object is placed in the centre
  gaussGridSize = 10
  gaussGridStep = 0.5
  x1 = 50
  y1 = 50
  coordXGauss,coordYGauss = np.mgrid[0:gaussGridSize:gaussGridStep,0:gaussGridSize:gaussGridStep]
  g = Gaussian2D(20,5,5,standardDev,standardDev)
  # Slicing a array and adding Objects into the picture
  slXstart = x1 - gaussGridSize
  slXstop = x1 + gaussGridSize
  slYstart = y1 - gaussGridSize
  slYstop = y1 + gaussGridSize
  array = np.nan_to_num(g(coordXGauss,coordYGauss))
  arr = np.nan_to_num(array)
  imageArray[slXstart:slXstop,slYstart:slYstop] += array
  if standardDev == 0:
    imageArray[50,50] = 20
  return imageArray 

def createDouble(x1,y1,x2,y2,image):  ### Creates 2 unresolved sources using image plane array and given coordinates
  # Creating Double source
  image[x1,y1] = 20
  image[x2,y2] = 20
  return image


def rCalcs(imageArray, size):   ## Calculates r-distances array and normalizes it. 
  M = imageArray.shape[0]
  coords_x, coords_y = np.meshgrid(np.arange(size), np.arange(size))
  r = np.float32(np.hypot(coords_x-(size-1)/2, coords_y-(size-1)/2))
  rmax = r.max()
  return (r,rmax)
 

def getQ(radius, wavelenght, x=5):       ### Obtains q-array, x - n.pixels at first r(r-in pixels) / first zero
  z = (1.496*(10)**11) # metres  =  1 AU
  r_F = np.sqrt((wavelenght*z)/(2*pi)) # np.sqrt((z/k)=(wavelenght*z)/(2*pi))
  q = x*radius/radius.max()/r_F

  return q
  


  #theta_small = 1/3600/100 # lenght per pixel : 10**(-6) # radians  = 0.2arcseconds
  #d_capital = pi*wavelenght/(2*theta_small)
  #theta_large = 1/3600 # Num of pix * angle of each pixel    ; 1 ARCsecond
  #d = wavelenght / theta_large # Resolution
  #D = wavelenght / theta_small
  #arcsec_per_pix = 0.03 # since 3 arcseconds across
  

def FresnelFilterNew(wavelenght, q_array): ## Fresnel Filter according to JP's formula, returns FF multiplied by Kolmogorov turbulence. Requires q as an argument 
        ### B - Constructed FF
  #x = 5  # n.pixels at first r(r- in pixels)
  z = (1.496*(10)**11) # metres  =  1 AU
  k = 2*pi / wavelenght
  z_k = z/k # z/k) fresnel scale   should be 1; whatever is inside the sin is = 1 so can calc the q
  r_F = np.sqrt((wavelenght*z)/(2*pi)) # np.sqrt((z/k)=(wavelenght*z)/(2*pi))
  r_F_u = 1/r_F
  #q = x*radius/radius.max()/r_F
  ff_jp = np.sin((q_array)**(2)*z/(2*k))**2
  turbulence = np.nan_to_num((q_array/r_F_u)**(-11/6.))
  ff_turb = ff_jp*turbulence
  ff_turb_sum = np.sum(ff_turb, axis=1)
  
  return (ff_turb)

def ScintilationIndexNew(visibility, fresnel_filter):  ## Takes in visibilit array and FF*Turbulence to calculate Scintilation Index
  ABC = visibility*fresnel_filter
  ABC_doubleIntegral = sum(list(map(sum, ABC))) # S. Index
  return (ABC_doubleIntegral)


def ScintilationIndexNormalisation(fresnel_filter, scintilation_index): ## Normalizes scintillation index by dividing it by an integral of FF and taking a sqaure root out of the result
      ## Scintilation Index: Normalised by ∫ ∫ BC === sqrt((∫ ∫ ABC)/(∫ ∫ BC))
  ff_turb_sum = np.sum(fresnel_filter, axis=1)
  BC = ff_turb_sum
  norm = np.sum(fresnel_filter)#, axis=0)
  siNorm = np.sqrt(scintilation_index/norm)
  return siNorm

def get_gaussian_visibilities(theta_s, b, w=2): # theta_s = hpbw
    """
    produce 1D gaussian for
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


def get_doubles_visibilities2D(theta_s):
    A = np.radians(1)/3600.
    '''
    produce double sources for
    - separation distance theta_s

    '''
     
def get_q_gaussian(wavelenght):
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
  r_F = np.sqrt((wavelenght*z)/(2*pi))

  q_double = np.linspace(10/(r_F), 0.001/r_F, 50)

  q = np.hypot(q_double[None, :], q_double[:, None])                ## q - main for doubles

  qx = q_double[None, :] * np.ones([len(q_double), 1])
  #plt.title('q')
  #plt.imshow(q)
  #plt.show()
  #plt.title('qx')
  #plt.imshow(qx)
  #plt.show()
  #plt.title('cos(4*pi*qx/100)')
  #plt.imshow(np.cos(4*np.pi*qx/100))
  #plt.show()

  arr = q

  arr2 = arr[::-1, :]

  arr2 = arr[::-1, :]
  arr3 = arr[:,::-1]
  arr4 = arr[::-1, ::-1]

  arr12 = np.concatenate((arr, arr2))
#plt.imshow(arr12)
#plt.show()
  arr34 = np.concatenate((arr3, arr4))
#plt.imshow(arr34)
#plt.show()
  arrF = np.concatenate((arr12, arr34), axis=1)
  return arrF, qx

def get_qx_final_double(qx):
  qx_neg_upper = qx[:,::-1]
  qx_neg_lower = qx_neg_upper
  qx_pos_lower = qx
  qx_neg = np.concatenate((qx_neg_upper, qx_neg_lower))
  qx_pos = np.concatenate((qx, qx_pos_lower))
  #plt.title('qx joint pos')
  #plt.imshow(qx_pos)
  #plt.show()
  qx_fin = np.concatenate((qx_neg, qx_pos), axis=1)
  return qx_fin

def get_double_ps(width, angle):     # width = separation dist in " , angle = orientation in degrees
  A = np.radians(1)/3600.
  F = 2*np.sqrt(2*np.log(2)) 
  C = 2*A*F*np.pi
  q, qx = get_q_doubles(2)
  qx_fin = get_qx_final_double(qx)

  qx_final_rot = rotate(qx_fin, angle, reshape=False)
  qx_rot_trunc = qx_final_rot[15:85,15:85]
  
  b = qx_rot_trunc * z/k
  q_reduced = (q[15:85,15:85])
  array_x = np.arange(-35,35, 1)

  separation = width
  Visibility = np.abs(np.cos(np.pi*b*separation*A))
  FF = FresnelFilterNew(2, q_reduced)                       #### q_reduced for roated and q for 100x100 unrotated array 
  #SI = ScintilationIndexNew(Visibility**2,FF)
  #SInorm = ScintilationIndexNormalisation(FF,SI)

  powerSpectrum_s = np.sum(Visibility**2*FF, axis=1)
  powerSpectrum_source = []
  for i in powerSpectrum_s: 
    if i > 0:
      powerSpectrum_source.append(i)
  powerSpectrum_p = np.sum(FF, axis=1)
  powerSpectrum_point = []
  for i in powerSpectrum_p: 
    if i > 0:
      powerSpectrum_point.append(i)
  plt.title(f'Separation: {separation:.2f}')
  plt.loglog(array_x, powerSpectrum_source)
  plt.loglog(array_x, powerSpectrum_point)
  plt.show()
  plt.title('PS')
  plt.plot(array_x, powerSpectrum_source)
  plt.show()


def get_gaussian_ps(HalfPowerBeamWidth, Orientation):

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
  
  #SInorm = ScintilationIndexNormalisation(FF,SI)


def get_gauss_SI_for_NSI(WAVELENGTH):  # get values for gaussian curve for NSI plot, returns x-values, then y-values
  
    arrF = get_q_gaussian(WAVELENGTH)
    standardDev_range = GaussHPBWVariationsLog(-1,0.5,50)
    array_x = np.arange(-50,50, 1)
    y_Gvalues = []
    x_Gvalues = []
    for i in range(len(standardDev_range)):
      standardDev = standardDev_range[i]
      A = np.radians(1)/3600.
      F = 2*np.sqrt(2*np.log(2)) 
      C = 2*A*F*np.pi
      q_gauss = arrF                           ###           setting q
      b = q_gauss*z/k
      Vis = get_gaussian_visibilities(standardDev, b, 2)
      FF = FresnelFilterNew(2, q_gauss)
      SI = ScintilationIndexNew(Vis,FF)
      powerSpectrum_point = np.sum(FF, axis=0)
      powerSpectrum_source = np.sum(Vis*FF, axis=0)
      #plt.loglog(array_x, powerSpectrum_source)
      #plt.show()
      #plt.plot(array_x, powerSpectrum_source)
      #plt.show()
  
      SInorm = ScintilationIndexNormalisation(FF,SI)

      y_Gvalues.append(SInorm)
      x_Gvalues.append(standardDev)
    return (x_Gvalues,y_Gvalues)

def get_doubles_SI_for_NSI(WAVELENGTH):   # get values for doubles' curve for NSI plot, returns x-values, then y-values

    A = np.pi/(180*3600) # arcseconds to radians
    r_F = np.sqrt((WAVELENGTH*z)/(2*pi))
    q, qx = get_q_doubles(2)
    qx_fin = get_qx_final_double(qx)
    angle = 0
    qx_final_rot = rotate(qx_fin, angle, reshape=False)
    qx_rot_trunc = qx_final_rot[15:85,15:85]
    b = qx_rot_trunc * z/k
    q_reduced = (q[15:85,15:85])
    array_x = np.arange(-35,35, 1)

    y_Dvalues = []
    x_Dvalues = []
    separation_range = np.logspace(-1,1, 50) 
    for i in range(len(separation_range)):
      separation = separation_range[i]
      Visibility = np.abs(np.cos(np.pi*b*separation*A))
      FF = FresnelFilterNew(2, q_reduced)                       #### q_reduced for roated and q for 100x100 unrotated array 
      SI = ScintilationIndexNew(Visibility**2,FF)
      SInorm = ScintilationIndexNormalisation(FF,SI)
      #powerSpectrum_source = np.sum(Visibility**2*FF, axis=1)
      #plt.title(f'Separation: {separation:.2f}')
      #plt.loglog(array_x, powerSpectrum_source)
      #plt.loglog(array_x, np.sum(FF, axis=1))
      #plt.show()
      #plt.title('PS')
      #plt.plot(array_x, powerSpectrum_source)
      #plt.show()
      y_Dvalues.append(SInorm)
      x_Dvalues.append(separation)
    
    return (x_Dvalues, y_Dvalues)


def get_narayan_for_NSI(): # Narayan curve for NSI plot, returns x-values, y-values
                           ### x-values for Narayan curve ###
    x_values = []
    separation_range = np.logspace(-1,0.5, 50)
    for i in range(len(separation_range)):
      x_values.append(separation_range[i])
    y_narayan = np.where(np.array(x_values) < 0.3, 1, (0.3/np.array(x_values))**(7/6.))
    return (x_values, y_narayan)

def get_nsi_plot(WAVELENGTH):
    gauss_x,gauss_y = get_gauss_SI_for_NSI(WAVELENGTH)
    print(gauss_x)
    print(gauss_y)
    double_x, double_y = get_doubles_SI_for_NSI(WAVELENGTH)
    narayan_x, narayan_y = get_narayan_for_NSI()

    plt.figure(figsize=(12,12))
    plt.title('qx')
 
    plt.loglog(gauss_x,gauss_y , 'o', label="Gaussian")
    plt.loglog(double_x, double_y, '+', label="Doubles")
    plt.loglog(narayan_x, narayan_y, label="Narayan")

    plt.axhline(y=(1/np.sqrt(2)))
    plt.xlabel("Gauss hpbw/Doubles separation")
    plt.ylabel("Normalised Scintilation Index")
    plt.legend()
    plt.show()

def xy_offset_to_gaussian(offset, sigma_a, sigma_b, pa=0.):
    """
    return offset in sigma from an elliptical Gaussian centred on (x0, y0). 
    offset:           numpy array of coordinates ( e.g. array([[x1, y1],[x2, y2]]))
    sigma_a, sigma_b: major and minor axes of gaussian
    pa:               position angle of guassian measured from y axis
                      defined by sigma_a, sigma_b and pa (in radians),
    return normalised offset in standard deviations (sigma) for each coordinate
    xy_offset_to_gaussian(((1.0, 0.0),), 1.0, 0.5, 0.0)
    >>> array([ 1.]) 
    """
    sigma_array = np.zeros(np.shape(offset)[0], dtype=complex)
    sigma_array.real = offset[..., 0]
    sigma_array.imag = offset[..., 1]
    
    rot = complex(cos(-pa), sin(-pa))
    
    sigma_array *= rot
    sigma_array.real /= sigma_a
    sigma_array.imag /= sigma_b

    return abs(sigma_array).astype(float)