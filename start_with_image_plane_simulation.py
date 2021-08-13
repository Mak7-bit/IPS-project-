
from matplotlib import pyplot as plt
from image_plane import *


## Visibility Test ##

y_values = []
x_values = []
s=100
imageArray = np.zeros((s, s))
standardDev_range = np.arange(0,3.2,0.2)
for i in range(len(standardDev_range)):
  standardDev = standardDev_range[i]
  Vis = createCircular_gaussian(standardDev, imageArray)
  plt.imshow(Vis)
  plt.show()



## Test: Adding a gaussian to an image plane and plotting it ##


import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models
from astropy.modeling.models import Rotation2D, Gaussian2D

s = 100 #512 so far the source is a single pixel, so assuming the size of the image is 100x100 pixels 
r = s # rows
c = s # columns
t = 10 # diffusion size 
a = np.zeros((r, c)) # empty array to plot

#x1 = np.linspace(-5,5,100)
#x2 = np.linspace(-5,5,100)
#print(x1)
x2=50
y2=50
x,y = np.mgrid[-5:5:0.1,-5:5:0.1]
#plt.figure(figsize=(12, 12))
standardDev_range = np.arange(0,3.2,0.2)
for i in range(len(standardDev_range)):
  standardDev = standardDev_range[i]
  print(f"standardDev={standardDev}")
  g = Gaussian2D(20,0,0,standardDev,standardDev)
  array = np.nan_to_num(g(x,y))
  print(array)
  print(f"Dimensions of array ={array.shape}")
  #a = a(row,col)+ array(row,col)
  plt.imshow(array)#, extent=(-50,50,-50,50))
  plt.show()

  gaussGridSize = 10
  gaussGridStep = 0.5
  coordXGauss,coordYGauss = np.mgrid[0:gaussGridSize:gaussGridStep,0:gaussGridSize:gaussGridStep]
  theta = 30
  theta_rad = radians(theta)
  g = Gaussian2D(20,5,5,standardDev,standardDev,theta_rad)
# Slicing a array and adding Objects into the picture
  slXstart = x2 - gaussGridSize
  slXstop = x2 + gaussGridSize
  slYstart = y2 - gaussGridSize
  slYstop = y2 + gaussGridSize
  array = np.nan_to_num(g(coordXGauss,coordYGauss))
  #array = np.nan_to_num(g(x,y))
  #print(g(coordXGauss,coordYGauss))
  #if standardDev == 0:
  #  print('Cought it')
  #  a = np.nan_to_num(a)
  #else:
  a[slXstart:slXstop,slYstart:slYstop] += array
  #a = np.nan_to_num(a)
  #plt.imshow(g(coordXGauss,coordYGauss))#, extent=(-50,50,-50,50))
  #plt.show()
  plt.imshow(a)#, extent=(-50,50,-50,50))
  plt.show()


#plt.imshow(g(x, y), origin='lower')
#plt.colorbar()


#a[56:76,9:29] += g(x,y)

plt.figure(figsize=(12, 12))
#plt.imshow(a)
#plt.colorbar()
#theta = 30
#theta_rad = radians(theta)
#x,y = np.mgrid[-30:30:0.5,-30:30:0.5]
#g = Gaussian2D(10,0,0,1.3,1.3,theta_rad)
plt.imshow(circular_gaussian(1.2))



        ## Simulation Test ##
          # Output: 2 Gaussian sources, their coordinates, image Plane, FT #



s = 100 #512 so far the source is a single pixel, so assuming the size of the image is 100x100 pixels 
r = s # rows
c = s # columns
t = 10 # diffusion size 
a = np.zeros((s, s)) # empty array to plot
a11 = np.zeros((s, s)) # arrays for diffusion; NAME: first number = Source number(1st or 2nd); 2nd Number is a version of the array
a12 = np.zeros((s, s))
a21 = np.zeros((s, s))
a22 = np.zeros((s, s))
#pix_per_mas = 1.0 # image size / number of pixels  : pixel_size = image size(in arcsec) / number_of_pixels.    ASSUMING THE PICTURE SIZE TO BE 100 mas : 1mas/px or 1px/mas
#pix_per_beam = 1.0 # usually x3 of pix_size  ???

x1 = 50#int(input('Enter x coordinate of Source One: '))
y1 = 50#int(input('Enter y coordinate of Source One: '))
separation_dist = 9#int(input('Enter distance between 2 sources: '))
theta = 0#int(input('Enter an angle of orientation of double(in degrees; counting from eastward dir.): '))

x2,y2 = point_pos(x1,y1,separation_dist,theta)

#d = 0 #int(input('Enter distance between 2 sources: '))
#theta = 0 #int(input('Enter an angle of orientation of double(in degrees; counting from eastward dir.): '))
ratio = 1#float(input('Enter a ratio of brightnesses of the Double Source'))
#x2,y2 = point_pos(x1,y1,d,theta)
print('Coordiantes of Source 1: ', x1,y1)
print('Coordiantes of Source 2: ', x2,y2)
#LS1, LS2 = brightnessRatio(ratio)
#LS2 = 100
#LS1 = LS2

# Gaussian Objects

### 1 
gaussGridSize = 10
gaussGridStep = 0.5
coordXGauss,coordYGauss = np.mgrid[0:gaussGridSize:gaussGridStep,0:gaussGridSize:gaussGridStep]
theta = 0
theta_rad = radians(theta)
g = Gaussian2D(20,5,5,0.8,0.8,theta_rad)
# Slicing a array and adding Objects into the picture
slXstart = x1 - gaussGridSize
slXstop = x1 + gaussGridSize
slYstart = y1 - gaussGridSize
slYstop = y1 + gaussGridSize
a[slXstart:slXstop,slYstart:slYstop] += g(coordXGauss,coordYGauss)
plt.imshow(g(coordXGauss,coordYGauss))
plt.show()
print(a)
print(a.max())

### 2 
gaussGridSize = 10
#G1 = -gaussGridSize/2
#G2 = gaussGridSize/2
gaussGridStep = 0.5
coordXGauss,coordYGauss = np.mgrid[0:gaussGridSize:gaussGridStep,0:gaussGridSize:gaussGridStep]
theta = 30
theta_rad = radians(theta)
g = Gaussian2D(20,5,5,0.8,0.8,theta_rad)
# Slicing a array and adding Objects into the picture
slXstart = x2 - gaussGridSize
slXstop = x2 + gaussGridSize
slYstart = y2 - gaussGridSize
slYstop = y2 + gaussGridSize
a[slXstart:slXstop,slYstart:slYstop] += g(coordXGauss,coordYGauss)

# noise
#for noisesize in range(100):
  #x = random.randint(1,s-1)
  #y = random.randint(1,s-1)
  #noisevalue = np.random.normal(0,5)
  #if noisevalue >= 0:
    #a[x,y] = noisevalue
  #else:
    #noisesize += 1
plt.figure(figsize=(12, 12))
plt.title('Resolved in 100x100 px, bilinear smoothening applied')
plt.xlabel('Right Ascension, per metre ') # = 1/3600 *10^-6 degrees  # ("RA / mas")
plt.ylabel('Declination , per metre ')  # ("Decl. / mas")
plt.imshow(a, interpolation='bilinear', extent=(-50,50,-50,50))
plt.colorbar()

plt.figure(figsize=(12, 12))
b = np.fft.fft2(a)
b = np.abs(b)
plt.title(' FFT of the Source, DC frequency in the middle')
plt.imshow(fftshift(b),interpolation='bilinear',extent=(-50,50,-50,50))
plt.colorbar()



      ## Normalisation ##
plt.figure(figsize=(12, 12))
b = np.fft.fft2(a)
b_noshift = np.abs(b)
b = np.fft.fftshift(b_noshift)
b=np.abs(b)
b=b/b.max()
#b1 = np.fft.ifft2(ifftshift(b))
plt.xlabel('metres^-1 ')
plt.ylabel('metres^-1 ')
plt.imshow(np.abs(b),interpolation='bilinear', extent=(-50,50,-50,50))

#print(a.shape)
#N = a.shape[-1]
#print(N)
#half_width = N/2/pix_per_mas
#a1 = a[::-1] # - inverse of the array 
#plt.imshow(a1, extent=(-half_width,half_width,-half_width,half_width))
plt.colorbar() 
plt.show()

b_noshift = b_noshift/b_noshift.max()
plt.figure(figsize=(12, 12))
plt.imshow(np.abs(b_noshift),interpolation='bilinear', extent=(-50,50,-50,50))
plt.colorbar() 
plt.show()
np.set_printoptions(threshold=False)
print(b_noshift)
print(b_noshift.min())
for r in range(1,s-1):
  c=1
  if b_noshift[r,c] == 0:
    print('Coordinates: ',r,c)
## FOR SOME REASON MINE IMAGE IS IN FACT GETTING REVERSED, UNLIKE JOHNS


        ## Tukey Filter ##


ty = tukey(100, alpha = 0.25) # 0.15 is the value taken from original code
tukey2 = ty[:, np.newaxis]*ty[np.newaxis, :]
plt.figure(figsize=(12, 12))
plt.imshow(a[::-1]*tukey2, extent=(-50,50,-50,50)) # theoretically, the image should be cleanned up from lower band wavelenghts
plt.show()

theta_large = 1/3600 # Num of pix * angle of each pixel
plot_uv_lim = theta_large    # 1m / 1mas
#plot_uv_lim = 0.01 # wHAT IS IT ACTUALLY REPRESENTS ???
#sl1 = M//2-M*(plot_uv_lim*pix_per_beam/pix_per_mas) # =490
#sl2 = int(M//2+M*(plot_uv_lim*pix_per_beam/pix_per_mas)) # =510



####
# Distance Clalcs


M = a.shape[0]
coords_x, coords_y = np.meshgrid(np.arange(s), np.arange(s))
#print(f"coords_x={coords_x}")
r = np.float32(np.hypot(coords_x-(s-1)/2, coords_y-(s-1)/2))
#print(f"M={M}")
plt.imshow(r)#, extent=(-plot_uv_lim, plot_uv_lim, -plot_uv_lim, plot_uv_lim))
plt.colorbar()
plt.show()
print(r)

##print(r.ndim)
rmax = r.max()
print(rmax)
#print(r.min())


theta_small = 1/3600/100 # lenght per pixel = 100metres/100 pixels
wavelenght =  2
theta_smth_rad = separation_dist * theta_small      ####10**(-6) # radians  = 0.2arcseconds
z = (1.496*(10)**11) # metres  =  1 AU
x = 5  # n.pixels at first r(r- in pixels)
q_r = pi*wavelenght/(2*theta_smth_rad*x)
#print(q_r)
q = q_r * r
k = 2*pi/ wavelenght
ff_jp = np.nan_to_num(np.sin(q**2*z/(2*k))**2)
ff_jp_Sum = np.sum(ff_jp, axis=1)
plt.figure(figsize=(12, 12))
#plt.imshow(ff_jp,interpolation='bilinear', extent=(-50,50,-50,50))
plt.plot(ff_jp_Sum)
#plt.colorbar() 
plt.show()
d_capital = pi*wavelenght/(2*theta_smth_rad)

z_k = z/k # z/k)

#theta_small = 1/3600/100 # lenght per pixel = 100metres/100 pixels
theta_large = 1/3600 # Num of pix * angle of each pixel    ; 1 ARCsecond
#wavelenght = 2 # metres
d = wavelenght / theta_large # Resolution
D = wavelenght / theta_small
arcsec_per_pix = 0.03 # since 3 arcseconds across



