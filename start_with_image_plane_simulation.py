def point_pos(x1, y1, d, theta):
    """
    Calculates position of Source 2 for a Double,
    taking x and y coordinates of Source 1,
    separation and orientation angle
    - theta in radians
    - distance in arcseconds
    """ 
    theta_rad = radians(theta) + pi/2 #pi/2 - radians(theta)          Applicable for visib
    x2 = int(x1 + d*cos(theta_rad))
    y2 = int(y1 + d*sin(theta_rad))
    return x2,y2
    
def brightnessRatio(ratio):
    '''
     A fucntion that defines brightness for 2 sources, using random distribution and brigh
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
    Creates a source image using gaussian distribution and plugs it into an array that rep
    - amplitude
    - standard Deviation 
    - x coordinate of source
    - y coordinate
    - array in image plane
    '''
    gaussGridSize = 10                                                  ##      Guassian i
    gaussGridStep = 0.5
    coordXGauss,coordYGauss = np.mgrid[0:gaussGridSize:gaussGridStep,0:gaussGridSize:gauss
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

def circular_gaussian(stddev): 
    '''       
    Creates circular gaussian on a small grid 
    returns array to be plt.imshow-ed
    '''

    #theta = 30
    #theta_rad = radians(theta)
    x,y = np.mgrid[-30:30:0.5,-30:30:0.5]
    g = Gaussian2D(10,0,0,stddev,stddev)#,theta_rad)
    return (g(x,y))

def createCircular_gaussian(standardDev,imageArray):
    '''        
    Creates a circular gaussia
    returns Image plane
    '''
    gaussGridSize = 10
    gaussGridStep = 0.5
    x1 = 50
    y1 = 50
    coordXGauss,coordYGauss = np.mgrid[0:gaussGridSize:gaussGridStep,0:gaussGridSize:gaussGr
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

def createDouble(x1,y1,x2,y2,image):
    '''
    Creates 2 unresolved sources using image plane a
    '''
    image[x1,y1] = 20
    image[x2,y2] = 20
    return image


def rCalcs(imageArray, size):   ## Calculates r-distances array and normalizes it. 
  M = imageArray.shape[0]
  coords_x, coords_y = np.meshgrid(np.arange(size), np.arange(size))
  r = np.float32(np.hypot(coords_x-(size-1)/2, coords_y-(size-1)/2))
  rmax = r.max()
  return (r,rmax)