def xy_rotate_scale(offset, sigma_a, sigma_b, pa):
    """
    return offset in sigma from an elliptical Gaussian centred on (x0, y0). 
    offset:           numpy array of coordinates ( e.g. array([[x1, y1],[x2, y2]]))
    sigma_a, sigma_b: major and minor axes of gaussian
    pa:               position angle of guassian measured from y axis
                      defined by sigma_a, sigma_b and pa (in degrees),
    return normalised offset in standard deviations (sigma) for each coordinate
    xy_offset_to_gaussian(((1.0, 0.0),), 1.0, 0.5, 0.0)
    >>> array([ 1.])
    """

    pa = pa*np.pi/180
    
    sigma_array = np.zeros(offset[..., 0].shape, dtype=complex)
    sigma_array.real = offset[..., 0]
    sigma_array.imag = offset[..., 1]
    rot = complex(np.cos(-pa), np.sin(-pa))
    sigma_array *= rot
    sigma_array.real /= sigma_a
    sigma_array.imag /= sigma_b
    final = np.abs(sigma_array).astype(float)
    return final

def gauss_for_rotated(offset_sigma):
    '''
    Auxilary function that does not require normalisation 
    '''
    F = 2* np.sqrt(2*np.log(2))
    return np.exp(-0.5*(F*offset_sigma)**2)

# creating the baselines q array
F = 2*np.sqrt(2*np.log(2))
wavelenght =2
r_F = np.sqrt((wavelenght*z)/(2*pi)) # np.sqrt((z/k)=(wavel
q_new_values = np.linspace(-10/r_F, 10/r_F, 50) #np.linspac
qx, qy = np.meshgrid(q_new_values, q_new_values)
q_hypot = np.stack((qx, qy), axis=2)
offset = q_hypot




wavelenght =2
r_F = np.sqrt((wavelenght*z)/(2*pi)) # np.sqrt((z/k)=(wavelenght*z)/(2*pi))
q_new_values = np.linspace(-10/r_F, 10/r_F, 50) #np.linspace(10/(r_F), 0.001/r_F, 50)

qx, qy = np.meshgrid(q_new_values, q_new_values)
q_hypot = np.stack((qx, qy), axis=2)

offset = q_hypot # is the offset

offset_sigma = rotatedGaussian(q_hypot, 0.00002, 0.00002/2, 30)



print('sigma_a: ', sigma_a,' , sigma_b: ',sigma_b)
plt.figure()
plt.imshow(gauss_for_rotated(offset_sigma))
plt.colorbar()