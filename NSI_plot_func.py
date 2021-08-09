def get_gauss_SI_for_NSI(WAVELENGTH):
    '''
     get values for gaussian curve for NSI plot,
      returns x-values, then y-values - lists
     
    '''  
  
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

def get_doubles_SI_for_NSI(WAVELENGTH):
    '''
    Get values for doubles' curve for NSI plot,
    returns x-values, then y-values - as lists
    '''

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

def get_narayan_for_NSI():
    ''' 
    Narayan curve for NSI plot, returns x-values, y-values - as lists
    '''
    x_values = []
    separation_range = np.logspace(-1,0.5, 50)
    for i in range(len(separation_range)):
      x_values.append(separation_range[i])
    y_narayan = np.where(np.array(x_values) < 0.3, 1, (0.3/np.array(x_values))**(7/6.))
    return (x_values, y_narayan)

def get_nsi_plot(WAVELENGTH):
    '''
    Plots NSI, requires wavelenght consistent with previous functions called
    '''
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