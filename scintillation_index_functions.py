def FresnelFilterNew(wavelenght, q_array):
    ''' 
    Fresnel Filter according to JP's formula, returns FF multiplied by Kolmogorov turbulence. Requires q as an argument 
    '''
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

def ScintilationIndexNew(visibility, fresnel_filter):
    '''
    Takes in visibilit array and FF*Turbulence to calculate Scintilation Index
    '''
    ABC = visibility*fresnel_filter
    ABC_doubleIntegral = sum(list(map(sum, ABC))) # S. Index
    return (ABC_doubleIntegral)


def ScintilationIndexNormalisation(fresnel_filter, scintilation_index):
    ''' 
    Normalizes scintillation index by dividing it by an integral of FF and taking a sqaure root out of the result
    Scintilation Index: Normalised by ∫ ∫ BC === sqrt((∫ ∫ ABC)/(∫ ∫ BC))
    '''
    ff_turb_sum = np.sum(fresnel_filter, axis=1)
    BC = ff_turb_sum
    norm = np.sum(fresnel_filter)#, axis=0)
    siNorm = np.sqrt(scintilation_index/norm)
    return siNorm