"""Calculation of quasi-Bessel beams after an ideal axicon

Python module for the calculation of quasi-Bessel beams generated via an ideal Gaussian beam and an ideal axicon as per
P. Wu et al., Theoretical analysis of a quasi-Bessel beam for laser ablation, Photon. Res. / Vol. 2, No. 3 / June 2014
http://dx.doi.org/10.1364/PRJ.2.000082

Implemented by Oskar Hofmann, 2021    
"""

import numpy as np
import numpy.typing as ntp
from math import pi
from scipy.special import jv # v-th order Bessel function of the first kind


def radial_intensity(rho: ntp.ArrayLike, z: ntp.ArrayLike, w: float = 1, I_0: float = 1, n: float = 1.5, alpha: float = 5, wavelength: float = 1E-6) -> ntp.ArrayLike:
    """Create the intensity profile of a quasi-Bessel beam behind an ideal axicon illuminated by an ideal Gaussian beam

    Args:
        rho (numpy.ndarray): radial distance from the optical axis
        z (numpy.ndarray): position from the axicon along the optical axis. Units are arbitrary but must be consistent with units for rho.
        w (float, optional): beam waist (radius) of the Gaussian beam before the axicon. Units are arbitrary but must be consistent with units for rho. Defaults to 1. 
        I_0 (float, optional): Peak intensity of the Gaussian beam before the axicon. Units are arbitrary and determine the units for the outout. Defaults to 1. 
        n (float, optional): Refractive index of the axicon at the chosen wavelength. Defaults to 1.5.
        alpha (float, optional): Base angle of the axicon in degrees. Defaults to 5.
        wavelength (float, optional): Wavelength of the beam. Units are arbitrary but must be consistent with units for rho. Defaults to 1E-6.

    Returns:
        numpy.ndarray: radial intensity of the quasi-Bessel beam at the given (rho, z) positions
    """
    
    alpha = alpha*pi/180 #convert from degree to rad
    beta = 2*pi*(n-1)*alpha/wavelength

    z_0 = (n-1)*alpha*z/w

    def F1(x: ntp.ArrayLike) -> ntp.ArrayLike: 
        #equation 2
        return np.sqrt(z_0+x)*np.exp(-(z_0+x)**2)
    
    def F2(x: ntp.ArrayLike) -> ntp.ArrayLike: 
        #equation 3

        #implement heaviside function via np.maximum of the argument and 0, otherwise np.sqrt() would lead to NaN for negative values (NaN*0 = NaN)
        return np.sqrt(np.maximum(z_0-x,0))*np.exp(-(z_0-x)**2)
        
    #equation 1
    intensity = I_0*pi*beta*w/2*(((F1(rho/w)+F2(rho/w))*jv(0,rho*beta))**2+((F1(rho/w)-F2(rho/w))*jv(1,rho*beta))**2)

    return intensity


def dof(w: float, n: float = 1.5, alpha: float = 5) -> float:
    """Calculate the depth of field of a quasi-Bessel beam behind an ideal axicon illuminated by an ideal Gaussian beam

    Args:
        w (float): beam waist (radius) of the Gaussian beam before the axicon. Units are arbitrary and determine the units for the outout.
        n (float, optional): Refractive index of the axicon at the chosen wavelength. Defaults to 1.5.
        alpha (float, optional):  Base angle of the axicon in degrees. Defaults to 5.

    Returns:
        float: depth of field of the quasi-Bessel beam as defined in the aforementioned paper
    """
    alpha = alpha*pi/180 #convert from degree to rad
    
    #z_max defined below equation 5
    return w/((n-1)*alpha)



def beam_diameter(wavelength: float = 1E-6, n: float = 1.5, alpha: float = 5) -> float:
    """Calculate the approx. beam diameter of a quasi-Bessel beam behind an ideal axicon illuminated by an ideal Gaussian beam defined via the first minimum

    Args:
        wavelength (float, optional): Wavelength of the beam. Units are arbitrary but must be consistent with units for rho. Defaults to 1E-6.
        n (float, optional): Refractive index of the axicon at the chosen wavelength. Defaults to 1.5.
        alpha (float, optional):  Base angle of the axicon in degrees. Defaults to 5.

    Returns:
        float: beam diamaeter of the QBB approximated via the first root of Bessel function J_0 (equation 11)
    """    

    alpha = alpha*pi/180 #convert from degree to rad

    #equation 11
    return 1.2*wavelength/(pi*(n-1)*alpha)