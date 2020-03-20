'''
File containing fit functions for distributions/resonances
'''
import numpy as np


def gaussFit(x,x0,A,sig):
    """
    Normalized gaussian that is used to fit resonance:
        x   - position
        x0  - centroid position
        A   - peak area
        sig - width
    """
    coef = A / (sig*np.sqrt(2*np.pi))
    gaussPart = np.exp(-.5*((x-x0)/sig)**2)

    return coef * gaussPart

def gaussFitPlusBack(x,x0,A,sig,c0,c1):
    """
    Normalized gaussian that is used to fit resonance:
        x   - position
        x0  - centroid position
        A   - peak area
        sig - width

    Linear background:
        c0 + c1*x
    """
    coef = A / (sig*np.sqrt(2*np.pi))
    gaussPart = np.exp(-.5*(x-x0)*(x-x0)/(sig*sig))
    back = c0 + c1*x

    return (coef * gaussPart) + back

def lorenFit(x,x0,A,gam):
    """
    Pseudo-Lorentz/Cauchy distribution
        x   - position
        x0  - centroid position
        gam - HWHM (half width half max)
        A   - peak area
    """
    coef = A*.5*gam/np.pi
    lorenPart = ( (x-x0)**2 + .25*gam**2)**(-1)

    return coef * lorenPart

def lorenFitPlusBack(x,x0,A,gam,c0,c1):
    """
    Pseudo-Lorentz/Cauchy distribution
        x   - position
        x0  - centroid position
        gam - HWHM (half width half max)
        A   - peak area

    Linear background:
        c0 + c1*x
    """
    coef = A*.5*gam/np.pi
    lorenPart = ( (x-x0)*(x-x0) + .25*gam*gam)**(-1)
    back = c0 + c1*x

    return (coef * lorenPart) + back
