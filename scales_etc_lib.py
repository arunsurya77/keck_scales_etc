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
from scipy.interpolate import interp1d

def resample_3d_cube(cube, x_array, x_array_new):
    """
    Resample a 3D cube array along its first dimension based on a new x_array_new.

    Parameters:
        cube (np.ndarray): 3D numpy array of shape [341, 108, 108].
        x_array (np.ndarray): Original x values corresponding to the first dimension of the cube.
        x_array_new (np.ndarray): New x values within the range of x_array.

    Returns:
        np.ndarray: Resampled 3D numpy array with shape [len(x_array_new), 108, 108].
    """
    if cube.shape[0] != len(x_array):
        raise ValueError("The length of x_array must match the first dimension of the cube.")

    if np.any(x_array_new < x_array.min()) or np.any(x_array_new > x_array.max()):
        raise ValueError("x_array_new values must be within the range of x_array.")

    # Interpolator function for each 2D frame
    interpolator = interp1d(x_array, cube, axis=0, kind='linear', bounds_error=True)

    # Resample the cube based on the new x_array_new
    resampled_cube = interpolator(x_array_new)

    return resampled_cube
    
def resample_data(x, y, new_x):
    """
    Resample x and y data to a new x array and compute the corresponding y values.

    Parameters:
        x (array-like): Original x data.
        y (array-like): Original y data corresponding to x.
        new_x (array-like): New x array to resample data.

    Returns:
        new_y (numpy.ndarray): Resampled y values corresponding to new_x.
    """
    # Ensure inputs are numpy arrays
    x = np.array(x)
    y = np.array(y)
    new_x = np.array(new_x)

    # Perform linear interpolation
    new_y = np.interp(new_x, x, y)

    return new_y

def get_model_spectra(string):
    print(string + 'loaded' )
    ext=0
    spec_file = os.path.expanduser('vega_all.fits')
    pf = fits.open(spec_file)
    spec = pf[ext].data
    head = pf[ext].header
    cdelt1 = head["cdelt1"]
    crval1 = head["crval1"]
    nelem = spec.shape[0]
    specwave = (np.arange(nelem))*cdelt1 + crval1  # Angstrom
    specwave /= 1e4 # -> microns
    return spec,specwave
