"""
All the constants should be accessed through the methods here, so that they would not be confused.
"""

from astropy.coordinates import SkyCoord

def return_zeropoint():
    """
    returns zeropoint for the SDSS magnitude system
    """
    return 22.5

def return_obs_RA_DEC():
    """
    Return OM10's lensed system's position
    """
    return SkyCoord('03h 32m 30s', '10d 00m 24s')

def get_flux_err():
    return 0.01 # unit : percentage

def get_flux_err_std():
    return 0.005 # unit : percentage

def get_first_moment_err():
    return 0.01 # unit : percentage

def get_first_moment_err_std():
    return 0.005

def get_second_moment_err():
    return 0.05

def get_second_moment_err_std():
    return 0.01

def get_x_min():
    return -5.0

def get_y_min():
    return -5.0

def get_x_max():
    return 5.0

def get_y_max():
    return 5.0

def get_distance():
    return 0.2

def get_pixel_arcsec_conversion():
    # based on the distance, x_max, x_min, y_max, and y_min, we can know how many arcsec corresponds to one pixel
    return 0.1 # fix later : it should be the same as get_distance()
