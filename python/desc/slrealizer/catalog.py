#==============================================

import numpy as np
import pandas
import desc.slrealizer
from fractions import Fraction
import math

#==============================================

"""
This file contains methods that help generate the toy catalog for the SLRealizer.
"""


def generate_data(curr_lens, curr_obs):    
    """
    Given the current observation detail and current lensed system, generate a mock catalog.
    """

    MJD, filter, PSF_HWHM, sky_mag = curr_obs[0], curr_obs[1], curr_obs[2], curr_obs[3]
    processed_image = desc.slrealizer.plot_all_objects(curr_obs, curr_lens)
    flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend(processed_image)
    flux += flux*noissify_data(desc.slrealizer.get_flux_err(), desc.slrealizer.get_flux_err_std())
    flux_err_calc = math.pow(10, (22.5 - sky_mag)/2.5)/5 # because Fb = 5 \sigma_b
    first_moment_x += noissify_data(desc.slrealizer.get_first_moment_err(), desc.slrealizer.get_first_moment_err_std()) * first_moment_x
    first_moment_y += noissify_data(desc.slrealizer.get_first_moment_err(), desc.slrealizer.get_first_moment_err_std()) * first_moment_y
    first_moment_x_err_calc, first_moment_y_err_calc = PSF_HWHM / math.pow(flux/flux_err_calc, 0.5), PSF_HWHM / math.pow(flux/flux_err_calc, 0.5)


    """
    REAL BUG: RA, DEC, 2nd moment error values are not yet valid
    """
    
    
    RA, RA_err, DEC, DEC_err = return_coordinate(first_moment_x, first_moment_y)
    I_xx, I_xy, I_yy = covariance_matrix[0][0], covariance_matrix[0][1], covariance_matrix[1][1]
    I_xx += I_xx * noissify_data(desc.slrealizer.get_second_moment_err(), desc.slrealizer.get_second_moment_err_std())
    I_xy += I_xy * noissify_data(desc.slrealizer.get_second_moment_err(), desc.slrealizer.get_second_moment_err_std())
    I_yy += I_yy * noissify_data(desc.slrealizer.get_second_moment_err(), desc.slrealizer.get_second_moment_err_std())
    size = I_xx + I_yy
    size_err = 0 # for now
    I_xx_err_calc, I_xy_err_calc, I_yy_err_calc = 0, 0, 0
    lensID = curr_lens[0]['LENSID']

    e = 0 # we do not have to calculate here, just for galsim
    return np.array([MJD, filter, RA, RA_err, DEC, DEC_err, first_moment_x, first_moment_x_err_calc, first_moment_y, first_moment_y_err_calc, flux, flux_err_calc, size, size_err, e, PSF_HWHM, sky_mag, lensID])

def noissify_data(mean, stdev):
    """
    Given a mean and a standard deviation of a measurement, add a noise to the data
    """
    return np.random.normal(loc=mean, scale=stdev)

def return_coordinate(first_moment_x, first_moment_y):
    """
    This method returns the RA and DEC values that take the x and y offset into account.
    RA and DEC values are provided in the observational history,
    and this method also assumes a random error in the measurement (+1/-1 deg)
    """

###                                                                         
    pos_err = 0.0 # unit : degree                                               
    pos_err_std = Fraction(1, 3) # Assuming that the three sigma is one degree, one sigma is 0.3333 degree                                                     
    ###      

    real_coordinate = desc.slrealizer.return_obs_RA_DEC()
    RA = real_coordinate.ra.deg
    DEC = real_coordinate.dec.deg
    # add the offset by the first moment
    RA += first_moment_x/(3600*np.cos(DEC))
    DEC += first_moment_y/3600
    # draw random number to set position in the FOV
    RA += np.random.uniform(-1.75, 1.75)
    DEC += np.random.uniform(-1.75, 1.75)
    RA += noissify_data(pos_err, pos_err_std)
    DEC += noissify_data(pos_err, pos_err_std)
    RA_err = 0
    DEC_err = 0
    return RA, RA_err, DEC, DEC_err

def return_mean_properties(lens_array):
    """
    returns the mean value of flux, 1st moment of x, 1st moment of y, qzz, qxy, qyy, flux_err, 1st moment x error, 1st moment y error, qxx error, qxy error, and qyy error
    """

    return lens_array['flux'].mean(), lens_array['x'].mean(), lens_array['y'].mean(), lens_array['size'].mean(), lens_array['flux_err'].mean(), lens_array['x_com_err'].mean(), lens_array['y_com_err'].mean(), lens_array['size_err'].mean(), lens_array['e'].mean()
