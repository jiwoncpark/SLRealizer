#==============================================

import numpy as np
import pandas
import desc.slrealizer

#==============================================

def generate_data(curr_lens, curr_obs):    
    """
    Given the current observation detail and current lensed system, generate a mock catalog.
    """

    ###
    first_moment_err = 0.01 # unit : percentage
    first_moment_err_std = 0.005
    flux_err = 0.01 # unit : percentage                                                                            
    flux_err_std = 0.005
    second_moment_err = 0.05 # unit : percentage
    second_moment_err_std = 0.01
    ###

    processed_image = desc.slrealizer.plot_all_objects(curr_obs, curr_lens, False)
    flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend_v2(processed_image)
    flux_err_calc = flux * noissify_data(flux_err, flux_err_std)
    first_moment_x_err_calc = noissify_data(first_moment_err, first_moment_err_std) * first_moment_x
    first_moment_y_err_calc = noissify_data(first_moment_err, first_moment_err_std) * first_moment_y
    MJD, filter, PSF_HWHM, sky = curr_obs[0], curr_obs[1], curr_obs[2], curr_obs[3]
    I_xx, I_xy, I_yy = covariance_matrix[0][0], covariance_matrix[0][1], covariance_matrix[1][1]
    I_xx_err_calc = I_xx * noissify_data(second_moment_err, second_moment_err_std)
    I_xy_err_calc = I_xy * noissify_data(second_moment_err, second_moment_err_std)
    I_yy_err_calc = I_yy * noissify_data(second_moment_err, second_moment_err_std)
    lensID = curr_lens[0]['LENSID']
    return np.array([MJD, filter, first_moment_x, first_moment_x_err_calc, first_moment_y, first_moment_y_err_calc, flux, flux_err_calc, I_xx, I_xx_err_calc, I_yy, I_yy_err_calc, I_xy, I_xy_err_calc, PSF_HWHM, sky, lensID])

def noissify_data(mean, stdev):
    """
    Given a mean and a standard deviation of a measurement, add a noise to the data
    """
    return np.random.normal(loc=mean, scale=stdev)
