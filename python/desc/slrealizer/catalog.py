import numpy as np
import pandas
import desc.slrealizer

#  MJD, filter, ra, ra_err, dec, dec_err, flux, flux_err, qxx, qxx_err, qyy, qyy_err, qxy, qxy_err, psf_sigma, sky
def make_lens_catalog(curr_lens, curr_obs):
    data = generate_data(curr_obs, curr_lens)
    df = pandas.DataFrame(columns=['MJD', 'filter', 'x', 'x_com_err', 'y', 'y_com_err', 'flux', 'flux_err', 'qxx', 'qxx_err', 'qyy', 'qyy_err', 'qxy', 'qxy_err', 'psf_sigma', 'sky'])
    df.loc[len(df)]= data
    print(df)
    print('saving the table with the name catalog.csv. Check your data folder (../../../data/)')
    df.to_csv('../../../data/catalog.csv')
    #df = pandas.DataFrame(columns=['MJD', 'filter', 'x', 'x_com_err', 'y', 'y_com_err', 'flux', 'flux_err', 'qxx', 'qxx_err', 'qyy', 'qyy_err', 'qxy', 'qxy_err', 'psf_sigma', 'sky'])
    #df.append(data)

#obsHistIDexpMJDfilterFWHMefffiveSigmaDepth
def generate_data(curr_lens, curr_obs):    
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
    print curr_obs
    MJD = curr_obs[0]
    filter = curr_obs[1]
    sky = curr_obs[3]
    I_xx = covariance_matrix[0][0]
    I_xx_err_calc = I_xx * noissify_data(second_moment_err, second_moment_err_std)
    I_xy = covariance_matrix[0][1]
    I_xy_err_calc = I_xy * noissify_data(second_moment_err, second_moment_err_std)
    I_yy = covariance_matrix[1][1]
    I_yy_err_calc = I_yy * noissify_data(second_moment_err, second_moment_err_std)
    PSF_HWHM = curr_obs[2]
    print('calculation done')
    return np.array([MJD, filter, first_moment_x, first_moment_x_err_calc, first_moment_y, first_moment_y_err_calc, flux, flux_err_calc, I_xx, I_xx_err_calc, I_yy, I_yy_err_calc, I_xy, I_xy_err_calc, PSF_HWHM, sky])

def noissify_data(mean, stdev):
    return np.random.normal(loc=mean, scale=stdev)
