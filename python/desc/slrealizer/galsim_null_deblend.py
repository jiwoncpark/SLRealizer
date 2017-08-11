#hello

import galsim
import numpy as np
import matplotlib.pyplot as plt
import desc.slrealizer
import math

def galsim_make_galsim_catalog(currLens, currObs):
    image = galsim_plot_all_objects(currLens, currObs)

def galsim_plot_all_objects(currLens, currObs):
    MJD, filter, PSF_HWHM, sky_mag = currObs[0], currObs[1], currObs[2], currObs[3]
    # init galaxy array
    curr_galaxy_mag = currLens[currObs[1]+'_SDSS_lens'][0]
    galaxy_flux = np.power(2.5, desc.slrealizer.return_zeropoint()-curr_galaxy_mag)
    #Only one of sigma, fwhm, and half_light_radius may be specified for Gaussian, so we just add the values!
    convolved_galaxy_size=currLens['REFF'][0]+PSF_HWHM
    galaxy = galsim.Gaussian(half_light_radius=convolved_galaxy_size, flux=galaxy_flux)
    galaxy = galaxy.shear(e1=currLens['ELLIP'][0])
    phi_angle = currLens['PHIE'][0] * galsim.degrees
    galaxy = galaxy.rotate(theta=phi_angle)
    
    for i in xrange(currLens['NIMG']):
        filter_quasar = currLens[currObs[1]+'_SDSS_quasar']
        curr_lens_mag = -2.5*np.log10(abs(currLens['MAG'][0][i])) + filter_quasar[0]
        mag_ratio = np.power(2.5, desc.slrealizer.return_zeropoint()-curr_lens_mag)
        lens = galsim.Gaussian(flux=mag_ratio, sigma=PSF_HWHM)
        lens = lens.shift([currLens['XIMG'][0][i],currLens['YIMG'][0][i]])
        galaxy += lens

    img = galaxy.drawImage(scale=0.2)
    plt.imshow(img.array)
    return img

def galsim_generate_data(currLens, currObs):
    img = desc.slrealizer.galsim_plot_all_objects(currLens, currObs)
    MJD, filter, PSF_HWHM, sky_mag = currObs[0], currObs[1], currObs[2], currObs[3]
    RA, RA_err, DEC, DEC_err, first_moment_x_err_calc, first_moment_y_err_calc, I_xx_err_calc, I_yy_err_calc, I_xy_err_calc = 0, 0, 0, 0, 0, 0, 0, 0, 0
    flux_err_calc = math.pow(10, (22.5 - sky_mag)/2.5)/5 # because Fb = 5 \sigma_b
    shape_info = img.FindAdaptiveMom()
    first_moment_x, first_moment_y = (shape_info.moments_centroid.x-(len(img.array)/2.0))*0.2, (shape_info.moments_centroid.y-(len(img.array)/2.0))*0.2
    # total image intensity for best-fit elliptical Gaussian from adaptive moments. Normally, this field is simply equal to the image flux (for objects that follow a Gaussian light distribution, otherwise it is something approximating the flux). However, if the image was drawn using `drawImage(method='sb')` then moments_amp relates to the flux via flux=(moments_amp)*(pixel scale)^2.
    flux = shape_info.moments_amp # I unit tested, and this is right though
    # moment calculation needs to be changed, because 2d array is not returned
    I_xx = shape_info.moments_sigma * 0.2 * 0.2 # unit of pixels is returned, so change units for arcsec squared 
    I_yy, I_xy, lensID = 1, 0, currLens['LENSID'][0]
    sample_array =  [MJD, filter, RA, RA_err, DEC, DEC_err, first_moment_x, first_moment_x_err_calc, first_moment_y, first_moment_y_err_calc, flux, flux_err_calc, I_xx, I_xx_err_calc, I_yy, I_yy_err_calc, I_xy, I_xy_err_calc, PSF_HWHM, sky_mag, lensID]
    return np.array(sample_array)
