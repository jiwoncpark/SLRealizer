#hello

import galsim
import numpy as np
import matplotlib.pyplot as plt
import desc.slrealizer
import math

def make_catalog(currLens, currObs):
    image = plot_all_objects(currLens, currObs)

def plot_all_objects(currLens, currObs):
    print('currLens: ', currLens)
    print('currObs: ', currObs)
    MJD, filter, PSF_HWHM, sky_mag = currObs[0], currObs[1], currObs[2], currObs[3]
    filter_AB_offset = 0.00
    curr_galaxy_mag = currLens[currObs[1]+'_SDSS_lens'][0]
    #print 'galaxy magnitude: ', curr_galaxy_mag
    galaxy_flux = np.power(2.5, desc.slrealizer.return_zeropoint()-curr_galaxy_mag+filter_AB_offset)
    #Only one of sigma, fwhm, and half_light_radius may be specified for Gaussian, so we just add the values!
    galaxy = galsim.Gaussian(half_light_radius=currLens['REFF'][0], flux=galaxy_flux)
    #print 'galaxy flux:', galaxy_flux
    for i in xrange(currLens['NIMG']):
        filter_quasar = currLens[currObs[1]+'_SDSS_quasar']
        #print 'curr_quasar_mag: ', filter_quasar
        curr_lens_mag = -2.5*np.log10(abs(currLens['MAG'][0][i]) + filter_AB_offset) + filter_quasar[0] + filter_AB_offset
        mag_ratio = np.power(2.5, desc.slrealizer.return_zeropoint()-curr_lens_mag)
        lens = galsim.Gaussian(flux=mag_ratio, sigma=0.0) # we are going to convolve later, so let's say bc of atmosphere sigma for lens is 0.01"
        lens = lens.shift([currLens['XIMG'][0][i],currLens['YIMG'][0][i]])
        galaxy += lens
        #print('quasar')
        #print('X: ', currLens['XIMG'][0][i], 'Y: ', currLens['YIMG'][0][i], 'mag_ratio', mag_ratio, 'sigma', PSF_HWHM)
    psf = galsim.Gaussian(flux=1, sigma=PSF_HWHM)
    galaxy = galsim.Convolve(galaxy, psf)
    phi_angle = currLens['PHIE'][0] * galsim.degrees
    galaxy = galaxy.shear(e=currLens['ELLIP'][0], beta=phi_angle)
    img = galaxy.drawImage(scale=0.2)
    #return img
#img_copy = img.copy()
    #print(img.FindAdaptiveMom())
    #plt.figure()
    #plt.imshow(img.array)
    #rng = galsim.BaseDeviate(1)
    #noise = galsim.PoissonNoise(rng, sky_level=poisson_noise)
    #noise = galsim.GaussianNoise(sky_level=poisson_noise)
    #img_copy.addNoise(noise)
    #print('********************')
    #print(img.FindAdaptiveMom())
    #plt.figure()
    #plt.imshow(img_copy.array)
    #print('********************')
    #print(img.FindAdaptiveMom())               
#plt.figure()
    #plt.imshow(img.array)
    return img

def generate_data(currLens, currObs, manual_error=True):
    img = desc.slrealizer.plot_all_objects(currLens, currObs)
    MJD, filter, PSF_HWHM, sky_mag = currObs[0], currObs[1], currObs[2], currObs[3]
    RA, RA_err, DEC, DEC_err, first_moment_x_err_calc, first_moment_y_err_calc, size_err = 0, 0, 0, 0, 0, 0, 0
    flux_err_calc = math.pow(10, (22.5 - sky_mag)/2.5)/5 # because Fb = 5 \sigma_b
    try:
        #print('HERE')
        shape_info = img.FindAdaptiveMom()
    except:
        #print('I can\'t manage the expection')
        return
#shape_info = img.EstimateShear()
    #print(shape_info.corrected_shape_err)
    first_moment_x, first_moment_y = (shape_info.moments_centroid.x-(len(img.array)/2.0))*0.2, (shape_info.moments_centroid.y-(len(img.array)/2.0))*0.2
    # total image intensity for best-fit elliptical Gaussian from adaptive moments. Normally, this field is simply equal to the image flux (for objects that follow a Gaussian light distribution, otherwise it is something approximating the flux). However, if the image was drawn using `drawImage(method='sb')` then moments_amp relates to the flux via flux=(moments_amp)*(pixel scale)^2.
    flux = shape_info.moments_amp # I unit tested, and this is right though
    # moment calculation needs to be changed, because 2d array is not returned
    size = shape_info.moments_sigma * 0.2 # unit of pixels is returned, so change units for arcsec
    lensID = currLens['LENSID'][0]
    e1 = shape_info.observed_shape.e1
    e2 = shape_info.observed_shape.e2
    e_squared = (e1*e1+e2*e2)
    e = np.sqrt(e_squared)
    size_err = 0
    if manual_error:
        size += size * noissify_data(desc.slrealizer.get_second_moment_err(), desc.slrealizer.get_second_moment_err_std())
        first_moment_x += noissify_data(desc.slrealizer.get_first_moment_err(), desc.slrealizer.get_first_moment_err_std()) * first_moment_x
        first_moment_y += noissify_data(desc.slrealizer.get_first_moment_err(), desc.slrealizer.get_first_moment_err_std()) * first_moment_y
        flux += flux*noissify_data(desc.slrealizer.get_flux_err(), desc.slrealizer.get_flux_err_std())
        sample_array =  [MJD, filter, RA, RA_err, DEC, DEC_err, first_moment_x, first_moment_x_err_calc, first_moment_y, first_moment_y_err_calc, flux, flux_err_calc, size, size_err, e, PSF_HWHM, sky_mag, lensID]
    print('returned')
    return np.array(sample_array)

def noissify_data(mean, stdev):
    """
    Given a mean and a standard deviation of a measurement, add a noise to the data
    """
    return np.random.normal(loc=mean, scale=stdev)
