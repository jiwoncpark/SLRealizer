import galsim
import numpy as np
import matplotlib.pyplot as plt
import desc.slrealizer
import math

def make_catalog(currLens, currObs):
    image = plot_all_objects(currLens, currObs)

def plot_all_objects(currLens, currObs, save_dir=None):
    MJD, filter, PSF_HWHM, sky_mag = currObs[0], currObs[1], currObs[2], currObs[3]
    filter_AB_offset = 0.00
    curr_galaxy_mag = currLens[currObs[1]+'_SDSS_lens'][0] # 0-indexing just gets the magnitude value
    galaxy_flux = np.power(2.5, desc.slrealizer.return_zeropoint() - 
                                curr_galaxy_mag + 
                                filter_AB_offset) # Note mag diff ~ 2.5**(flux ratio)
    #Only one of sigma, fwhm, and half_light_radius may be specified for Gaussian, so we just add the values!
    # TODO unit of REFF?
    galaxy = galsim.Gaussian(half_light_radius=currLens['REFF'][0], flux=galaxy_flux)
    big_fft_params = galsim.GSParams(maximum_fft_size=10240)
    phi_angle = currLens['PHIE'][0] * galsim.degrees
    galaxy = galaxy.shear(e=currLens['ELLIP'][0], beta=phi_angle)
    for i in xrange(currLens['NIMG']):
        filter_quasar = currLens[currObs[1]+'_SDSS_quasar']
        curr_lens_mag = -2.5*np.log10(abs(currLens['MAG'][0][i]) + filter_AB_offset) + filter_quasar[0] + filter_AB_offset
        mag_ratio = np.power(2.5, desc.slrealizer.return_zeropoint() - curr_lens_mag)
        lens = galsim.Gaussian(flux=mag_ratio, sigma=0.0) # we are going to convolve later, so lens is a point source for now
        lens = lens.shift([currLens['XIMG'][0][i], currLens['YIMG'][0][i]])
        galaxy += lens
    psf = galsim.Gaussian(flux=1.0, sigma=PSF_HWHM)
    galaxy = galsim.Convolve(galaxy, psf, gsparams=big_fft_params)
    # TODO make pixel scale globally configurable
    # TODO since precision of HSM depends on pixel resolution, either don't pass in img to HSM OR manually define good resolution
    img = galaxy.drawImage(scale=0.2)
    if save_dir is not None:
        plt.imshow(img.array, interpolation='none', aspect='auto')
        plt.savefig(save_dir+'before_deblend.png')
        plt.close()
    # TODO need to close fig
    return img

def generate_data(currLens, currObs, manual_error=True):
    img = desc.slrealizer.plot_all_objects(currLens, currObs)
    MJD, band, PSF_HWHM, sky_mag = currObs[0], currObs[1], currObs[2], currObs[3]
    RA, RA_err, DEC, DEC_err, first_moment_x_err_calc, first_moment_y_err_calc, size_err = 0, 0, 0, 0, 0, 0, 0
    flux_err_calc = np.power(10, (22.5 - sky_mag)/2.5)/5. # because Fb = 5 \sigma_b
    try:
        shape_info = img.FindAdaptiveMom()
    except:
        print "HSM failed"
        return
    nx, ny = img.array.shape
    # TODO make pixel scale globally configurable
    pixel_scale = 0.2 # arcsec/pixel
    # Calculate the real position from the arbitrary pixel position
    first_moment_x, first_moment_y = (shape_info.moments_centroid.x - 0.5*nx)*pixel_scale, \
                                     (shape_info.moments_centroid.y - 0.5*ny)*pixel_scale 
    # Total image intensity for best-fit elliptical Gaussian from adaptive moments
    flux = shape_info.moments_amp
    size = shape_info.moments_sigma * pixel_scale # unit of pixels is returned, so change units for arcsec
    lensID = currLens['LENSID'][0]
    
    e = shape_info.observed_shape.e
    # FIXIT not sure if flooring necessary. ceiling?
    e = max(0.001, e)
    #if e == 0: e = 0.001 # prevent division by zero
    phi = shape_info.observed_shape.beta
    size_err = 0
    if manual_error:
        size += size * noissify_data(desc.slrealizer.get_second_moment_err(), desc.slrealizer.get_second_moment_err_std())
        first_moment_x += noissify_data(desc.slrealizer.get_first_moment_err(), desc.slrealizer.get_first_moment_err_std()) * first_moment_x
        first_moment_y += noissify_data(desc.slrealizer.get_first_moment_err(), desc.slrealizer.get_first_moment_err_std()) * first_moment_y
        flux += flux*noissify_data(desc.slrealizer.get_flux_err(), desc.slrealizer.get_flux_err_std())
    # TODO consider just returning shear object, etc.
    params_dict = {'MJD': MJD, 'band': band, 'RA': RA, 'RA_err': RA_err, 'DEC': DEC, 'DEC_err': DEC_err,
                   'first_moment_x': first_moment_x, 'first_moment_x_err_calc': first_moment_x_err_calc,
                   'first_moment_y': first_moment_y, 'first_moment_y_err_calc': first_moment_y_err_calc,
                   'flux': flux, 'flux_err_calc': flux_err_calc, 'size': size, 'size_err': size_err,
                   'e': e, 'phi': phi, 'PSF_HWHM': PSF_HWHM, 'sky_mag': sky_mag, 'lensID': lensID}
    return params_dict

def noissify_data(mean, stdev):
    """
    Given a mean and a standard deviation of a measurement, add a noise to the data
    """
    return np.random.normal(loc=mean, scale=stdev)

def main():
    #TODO make unit test
    import om10
    
    catalog_f = '../data/qso_mock.fits'
    db = om10.DB(catalog=catalog_f)
    
    lens_example = db.get_lens(db.sample[0]['LENSID'])
    imag = lens_example['i_SDSS_lens'][0]
    print "galaxy magnitude: ", imag
    
if __name__ == "__main__":
    main()