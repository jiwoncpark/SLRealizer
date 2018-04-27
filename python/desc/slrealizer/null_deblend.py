import galsim
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from utils import *

def draw_system(currLens, currObs, save_dir=None):
    '''
    Draws all objects of the given lens system
    in the given observation conditions using GalSim,
    and returns the tuple (Matplotlib image, GalSim object).
    
    Keyword arguments:
    currLens -- a row of the OM10 DB
    currObs -- a row of the observation history df
    save_dir -- directory in which to save the image
    '''
    MJD, band, PSF_FWHM, sky_mag = currObs['expMJD'], currObs['filter'], currObs['FWHMeff'], currObs['fiveSigmaDepth']
    ###############
    # Lens galaxy #
    ###############
    # Compute flux ratio
    gal_mag = currLens[band + '_SDSS_lens'][0]
    galaxy_flux = from_mag_to_flux(gal_mag\
                                      - return_zeropoint()\
                                      - get_filter_AB_offset())
    # Define GalSim object, give right ellipticity
    galaxy = galsim.Gaussian(half_light_radius=currLens['REFF_T'][0], flux=galaxy_flux)\
                   .shear(e=currLens['ELLIP'][0], beta=currLens['PHIE'][0] * galsim.degrees)
    #################
    # Lensed quasar #
    #################
    for i in xrange(currLens['NIMG']):
        # TODO ask how to factor in magnification?
        # Compute flux ratio
        curr_obj_mag = currLens[band + '_SDSS_quasar'][0]\
                       + from_flux_to_mag(currLens['MAG'][0][i] + get_filter_AB_offset())\
                       + get_filter_AB_offset()
        obj_flux = from_mag_to_flux(curr_obj_mag - return_zeropoint()) # don't subtract AB offset?
        # Define GalSim object, shift accordingly
        lens = galsim.Gaussian(flux=obj_flux, sigma=0.0)\
                     .shift([currLens['XIMG'][0][i], currLens['YIMG'][0][i]])
        galaxy += lens
    big_fft_params = galsim.GSParams(maximum_fft_size=10240)
    psf = galsim.Gaussian(flux=1.0, fwhm=PSF_FWHM)
    galsim_obj = galsim.Convolve(galaxy, psf, gsparams=big_fft_params)
    # TODO make pixel scale globally configurable
    # TODO since precision of HSM depends on pixel resolution, either don't pass in img to HSM OR manually define good resolution
    img = galsim_obj.drawImage(scale=get_pixel_scale())
    if save_dir is not None:
        plt.imshow(img.array, interpolation='none', aspect='auto')
        plt.savefig(save_dir+'before_deblend.png')
        plt.close()
    return img, galsim_obj

def generate_row(currLens, currObs, manual_error=True):
    '''
    Returns a dictionary of lens system's properties
    computed from one lens system and one observation,
    which makes up a row of the source table.
    
    Keyword arguments:
    currLens -- a row of the OM10 DB
    currObs -- a row of the observation history df
    manual_error -- if True, adds some predefined noise (default: True)
    '''
    lensID = currLens['LENSID'][0]
    img, obj = draw_system(currLens, currObs)
    histID, MJD, band, PSF_FWHM, sky_mag = currObs
    RA, RA_err, DEC, DEC_err, first_moment_x_err_calc, first_moment_y_err_calc, size_err = 0, 0, 0, 0, 0, 0, 0
    flux_err_calc = from_mag_to_flux(sky_mag-22.5)/5.0 # because Fb = 5 \sigma_b
    try:
        shape_info = img.FindAdaptiveMom()
    except:
        print "HSM failed"
        return None
    nx, ny = img.array.shape
    # Calculate the real position from the arbitrary pixel position
    first_moment_x, first_moment_y = (shape_info.moments_centroid.x - 0.5*nx)*get_pixel_scale(),\
                                     (shape_info.moments_centroid.y - 0.5*ny)*get_pixel_scale()
    # TODO just use known flux?
    flux = shape_info.moments_amp
    size = shape_info.moments_sigma * get_pixel_scale() # unit of pixels is returned, so change units for arcsec
    e1 = shape_info.observed_shape.e1
    e2 = shape_info.observed_shape.e2
    size_err = 0
    if manual_error:
        size += noissify_data(get_second_moment_err(), get_second_moment_err_std())
        first_moment_x += noissify_data(get_first_moment_err(), get_first_moment_err_std(), first_moment_x)
        first_moment_y += noissify_data(get_first_moment_err(), get_first_moment_err_std(), first_moment_y) 
        flux += noissify_data(get_flux_err(), get_flux_err_std(), flux)
    # TODO consider just returning shear object, etc.
    params_dict = {'MJD': MJD, 'filter': band, 'RA': RA, 'RA_err': RA_err, 'DEC': DEC, 'DEC_err': DEC_err,
                   'x': first_moment_x, 'x_com_err': first_moment_x_err_calc,
                   'y': first_moment_y, 'y_com_err': first_moment_y_err_calc,
                   'flux': flux, 'flux_err': flux_err_calc, 'size': size, 'size_err': size_err, 'e1': e1, 'e2': e2,
                   'psf_sigma': PSF_HWHM, 'sky': sky_mag, 'lensid': lensID}
    return params_dict

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