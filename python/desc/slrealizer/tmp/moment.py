from __future__ import print_function
import math
import numpy as np

magnitude_zeropoint = -10 # where we define the flux to be 1

def calculate_lens_zeroth_moment(currObs, currLens):
    """
       Given a lensed system, calculate the zeroth moment(total flux) of lensing galaxy.
    """
    filterLens = currObs[1] + '_SDSS_lens'
    lens_mag = currLens[filterLens]
    lensFlux = pow(2.5, magnitude_zeropoint-lens_mag)
    return lensFlux

def calculate_image_zeroth_moment(currObs, currLens):
    """
        Given a lensed system, calculate the zeroth moment(total flux) of quasar's images.
    """
    print('calculating images zeroth moment')
    # Initialize the flux
    imageFlux=0
    for i in range(int(currLens['NIMG'][0])):
#        print(currLens['MAG'])
        imageMag = currLens['MAG'][0][i]
        print(pow(2.5, magnitude_zeropoint-imageMag))
        imageFlux = imageFlux+ pow(2.5, magnitude_zeropoint-imageMag)
    return imageFlux

def return_image_first_moment(currObs, currLens):
    imageXMoment = []
    imageYMoment = []
    imageFlux = 0
    #filterQuasar = currObs[1] + '_SDSS_quasar'                                                                
    for i in xrange(currLens['NIMG']):
        #if int(currLens['MAG'][0][i]) is 0:                                                                   
        #    imageMag = currLens[filterQuasar]                                                                
        #else:                                                                                               
        imageMag= currLens['MAG'][0][i]
        imageFlux= imageFlux + pow(2.5, magnitude_zeropoint-imageMag)
        # position * flux = moment                                                                             
        imageCurrXMoment = currLens['XIMG'][0][i] * imageFlux
        imageCurrYMoment = currLens['YIMG'][0][i] * imageFlux
        imageXMoment.append(imageCurrXMoment)
        imageYMoment.append(imageCurrYMoment)
    return imageXMoment, imageYMoment

def calculate_image_first_moment(currObs, currLens):
    """
        Given a lensed system, calculate the first moment(total moment for x and y) of quasar's images.
    """
    imageXMoment, imageYMoment = return_image_first_moment(currObs, currLens)
    return sum(imageXMoment), sum(imageYMoment)

def calculate_image_second_moment(currObs, currLens):
    """                                                                                                        
    Given a lensed system, calculate the first moment(total moment for x and y) of quasar's images.          
    """
    imageXMoment = 0
    imageYMoment = 0
    imageFlux = 0
    for i in xrange(currLens['NIMG']):
        imageMag= currLens['MAG'][0][i]
        imageFlux=imageFlux+ pow(2.5, magnitude_zeropoint-imageMag)
        imageCurrXMoment = currLens['XIMG'][0][i] * currLens['XIMG'][0][i] * imageFlux
        imageCurrYMoment = currLens['YIMG'][0][i] * currLens['YIMG'][0][i] * imageFlux
        imageXMoment += imageCurrXMoment
        imageYMoment += imageCurrYMoment
    return imageXMoment, imageYMoment

def calculate_image_second_moment_nc(currObs, currLens):
    """
        Given a lensed system, calculate the first moment(total moment for x and y) of quasar's images.
        """
    imageXMoment = 0
    imageYMoment = 0
    PSF_HWHM = currObs[2]
    for i in xrange(currLens['NIMG']):
        imageXMoment += currLens['XIMG'][0][i] * currLens['XIMG'][0][i] * PSF_HWHM
        imageYMoment += currLens['YIMG'][0][i] * currLens['YIMG'][0][i] * PSF_HWHM
    return imageXMoment, imageYMoment

def zeroth_moment(currObs, currLens):
    """
        Given a lens system, returns the total flux (zeroth moment).
        This function calls the other method named 'calculate_image_zeroth_moment' to calculate the total flux of images. This also calls 'calculate_lens_zeroth_moment' to calculate the total flux of lens.
    """
    return calculate_lens_zeroth_moment(currObs, currLens) + calculate_image_zeroth_moment(currObs, currLens)

def first_moment(currObs, currLens):
    """
        Given a lens system, returns the center of flux (first moment) for each x and y axis
    """
    filterLens = currObs[1] + '_SDSS_lens'
    filterQuasar = currObs[1] + '_SDSS_quasar'
    # Initialize the moment
    yMoment = 0
    xMoment = 0

    # We calculate the contributions to the first moment by the lensing galaxy, and add it to each x and y moments.
    xMoment += calculate_lens_zeroth_moment(currObs, currLens) * 0 # because the galaxy is always centered
    yMoment += calculate_lens_zeroth_moment(currObs, currLens) * 0 # because the galaxy is always centered
    # Add moments calculated through `calculate_image_first_moment` method
    yMoment += calculate_image_first_moment(currObs, currLens)[1]
    xMoment += calculate_image_first_moment(currObs, currLens)[0]

    return xMoment, yMoment


def second_moment(currObs, currLens, convolved=False):
    """
        Given a lens system, returns the second moment for each x and y axis
    """

    filterLens = currObs[1] + '_SDSS_lens'
    filterQuasar = currObs[1] + '_SDSS_quasar'
    # Initialize the moment                                                                                      
    yMoment = 0
    xMoment = 0

    """
    # First approach, using the real definition of the second moment
    # We calculate the contributions to the first moment by the lensing galaxy, and add it to each x and y momen\ts.
    xMoment += calculate_lens_zeroth_moment(currObs, currLens) * currLens['XSRC'][0] * currLens['XSRC'][0]
    yMoment += calculate_lens_zeroth_moment(currObs, currLens) * currLens['YSRC'][0] * currLens['YSRC'][0]
    # Add moments calculated through `calculate_image_first_moment` method                                       
    yMoment += calculate_image_second_moment(currObs, currLens)[1]
    xMoment += calculate_image_second_moment(currObs, currLens)[0]
    """
    
    # Second approach, non-convolved galaxy with quasars using FWHM
    PSF_HWHM = currObs[2]
    xMoment += PSF_HWHM * currLens['XSRC'][0] * currLens['XSRC'][0]
    yMoment += PSF_HWHM * currLens['YSRC'][0] * currLens['YSRC'][0]
    yMoment += calculate_image_second_moment_nc(currObs, currLens)[1]
    xMoment += calculate_image_second_moment_nc(currObs, currLens)[0]
    
    return xMoment, yMoment    

#return xSecondMoment, ySecondMoment
"""
def covariance_matrix(currObs, currLens, convolved=False):
    x_minus_mu_array = []
    y_minus_mu_array = []
    flux = zeroth_moment(currObs, currLens)
    x_first_moment, y_first_moment = first_moment(currObs, currLens)
    x_mu, y_mu = x_first_moment/flux, y_first_moment/flux
    lens_flux = calculate_lens_zeroth_moment(currObs, currLens)
    x_minus_mu_array.append(-x_mu*lens_flux)
    y_minus_mu_array.append(-y_mu*lens_flux)
    for i in xrange(currLens['NIMG']):
        imageMag= currLens['MAG'][0][i]
        imageFlux= pow(2.5, magnitude_zeropoint-imageMag)
        x_minus_mu_array.append((currLens['XIMG'][0][i]-x_mu) * imageFlux)
        y_minus_mu_array.append((currLens['YIMG'][0][i]-y_mu) * imageFlux)
    data_minus_mu = np.squeeze(np.asarray([x_minus_mu_array, y_minus_mu_array]))
    print(data_minus_mu)
#    print(data_minus_mu.shape)
    covariance_matrices = np.matmul(np.transpose(data_minus_mu), data_minus_mu)
    print(covariance_matrices)
    return covariance_matrices
"""

def covariance_matrix(currObs, currLens, convolved=False):
    flux = zeroth_moment(currObs, currLens)
    x_first_moment, y_first_moment = first_moment(currObs, currLens)
    x_mu, y_mu = x_first_moment/flux, y_first_moment/flux
    lens_flux = calculate_lens_zeroth_moment(currObs, currLens)
    lens_array = [-x_mu, -y_mu]
    data_minus_mu = np.matmul(np.asarray(lens_array), np.transpose(np.asarray(lens_array))) * lens_flux
    #data_minus_mu.append([-x_mu*lens_flux, -y_mu*lens_flux])
    #x_minus_mu_array.append(-x_mu*lens_flux)
    #y_minus_mu_array.append(-y_mu*lens_flux)
    for i in xrange(currLens['NIMG']):
        imageMag= currLens['MAG'][0][i]
        imageFlux= pow(2.5, magnitude_zeropoint-imageMag)
        image_array = [(currLens['XIMG'][0][i]-x_mu), (currLens['YIMG'][0][i]-y_mu)]
        data_minus_mu = data_minus_mu + (np.matmul(np.asarray(image_array), np.transpose(np.asarray(image_array))) * imageFlux)

    #print(data_minus_mu)
    #    print(data_minus_mu.shape)
    covariance_matrices = data_minus_mu / flux
    print(covariance_matrices)
    return covariance_matrices

