from __future__ import print_function
import math
import numpy as np

def calculate_lens_zeroth_moment(currObs, currLens):
    """
       Given a lensed system, calculate the zeroth moment(total flux) of lensing galaxy.
    """
    filterLens = currObs[1] + '_SDSS_lens'
    lens_mag = currLens[filterLens]
    lensFlux = pow(2.5, lens_mag)
    return lensFlux

def calculate_image_zeroth_moment(currObs, currLens):
    """
        Given a lensed system, calculate the zeroth moment(total flux) of quasar's images.
    """
    # Initialize the flux
    filterQuasar = currObs[1] + '_SDSS_quasar'
    imageFlux = 0
    nImage = 0
    if int(currLens['NIMG'][0]) is 4:
        #print('here')
        nImage = 4
    if int(currLens['NIMG'][0]) is 2:
        #print('here')
        nImage = 2
    #print(nImage)
    for i in range(nImage):
        if int(currLens['MAG'][0][i]) is 0:
            imageMag = currLens[filterQuasar]
        else :
            #print(currLens[filterQuasar])
            #print(abs(currLens['MAG'][0][i]))
            imageMag= currLens[filterQuasar] - 2.5*np.log10(abs(currLens['MAG'][0][i]))
        #print('imageMag')
        #print(imageMag)
        imageFlux += pow(2.5, -imageMag)
        #print('imageFlux')
        #print(imageFlux)
    return imageFlux

def calculate_image_first_moment(currObs, currLens):
    """
        Given a lensed system, calculate the first moment(total moment for x and y) of quasar's images.
    """
    imageXMoment = 0
    imageYMoment = 0
    imageFlux = 0
    filterQuasar = currObs[1] + '_SDSS_quasar'
    for i in range(4):
        if int(currLens['MAG'][0][i]) is 0:
            imageMag = currLens[filterQuasar]
        else:
            imageMag= currLens[filterQuasar] - 2.5*np.log10(abs(currLens['MAG'][0][i]))
        imageFlux += pow(2.5, -imageMag)
        # position * flux = moment
        imageCurrXMoment = currLens['XIMG'][0][i] * imageFlux
        imageCurrYMoment = currLens['YIMG'][0][i] * imageFlux
        imageXMoment += imageCurrXMoment
        imageYMoment += imageCurrYMoment
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
    xMoment += calculate_lens_zeroth_moment(currObs, currLens) * currLens['XSRC'][0]
    yMoment += calculate_lens_zeroth_moment(currObs, currLens) * currLens['YSRC'][0]
    # Add moments calculated through `calculate_image_first_moment` method
    yMoment += calculate_image_first_moment(currObs, currLens)[1]
    xMoment += calculate_image_first_moment(currObs, currLens)[0]

    #print('Image X Moment')
    #print(xMoment)
    #print('Image Y Moment')
    #print(yMoment)

    return xMoment, yMoment


def second_moment(currObs, currLens, convolved=False):
    """
        Given a lens system, returns the second moment for each x and y axis
    """
    
    xMoment, yMoment = first_moment(currObs, currLens)
    
    
    xSecondMoment = xMoment * xMoment
    ySecondMoment = yMoment * yMoment
    
    yTotalSigma = 0
    xTotalSigma = 0
    if currLens['NIMG'][0] is 2:
        yTotalSigma += currObs[2] * 2
        xTotalSigma += currObs[2] * 2
    if currLens ['NIMG'][0] is 4:
        yTotalSigma += currObs[2] * 4
        xTotalSigma += currObs[2] * 4
    if not convolved:
        yTotalSigma += currObs[2]
        xTotalSigma += currObs[2]
    else:
        lens_FWHM = currLens['REFF']
        galaxy_HWHM = desc.slrealizer.convolve(PSF_HWHM, lens_FWHM)/scale_factor
        yTotalSigma += galaxy_HWHM
        xTotalSigma += galaxy_HWHM

    xSecondMoment += xTotalSigma * xTotalSigma
    ySecondMoment += yTotalSigma * yTotalSigma
    return xSecondMoment, ySecondMoment
