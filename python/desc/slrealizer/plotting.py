import plotting
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math

def determine_alpha(mag_ratio):
    """
    Bug: No docstring
    Bug: Each image should have its own alpha, because the different
         images have different magnifications (these are in the OM10 "MAG" columns, and need to be applied to each image's sourceX
         magnitude.)
    Bug: function name should be "determine_alpha" to be PEP8 compliant
    """
    if(mag_ratio>1):
        quasar_alpha = 1/mag_ratio
        lens_alpha = 1
    else:
        quasar_alpha = 1
        lens_alpha = mag_ratio
    return quasar_alpha, lens_alpha

def determineScale(lensX, sourceX, lensY, sourceY):
    if (abs(lensX[0]+sourceX))>(abs(sourceX)):
        plotX = abs(lensX[0]+sourceX)
    else:
        plotX = abs(sourceX)
    if (abs(lensY[0]+sourceY))>(abs(sourceY)):
        plotY = abs(lensY[0]+sourceY)
    else:
        plotY = abs(sourceY)
    return plotX, plotY

def draw_model(currObs, currLens, convolve=False):
    #obsHist has MJD Filter FWHM 5sigmag
    filterQuasar = currObs[1] + '_SDSS_quasar'
    filterLens = currObs[1] + '_SDSS_lens'
    
    if (currObs[1]=='u'):
        circleColor = 'violet'
    elif (currObs[1]=='g'):
        circleColor = 'blue'
    elif (currObs[1]=='r'):
        circleColor = 'green'
    elif (currObs[1]=='i'):
        circleColor = 'orange'
    elif (currObs[1]=='z'):
        circleColor = 'red'
    else:
        raise ValueError('Unknown filter name '+currObs[1])
    
    scale_factor = 2
    PSF_HWHM = currObs[2]/scale_factor
                    
    plt.axis('scaled')
    plt.ylim(-3, 3)
    plt.xlim(-3, 3)
    plt.xlabel('xPosition')
    plt.ylabel('yPosition')
    plt.title('Observation with ' + 'filter ' + currObs[1] + ' on MJD ' + str(currObs[0]))

    #Draw lens images:
    if(convolve):
    # For now, we just let lensFWHM to be zero.
    #lensFWHM = currLens['REFF']
        lensFWHM = 0.0
        galaxyHWHM = convolve(currObs[2], lensFWHM)/scale_factor
            # Bug: the LENS position is in the CENTER of the field.
            #      The SOURCE position is UNOBSERVABLE.
            #      The IMAGE positions are in teh OM10 catalog in the
            #      XIMG and YIMG columns.
        source = plt.Circle((0, 0),
                            radius=galaxy_HWHM, alpha=1.0, fc=circleColor, linewidth=0)
    else:
        source = plt.Circle((0, 0),
                              radius=PSF_HWHM,alpha=1.0,fc=circleColor, linewidth=0,edgecolor='b', linestyle='dashed')
    plt.gca().add_patch(source)
    filterLens = currObs[1] + '_SDSS_lens'
    lens_mag = currLens[filterLens]

    # Draw quasar images:
    for i in range(4):
        sourceX = currLens['XIMG'][0][i]
        sourceY = currLens['YIMG'][0][i]
        quasar_mag = currLens['MAG'][0][i]
        print quasar_mag
        print lens_mag
        mag_ratio = math.pow(2.5, -lens_mag+quasar_mag)
        quasar_alpha, lens_alpha = determine_alpha(mag_ratio)
        if determine_alpha < 0.1:
            quasar_alpha = 0.1
        print "In 'draw_model', mag_ratio, quasar_alpha, lens_alpha =", \
            mag_ratio, quasar_alpha, lens_alpha
        image = plt.Circle((sourceX, sourceY),
                            radius=PSF_HWHM,
                            alpha=quasar_alpha,
                            fc=circleColor, linewidth=0)
                            # Draw lens galaxy:
        plt.gca().add_patch(image)
                                            
            #seeing = plt.Circle(((-plotY-2*currObs[2])/scale_factor, (plotX-2*currObs[2])/scale_factor), radius=PSF_HWHM, alpha=0.1, fc='k')
    seeing = plt.Circle((-2.5, -2.5),
                radius=PSF_HWHM,
                alpha=0.1,
                fc='black')
    plt.legend((source, seeing, image), ('QSO images', 'PSF', 'Lens galaxy'), fontsize=10)
    plt.gca().add_patch(seeing)


def convolve(obsFWHM, initialFWHM=0.0):
    seeing = obsFWHM
    seeingsigma = seeingToSig(seeing)
    initSigma = seeingToSig(initialFWHM)
    convolveSigma = seeingsigma + initSigma
    return convolveSigma

def seeingToSig(seeing):
    return seeing / np.sqrt(8 * np.log(2))

def sigmaTofwhm(sigma):
    return sigma * np.sqrt(8 * np.log(2))

#def fail_to_deblend():
#    plot_fail_to_deblend()

#def plot_fail_to_deblend():
#    
