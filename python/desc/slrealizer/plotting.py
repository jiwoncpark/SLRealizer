# ====================================================================

import plotting
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math

# ======================================================================

x_min = -10
x_max = 10
y_min = -10
y_max = 10

def draw_model(curr_obs, curr_lens, convolve=False):
    
    """
        given a lensed system and a observation epoch, this visualize how the system would look like.
        """
    
    #obsHist has MJD Filter FWHM 5sigmag
    circle_color = choose_color(curr_obs)
    filter_lens = curr_obs[1] + '_SDSS_lens'
    #scale_factor = 2
    PSF_HWHM = curr_obs[2] #/scale_factor
    
    init_plot()
    
    #Draw lensing galaxy' images:
    if(convolve):
        # For now, we just let lensFWHM to be zero.
        #lensFWHM = curr_lens['REFF']
        lensFWHM = 0.0
        galaxyHWHM = convolve(curr_obs[2], lensFWHM) #/scale_factor
        galaxy = plt.Circle((0, 0),
                            radius=galaxy_HWHM, alpha=1.0, fc='c', linewidth=0)
    else:
        galaxy = plt.Circle((0, 0),
                            radius=PSF_HWHM,alpha=1.0,fc='c', linewidth=0,edgecolor='b', linestyle='dashed')
    plt.gca().add_patch(galaxy)
    lens_mag = curr_lens[filter_lens]

    #sort array to find brightest image
    brightest_image_magnitude = min(filter(lambda a: a != 2, curr_lens['MAG'][0]))
    print 'brightest_image_magnitude : ',
    print brightest_image_magnitude
    
    # Draw quasar images:
    for i in xrange(curr_lens['NIMG']):
        sourceX = curr_lens['XIMG'][0][i]
        sourceY = curr_lens['YIMG'][0][i]
        image_mag = curr_lens['MAG'][0][i]
        print 'quasar image mag:', image_mag, 'FWHM:', PSF_HWHM, 'Quasar X:', sourceX, 'Quasar Y: ', sourceY
        # image alpha no longer in use
        #mag_ratio = math.pow(2.5, -lens_mag+image_mag)
        #image_alpha, lens_alpha = determine_alpha(mag_ratio)
        #if image_alpha < 0.1:
        #    image_alpha = 0.1
        image_alpha = math.pow(2.5, brightest_image_magnitude-image_mag)
        #print "In 'draw_model', mag_ratio, quasar_alpha, lens_alpha =", mag_ratio, image_alpha, lens_alpha
        image = plt.Circle((sourceX, sourceY),
                           radius=PSF_HWHM,
                           alpha=image_alpha,
                           fc=circle_color, linewidth=0)
                           # Draw lens galaxy:
        plt.gca().add_patch(image)

    #seeing = plt.Circle(((-plotY-2*curr_obs[2])/scale_factor, (plotX-2*curr_obs[2])/scale_factor), radius=PSF_HWHM, alpha=0.1, fc='k')
    seeing = plt.Circle((-2.5, -2.5),
                        radius=PSF_HWHM,
                        alpha=0.1,
                        fc='black')
    plt.legend((galaxy, seeing, image), ('Lens Galaxy', 'PSF', 'QSO images'), fontsize=10)
    plt.gca().add_patch(seeing)


def determine_alpha(mag_ratio):
    """
        Given the magnitude ratio (which is given by 2.5^(m1-m2)), calculate the alpha values for each patch that represents a quasar image
        """
    if(mag_ratio>1):
        quasar_alpha = 1/mag_ratio
        lens_alpha = 1
    else:
        quasar_alpha = 1
        lens_alpha = mag_ratio
    return quasar_alpha, lens_alpha

def choose_color(curr_obs):
    """
        Given a filter, this method choose an appropriate color for the patches (which have a shape of circles)
        """
    if (curr_obs[1]=='u'):
        circle_color = 'violet'
    elif (curr_obs[1]=='g'):
        circle_color = 'blue'
    elif (curr_obs[1]=='r'):
        circle_color = 'green'
    elif (curr_obs[1]=='i'):
        circle_color = 'orange'
    elif (curr_obs[1]=='z'):
        circle_color = 'red'
    else:
        raise ValueError('Unknown filter name '+curr_obs[1])
    return circle_color

def init_plot():
    
    """
        initialize the plot by setting axis and labels
        """
    
    plt.axis('scaled')
    plt.ylim(y_min, y_max)
    plt.xlim(x_min, x_max)
    plt.xlabel('xPosition')
    plt.ylabel('yPosition')
    plt.title('Observation with ' + 'filter ' + curr_obs[1] + ' on MJD ' + str(curr_obs[0]))

def convolve(obs_FWHM, initialFWHM=0.0):
    """
        Given observation date's FWHM and the galaxy's size, compute the convolved gaussian's FWHM.
        """
    seeing = obs_FWHM
    seeing_sigma = fwhm_to_sig(seeing)
    init_sigma = fwhm_to_sig(initialFWHM)
    convolve_sigma = seeing_sigma + init_sigma
    return convolve_sigma

def fwhm_to_sig(seeing):
    """
        Given the FWHM, return the sigma value
        """
    return seeing / np.sqrt(8 * np.log(2))

def sigma_to_fwhm(sigma):
    """
        Given the sigma value, return FWHM
        """
    return sigma * np.sqrt(8 * np.log(2))

"""
    Will be eventually deleted, keeping just in case
def determine_scale(lens_X, source_X, lens_Y, source_Y):
    if (abs(lens_X[0]+source_X))>(abs(source_X)):
        plot_X = abs(lens_X[0]+source_X)
    else:
        plot_X = abs(source_X)
    if (abs(lens_Y[0]+source_Y))>(abs(source_Y)):
        plot_Y = abs(lens_Y[0]+source_Y)
    else:
        plot_Y = abs(source_Y)
    return plot_X, plot_Y
"""
