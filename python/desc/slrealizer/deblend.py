# ======================================================================
from __future__ import print_function
from scipy.stats import multivariate_normal
import plotting
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
import matplotlib
import scipy
import scipy.optimize as opt
import math
import photutils
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
from photutils import DAOStarFinder

# ======================================================================

"""
Given a specific date and the OM10 catalog, this deblends the sources that are on the catalog.
Assumes null deblender where all the sources are assumed to be observed as single objects.
All the sources are assumed to have Gaussian PSFs.
"""

#global variable that controls the size of the plot
x_min = -10
x_max = 10
y_min = -10
y_max = 10
distance = 0.01

number_of_rows = int((x_max - x_min)/distance)
number_of_columns = int((y_max - y_min)/distance)

def deblend(currObs, currLens, debug, null_deblender):
    """
    If the user wants to see the plot drawn by plotting.py in the debug mode, this code draws it.
    Otherwise, it acts like a wrapper method -- this just calls blend_all_objects.
    """
    print('Deblending starts.....')
    if debug:
        print('This is the simple plot of the system')
        plotting.draw_model(currObs, currLens, debug)
        #plt.clf()
    if null_deblender:
        null_deblending()
    else:
        blend_all_objects(currObs, currLens, debug)

def null_deblending():


def blend_all_objects(currObs, currLens, debug):
    """
    Given a current observation epoch details and a lensed system, this method blends all the light sources, assuming the gaussian PSF. 
    The blending is done by generating a 2d-gaussian and adding the values to the 2d array.
    Then, that 2d array is passed to the deblending method to find the discrete sources.
    """
    #obsHist has MJD Filter FWHM 5sigmag
    filterLens = currObs[1] + '_SDSS_lens'
    lens_mag = currLens[filterLens]
    #galaxy_x, galaxy_y, PSF_HWHM = currLens['XSRC'][0], currLens['YSRC'][0], currObs[2]
    galaxy_x, galaxy_y, PSF_HWHM = 0, 0, currObs[2] # lens centered at 0,0
    if debug:
        print ('galaxy_x, galaxy_y, PSF_HWHM:'), galaxy_x, galaxy_y, PSF_HWHM
    x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
    pos = np.dstack((x, y))
    rv = scipy.stats.multivariate_normal([galaxy_x,galaxy_y], [[PSF_HWHM*PSF_HWHM, 0], [0, PSF_HWHM*PSF_HWHM]])
    global image
    image = [[0]*number_of_rows for _ in range(number_of_columns)]
    image = image + rv.pdf(pos)
    # iterate for the lens
    for i in xrange(currLens['NIMG']):
        if debug:
            print ('XIMG, YIMG, MAG: ', currLens['XIMG'][0][i], currLens['YIMG'][0][i], currLens['MAG'][0][i])
            mag_ratio = math.pow(2.5, currLens[filterLens]-currLens['MAG'][0][i])
            print ('Magnitude ratio is : ', mag_ratio)
        rv = scipy.stats.multivariate_normal([currLens['XIMG'][0][i],currLens['YIMG'][0][i]], [[PSF_HWHM*PSF_HWHM, 0], [0, PSF_HWHM*PSF_HWHM]]) #, [[PSF_HWHM*PSF_HWHM, 0], [0, PSF_HWHM*PSF_HWHM]])
        image = image + rv.pdf(pos)*math.pow(2.5, currLens[filterLens]-currLens['MAG'][0][i]) #scale
    if debug:
        show_color_map()
    # detect sources using IRAF routine
    segm = photutils.detect_sources(image, 0.1, 8)
    if debug:
        plt.imshow(segm.data, origin='lower', interpolation='nearest')
    daofind = DAOStarFinder(fwhm=PSF_HWHM, threshold=0.1)
    sources = daofind(image)
    if debug:
        show_source_position(sources)
    print ('these are the objects that I identified: ', sources)

def show_color_map():
    """
    Only called in debug mode. Draws the color map to represent the 2d array of the mock lensed system.
    """
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap', ['black', 'green', 'yellow'], 256)
    img2 = plt.imshow(image, interpolation='nearest', cmap = cmap2, origin='lower', extent=[x_min, x_max, y_min, y_max], aspect = "auto")
    plt.colorbar(img2,cmap=cmap2)
    plt.show()

def show_source_position(sources):
    """
    Only called in debug mode. Draws the 2d image and mark where the sources are detected.
    """
    positions = (sources['xcentroid'], sources['ycentroid'])
    apertures = CircularAperture(positions, r=4.)
    norm = ImageNormalize(stretch=SqrtStretch())
    plt.imshow(image, cmap='Greys', origin='lower', norm=norm)
    apertures.plot(color='red', lw=3.0, alpha=1.0)

# my guess is : https://www.google.de/search?q=2d+gaussian+area&start=10&sa=N&tbm=isch&imgil=2DtvXl3-AibpvM%253A%253Bkgo2jfIP68y8RM%253Bhttps%25253A%25252F%25252Fstackoverflow.com%25252Fquestions%25252F13658799%25252Fplot-a-grid-of-gaussians-with-matlab&source=iu&pf=m&fir=2DtvXl3-AibpvM%253A%252Ckgo2jfIP68y8RM%252C_&usg=__RmeCJMcLu03ro5YjgKk9fGZ53U8%3D&biw=1183&bih=588&ved=0ahUKEwj6wur9yZTVAhXrhlQKHT54DXo4ChDKNwgz&ei=UuluWfrRLOuN0gK-8LXQBw#imgrc=Puj8GXmbAPS6nM: so amplitude is proportional to the flux...?

#https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
