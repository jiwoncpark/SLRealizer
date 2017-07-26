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
import moment
import desc.slrealizer
import scipy.stats
from scipy.stats import moment
from scipy import ndimage
import skimage
import skimage.measure
from skimage.measure import moments
from moment import covariance_matrix
# ======================================================================

"""
Given a specific date and the OM10 catalog, this deblends the sources that are on the catalog.
Assumes null deblender where all the sources are assumed to be observed as single objects.
All the sources are assumed to have Gaussian PSFs.
"""

#global variable that controls the size of the plot
x_min = -5.0
x_max = 5.0
y_min = -5.0
y_max = 5.0
distance = 0.01

number_of_rows = int((x_max - x_min)/distance)
number_of_columns = int((y_max - y_min)/distance)

def deblend_test(currObs, currLens, null_deblend = False, debug=False):
    """
    If the user wants to see the plot drawn by plotting.py in the debug mode, this code draws it.
    Otherwise, it acts like a wrapper method -- this just calls blend_all_objects.
    """
    #print('Deblending starts.....')
    if debug:
        pass
        #print('This is the simple plot of the system')
        #plotting.draw_model(currObs, currLens, debug)
        #plt.clf()
    if null_deblend:
        #print('null deblender')
        image = null_deblending(currObs, currLens, debug)
        show_color_map(image)
    else:
        image = plot_all_objects(currObs, currLens, debug)
        blend_all_objects(currObs, currLens, debug, image)

def null_deblending(moment_matrix, image2, debug):
    x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
    pos = np.dstack((x, y))
    print('moment_matrix: ', moment_matrix)
    zeroth_moment = moment_matrix[0][0]
    first_moment_x = moment_matrix[1][0] / zeroth_moment
    first_moment_y = moment_matrix[0][1] / zeroth_moment
    covariance_matrix = [[moment_matrix[2][0], moment_matrix[1][1]], [moment_matrix[1][1], moment_matrix[0][2]]]
    covariance_matrix /= zeroth_moment**2
#    covariance_matrix = covariance_matrix/(zeroth_moment**2)
    print(zeroth_moment, first_moment_x, first_moment_y, covariance_matrix)

    ## TESTING
    zeroth_moment = np.sum(image2)
    first_moment_x = 0
    first_moment_y = 0
    print(zeroth_moment, 'sum of the image values')
    _, _, covariance_matrix = intertial_axis(image2)
    covariance_matrix = covariance_matrix
    #covariance_matrix = [[0.5, 0], [0, 0.5]]
    print(covariance_matrix, 'different covariance matrix')
    rv = scipy.stats.multivariate_normal([first_moment_x,first_moment_y], covariance_matrix, allow_singular=True) #FIX BUG     
    image = [[0]*number_of_rows for _ in range(number_of_columns)]
    image = image + rv.pdf(pos)*zeroth_moment
    fig, ax = plt.subplots()
    ax.imshow(image)
    #plot_bars(xbar, ybar, cov, ax)
    plt.show()
    return image

def plot_all_objects(currObs, currLens, debug):
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
    #PSF_sigma = desc.slrealizer.fwhm_to_sig(PSF_HWHM)
    ####BUG! DECIDE WHICH ONE IS CORRECT
    PSF_sigma = PSF_HWHM
    rv = scipy.stats.multivariate_normal([galaxy_x,galaxy_y], [[PSF_sigma*PSF_sigma, 0], [0, PSF_sigma*PSF_sigma]], allow_singular=True) #FIX BUG
    image = [[0]*number_of_rows for _ in range(number_of_columns)]
    image = image + rv.pdf(pos)*math.pow(2.5, desc.slrealizer.return_zeropoint() - currLens[filterLens])
    print('multiplication factor : ', math.pow(2.5, desc.slrealizer.return_zeropoint() - currLens[filterLens]))
    # iterate for the lens
    for i in xrange(currLens['NIMG']):
        if debug:
            #print ('XIMG, YIMG, MAG: ', currLens['XIMG'][0][i], currLens['YIMG'][0][i], currLens['MAG'][0][i])
            pass
        mag_ratio = math.pow(2.5, desc.slrealizer.return_zeropoint() -currLens['MAG'][0][i])
        print(mag_ratio)
        print('PSF_sigma: ', PSF_sigma)
            #print ('Magnitude ratio is : ', mag_ratio)
        print('currLensX: ', currLens['XIMG'][0][i], 'currLensY: ', currLens['YIMG'][0][i])
        rv = scipy.stats.multivariate_normal([currLens['XIMG'][0][i],currLens['YIMG'][0][i]], [[PSF_sigma*PSF_sigma, 0], [0, PSF_sigma*PSF_sigma]], allow_singular=True) #, [[PSF_HWHM*PSF_HWHM, 0], [0, PSF_HWHM*PSF_HWHM]])
        image = image + rv.pdf(pos)*mag_ratio #scale
        print('multiplication factor : ', mag_ratio)
    return image

def blend_all_objects(currObs, currLens, debug, input_image):
    if debug:
        show_color_map(input_image)
    # detect sources using IRAF routine
    segm = photutils.detect_sources(input_image, 0.1, 8)
    if debug:
        plt.imshow(segm.data, origin='lower', interpolation='nearest')
    daofind = DAOStarFinder(fwhm=currObs[2], threshold=0.1)
    sources = daofind(input_image)
    if debug:
        show_source_position(sources, input_image)
    # print ('these are the objects that I identified: ', sources)

def show_color_map(input_image):
    """
    Only called in debug mode. Draws the color map to represent the 2d array of the mock lensed system.
    """
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap', ['black', 'green', 'yellow'], 256)
    img2 = plt.imshow(input_image, interpolation='nearest', cmap = cmap2, origin='lower', extent=[x_min, x_max, y_min, y_max], aspect = "auto")
    plt.colorbar(img2,cmap=cmap2)
    plt.show()

def show_source_position(sources, input_image):
    """
    Only called in debug mode. Draws the 2d image and mark where the sources are detected.
    """
    positions = (sources['xcentroid'], sources['ycentroid'])
    apertures = CircularAperture(positions, r=4.)
    norm = ImageNormalize(stretch=SqrtStretch())
    plt.imshow(input_image, cmap='Greys', origin='lower', norm=norm)
    apertures.plot(color='red', lw=3.0, alpha=1.0)

# my guess is : https://www.google.de/search?q=2d+gaussian+area&start=10&sa=N&tbm=isch&imgil=2DtvXl3-AibpvM%253A%253Bkgo2jfIP68y8RM%253Bhttps%25253A%25252F%25252Fstackoverflow.com%25252Fquestions%25252F13658799%25252Fplot-a-grid-of-gaussians-with-matlab&source=iu&pf=m&fir=2DtvXl3-AibpvM%253A%252Ckgo2jfIP68y8RM%252C_&usg=__RmeCJMcLu03ro5YjgKk9fGZ53U8%3D&biw=1183&bih=588&ved=0ahUKEwj6wur9yZTVAhXrhlQKHT54DXo4ChDKNwgz&ei=UuluWfrRLOuN0gK-8LXQBw#imgrc=Puj8GXmbAPS6nM: so amplitude is proportional to the flux...?

#https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m


## copying stackexchange
import numpy as np
import matplotlib.pyplot as plt

def please_work(data):
    xbar, ybar, cov = intertial_axis(data)

    fig, ax = plt.subplots()
    ax.imshow(data)
    plot_bars(xbar, ybar, cov, ax)
    plt.show()

def raw_moment(data, iord, jord):
    nrows, ncols = data.shape
    y, x = np.mgrid[:nrows, :ncols]
    data = data * x**iord * y**jord
    return data.sum()

def intertial_axis(data):
    """Calculate the x-mean, y-mean, and cov matrix of an image."""
    data_sum = data.sum()
    m10 = raw_moment(data, 1, 0)
    m01 = raw_moment(data, 0, 1)
    x_bar = m10 / data_sum
    y_bar = m01 / data_sum
    u11 = (raw_moment(data, 1, 1) - x_bar * m01) / data_sum
    u20 = (raw_moment(data, 2, 0) - x_bar * m10) / data_sum
    u02 = (raw_moment(data, 0, 2) - y_bar * m01) / data_sum
    cov = np.array([[u20, u11], [u11, u02]])
    return x_bar, y_bar, cov

def plot_bars(x_bar, y_bar, cov, ax):
    """Plot bars with a length of 2 stddev along the principal axes."""
    def make_lines(eigvals, eigvecs, mean, i):
        """Make lines a length of 2 stddev."""
        std = np.sqrt(eigvals[i])
        vec = 2 * std * eigvecs[:,i] / np.hypot(*eigvecs[:,i])
        x, y = np.vstack((mean-vec, mean, mean+vec)).T
        return x, y
    print("This is the covariance I calculated: ", cov)
    print("This is the eigenvalue I have", np.linalg.eigh(cov))
    mean = np.array([x_bar, y_bar])
    eigvals, eigvecs = np.linalg.eigh(cov)
    ax.plot(*make_lines(eigvals, eigvecs, mean, 0), marker='o', color='white')
    ax.plot(*make_lines(eigvals, eigvecs, mean, -1), marker='o', color='red')
    ax.axis('image')
