#import numpy
from scipy.stats import multivariate_normal
import plotting
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
import matplotlib

#from matplotlib import mpl

"""

deblend.py deblends the light sources for OM10. This assumes null deblender where all the sources are observed as a single object. All the sources are assumed to have Gaussian PSFs.
"""

number_of_rows = 10000
number_of_columns = 10000
image = [[0]*number_of_rows for _ in range(number_of_columns)]

def deblend(currObs, currLens):
    print 'Deblending starts.....'
    print 'The system looks like this....'
    plotting.draw_model(currObs, currLens)
    blend_all_objects(currObs, currLens)

def blend_all_objects(currObs, currLens):
    fig_d = plt.figure()
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap', ['black', 'blue', 'green', 'yellow'], 256)
    x, y = np.mgrid[-3:3:.01, -3:3:.01]
    
def gaussian_fit():


#https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m

