from __future__ import print_function
import math
import numpy as np
import desc.slrealizer
import scipy
import scipy.ndimage
from scipy.ndimage import center_of_mass
import scipy.stats
#desc.slrealizer.return_zeropoint = -10 # where we define the flux to be 1
#desc.slrealizer.return_zeropoint = desc.slrealizer.return_zeropoint()

def zeroth_moment(image):
    return image.sum()
    
def first_moment(image):
    return scipy.ndimage.center_of_mass(image)

def second_moment(image):
    return scipy.stats.moment(image, moment=2, axis=None)

def return_zeropoint():
    return 0
