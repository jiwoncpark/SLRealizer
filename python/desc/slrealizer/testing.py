import desc.slrealizer
import deblend
import scipy
import scipy.stats

def chi_square_distance(image1, image2):
    return scipy.stats.chisquare(image1.ravel(), image2.ravel(), axis=None)


def KL_distance(image1, image2):
    return scipy.stats.entropy(image1.ravel(), image2.ravel())
