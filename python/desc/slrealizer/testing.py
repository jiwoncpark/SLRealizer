import desc.slrealizer
import deblend
import scipy
import scipy.stats

def chi_square_distance(image1, image2):
    chi2, p, dof, ex = scipy.stats.chi2_contingency(image1, image2)
    print('chi squared: ', chi2, 'p: ', p)
    return chi2, p


def KL_distance(image1, image2):
    return scipy.stats.entropy(image1.ravel(), image2.ravel())
