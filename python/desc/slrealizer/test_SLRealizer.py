# ======================================================================

import os, unittest
import desc.slrealizer
import null_deblend_v1
import null_deblend_v2
import null_deblend_v3 # doesn't work
import numpy as np
import scipy
import skimage

# ======================================================================

class TestCase(unittest.TestCase):
    
    """
    Tests SLRealizer class.
    Notes
    -----
    Execute these tests with:
        nosetests
    from anywhere in the module, provided you have run
        pip install nose
    """
    
    # ------------------------------------------------------------------

    def setUp(self):
        global x_min, x_max, y_min, y_max, distance
        x_min, x_max, y_min, y_max = 0, 10, 0, 10
        distance = 0.01
        x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
        pos = np.dstack((x, y))
        number_of_rows = int((x_max-x_min)/distance)
        number_of_columns = int((y_max-y_min)/distance)
        image = [[0]*number_of_rows for _ in range(number_of_columns)]
        return image

    def generate_gaussian(self):
        image = self.setUp()
        global x_com, y_com
        global covariance_matrix
        x_com, y_com = 5, 5
        covariance_matrix = [[0.1, 0.2], [0.2, 0.8]]
        center_of_mass = [x_com, y_com]
        rv = scipy.stats.multivariate_normal([5, 5], covariance_matrix)
        x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
        pos = np.dstack((x, y))
        image = image + rv.pdf(pos)
        return image

    def test_multivariate_v1(self):
        image = self.generate_gaussian()
        first_moment_x, first_moment_y, covariance_matric = desc.slrealizer.inertial_axis(image)
        print('First Moment for multivariate v1 ', first_moment_x, first_moment_y)
        print('Real valeud should be', x_com, y_com)
        self.assertAlmostEqual(distance*first_moment_x, x_com, places=5)
        self.assertAlmostEqual(distance*first_moment_y, y_com, places=5)
        #print('eigen vector of covariance matrix is: ', np.eig(cov))
        self.assertAlmostEqual(covariance_matric[0][0], covariance_matrix[0][0], places=5)
        return

    def test_multivariate_v2(self):
        image = self.generate_gaussian()
        moment_matrix = skimage.measure.moments(image)
        zeroth_moment = moment_matrix[0][0]
        first_moment_x = moment_matrix[1][0] / zeroth_moment
        first_moment_y = moment_matrix[0][1] / zeroth_moment
        covariance_matric = [[moment_matrix[2][0], moment_matrix[1][1]], [moment_matrix[1][1], moment_matrix[0][2]]]
        self.assertAlmostEqual(distance*first_moment_x, x_com, places=5)
        self.assertAlmostEqual(distance*first_moment_y, y_com, places=5)
        self.assertAlmostEqual(covariance_matric[0][0]-zeroth_moment*first_moment_x*first_moment_x, covariance_matrix[0][0], places=5)
        return

# ======================================================================

if __name__ == '__main__':

    unittest.main()

# ======================================================================
