# ======================================================================
import os, unittest
import numpy as np
import scipy
import skimage
import desc.slrealizer
import pandas
import galsim
import matplotlib.pyplot as plt
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
        # setup 2d array with given dimensions (declared as global variables)
        global x_min, x_max, y_min, y_max, distance
        x_min, x_max, y_min, y_max = desc.slrealizer.get_x_min(), desc.slrealizer.get_x_max(), desc.slrealizer.get_y_min(), desc.slrealizer.get_y_max()
        distance = desc.slrealizer.get_distance()
        x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
        pos = np.dstack((x, y))
        number_of_rows = int((x_max-x_min)/distance)
        number_of_columns = int((y_max-y_min)/distance)
        image = [[0]*number_of_rows for _ in range(number_of_columns)]
        # return initialized 2d array
        return image

    def generate_gaussian(self):
        image = self.setUp()
        # declare the center of moment and covariance matrix as a global variable, so that test_multivariate can access them
        global x_com, y_com
        global covariance_matrix
        x_com, y_com = 2.0, 3.0
        covariance_matrix = [[0.8, 0.4], [0.4, 0.6]]
        center_of_mass = [x_com, y_com]
        rv = scipy.stats.multivariate_normal([x_com, y_com], covariance_matrix)
        x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
        pos = np.dstack((x, y))
        global flux
        flux = 1
        # regardless of what flux is, we know that all x_com, y_com, and covariance_matrix should not change
        image = image + rv.pdf(pos)*flux
        return image

    def test_multivariate(self):
        # this method tests whether we can get the same moment values from the array produced in generate_gaussian method
        image = self.generate_gaussian()
        moment_matrix = skimage.measure.moments(image)
        zeroth_moment = moment_matrix[0][0]
        print('Zeroth moment: ', zeroth_moment)
        first_moment_x = moment_matrix[0][1] / zeroth_moment
        first_moment_y = moment_matrix[1][0] / zeroth_moment
        moment_matrix = skimage.measure.moments_central(image, first_moment_x, first_moment_y)
        covariance_matric = [[moment_matrix[2][0], moment_matrix[1][1]], [moment_matrix[1][1], moment_matrix[0][2]]]
        covariance_matric /= zeroth_moment
        covariance_matric /= 1/(desc.slrealizer.get_distance()*desc.slrealizer.get_distance())
        print('Covariance matrix is: ', covariance_matric)
        self.assertAlmostEqual(desc.slrealizer.get_x_min()+distance*first_moment_x, x_com, places=1)
        self.assertAlmostEqual(desc.slrealizer.get_y_min()+distance*first_moment_y, y_com, places=1)
        self.assertAlmostEqual(covariance_matric[1][1], covariance_matrix[0][0], places=1)
        self.assertAlmostEqual(covariance_matric[1][0], covariance_matrix[1][0], places=1)
        self.assertAlmostEqual(covariance_matric[0][0], covariance_matrix[1][1], places=1)
        self.assertAlmostEqual(covariance_matric[1][0], covariance_matrix[0][1], places=1)
        return

    def test_null_deblend(self):
        image = self.generate_gaussian()
        flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend(image)
        returned_image=desc.slrealizer.null_deblend_plot(flux, first_moment_x, first_moment_y, covariance_matrix)
        self.assertEqual(image.all(), returned_image.all())

    def test_galsim(self):
        galaxy = galsim.Gaussian(flux=10.0, sigma=10)
        galaxy = galaxy.shift(1, 2)
        galaxy = galaxy.shear(e1=0.2)
        #psf = galsim.Gaussian(flux=1, sigma=3.0)
        #galaxy = galsim.Convolve([galaxy,psf]) 
        img = galaxy.drawImage(scale=.2)
        shape_info = img.FindAdaptiveMom()
        first_moment_x, first_moment_y = (shape_info.moments_centroid.x-(len(img.array)/2.0))*0.2, (shape_info.moments_centroid.y-(len(img.array)/2.0))*0.2
        # total image intensity for best-fit elliptical Gaussian from adaptive moments. Normally, this field is simply equal to the image flux (for objects that follow a Gaussian light distribution, otherwise it is something approximating the flux). However, if the image was drawn using `drawImage(method='sb')` then moments_amp relates to the flux via flux=(moments_amp)*(pixel scale)^2.
        flux = shape_info.moments_amp
        # moment calculation needs to be changed, because 2d array is not returned                       
        I_xx = shape_info.moments_sigma # unit of pixels is returned, so change units for arcsec squared 
        print('I_xx', I_xx)
        print(shape_info.observed_shape.e1)
        print(shape_info.observed_shape.e2)
        print(shape_info.observed_shape.g1, shape_info.observed_shape.g2, shape_info.moments_n_iter)
        
        print('******************HEY****************')
        self.assertAlmostEqual(flux, 10)
        #self.assertAlmostEqual(first_moment_x, 1, places=0)
        #self.assertAlmostEqual(first_moment_y, 2, places=0)
        self.assertAlmostEqual(I_xx*0.2, 3.0, places=2)

# ======================================================================

if __name__ == '__main__':

    unittest.main()

# ======================================================================
