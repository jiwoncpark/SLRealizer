import desc.slrealizer
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math
import skimage

class SLRealizer(object):

    # Maybe make a separate method for it
    # Best OOP practice?

    def __init__(self, catalog=None, observation="../../data/twinkles_observation_history.csv"):
        """
            Reads in a lens sample catalog and observation data.
            We assume lenses are OM10 lenses and observation file is .csv file
        """
        self.catalog = catalog
        self.observation = pd.read_csv(observation,index_col=0).as_matrix()

    def plot_lens_random_date(self, lensID=None, debug=False, convolve=False):
        """
        Given a specific lens, this code plots a lens after choosing a random observation epoch.
        """

        if lensID is None:
            print 'No lens system selected for plotting.'
            return
        import random
        # Keep randomly selecting epochs until we get one that is not in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        # Now visualize the lens system at the epoch defined by the randomIndex:
        desc.slrealizer.draw_model(self.observation[randomIndex],
                                   self.catalog.get_lens(lensID),
                                   convolve, debug)
        return

    # For now set all to true so that we can debug easily
    def deblend(self, lensID=None, null_deblend=True, debug=False, show_plot=True, version=None, report_distance=True):
        if lensID is None:
            print('No lens system selected for calculating the statistics')
            return
        if null_deblend is False:
            print('Sorry, working deblender is currently not being supported.')
            return
        if version is None:
            print('Select either 1 or 2')
            return
        import random
        # Keep randomly selecting epochs until we get one that is not in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        image2 = desc.slrealizer.plot_all_objects(self.observation[randomIndex], self.catalog.get_lens(lensID), debug)
        print('#####################BEFORE DEBLEND PLOT LENSES##################################')
        desc.slrealizer.show_color_map(image2)
        if version is 1:
            desc.slrealizer.please_work(image2)
            image = desc.slrealizer.null_deblend_v1(image2)
        if version is 2:
            image = desc.slrealizer.null_deblend_v2(image2)
        if version is 3:
            image = desc.slrealizer.null_deblend_v3(image2)
        print('#####################PRINTING NULL DEBLENDER\'S PLOT###############################')
        desc.slrealizer.show_color_map(image)
        if report_distance:
            print('###############################################################################')
            print('Chi squared distance is : ', desc.slrealizer.chi_square_distance(image, image2))
            print('KL distance is : ', desc.slrealizer.KL_distance(image, image2))
