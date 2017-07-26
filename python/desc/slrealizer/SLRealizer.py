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
    def deblend(self, lensID=None, null_deblend=True, debug=True, show_plot=True, report_distance=True):
        if lensID is None:
            print('No lens system selected for calculating the statistics')
            return
        if null_deblend is False:
            print('Sorry, working deblender is currently not being supported.')
            return
        import random
        # Keep randomly selecting epochs until we get one that is not in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
            # Now visualize the lens system at the epoch defined by the randomIndex: 
        image2 = desc.slrealizer.plot_all_objects(self.observation[randomIndex], self.catalog.get_lens(lensID), debug)
        moment_matrix = skimage.measure.moments(image2)
        print(skimage.measure.moments(image2),'hello')
#        desc.slrealizer.please_work(image2)
        if show_plot:
            print('#####################BEFORE DEBLEND PLOT LENSES##################################')
#            desc.slrealizer.show_color_map(image2)
        if null_deblend:
            image1 = desc.slrealizer.null_deblending(moment_matrix, debug)
        if show_plot:
            print('#####################PRINTING NULL DEBLENDER\'S PLOT###############################')
            desc.slrealizer.show_color_map(image1)
        if report_distance:
            print('###############################################################################')
            print('Chi squared distance is : ', desc.slrealizer.chi_square_distance(image1, image2))
            print('KL distance is : ', desc.slrealizer.KL_distance(image1, image2))
