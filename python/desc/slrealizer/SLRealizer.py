import desc.slrealizer
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math

class SLRealizer(object):

    scale_factor = 2

    # Maybe make a separate method for it
    # Best OOP practice?

    def __init__(self, catalog=None, observation="../../data/twinkles_observation_history.csv"):
        """
            Reads in a lens sample catalog and observation data.
            We assume lenses are OM10 lenses and observation file is .csv file
        """
        self.catalog = catalog
        self.observation = pd.read_csv(observation,index_col=0).as_matrix()

    # def plot_lens_random_date(self, lensID = 7176527, convolve=False):
    #     plotting.plot_lens_random_date(self, lensID, convolve)

    def plot_lens_random_date(self, lensID=None, debug=False, convolve=False):
        if lensID is None:
            print 'No lens system selected for plotting.'
            return
        import random
        # Keep randomly selecting epochs until we get one that is not
        # in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        # Now visualize the lens system at the epoch defined by the
        # randomIndex:
        desc.slrealizer.draw_model(self.observation[randomIndex],
                                   self.catalog.get_lens(lensID),
                                   convolve, debug)
        return

    def deblend(self, lensID=None, null_deblend = False, debug=False):
        if lensID is None:
            print 'No lens system selected for plotting.'
            return
        import random
        # Keep randomly selecting epochs until we get one that is not           
        # in the 'y' filter:                                                    
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
            # Now visualize the lens system at the epoch defined by the             
            # randomIndex:                                                          
        desc.slrealizer.deblend_test(self.observation[randomIndex],
                                   self.catalog.get_lens(lensID), null_deblend, debug)
        return

    # For now set all to true so that we can debug easily
    def deblend_distance(self, lensID=None, null_deblend=True, debug=True, show_plot=True):
        if lensID is None:
            print 'No lens system selected for calculating the statistics'
            return
        import random
        # Keep randomly selecting epochs until we get one that is not                                       
        # in the 'y' filter:                                         
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
            # Now visualize the lens system at the epoch defined by the                                      
            # randomIndex:                                          
        #if show_plot:
            #deblend(lensID, null_deblend, debug, randomIndex)
        if null_deblend:
            image1 = desc.slrealizer.null_deblending(self.observation[randomIndex], self.catalog.get_lens(lensID), debug)
        if show_plot:
            desc.slrealizer.show_color_map(image1)
        image2 = desc.slrealizer.plot_all_objects(self.observation[randomIndex], self.catalog.get_lens(lensID), debug)
        if show_plot:
            desc.slrealizer.show_color_map(image2)
        print('Chi squared distance is : ', desc.slrealizer.chi_square_distance(image1, image2))
        print('KL distance is : ', desc.slrealizer.KL_distance(image1, image2))
