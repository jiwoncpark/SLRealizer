from slrealize_base import SLRealizer
import utils
import constants
import pandas
import numpy as np
import galsim

class SDSSRealizer(SLRealizer):

    def __init__(self, catalog, observation):
        super(SDSSRealizer, self).__init__(observation=observation)
        self._catalog = catalog
        self._num_systems = len(self._catalog.sample)
        
    def _get_lens(self, lensID=None, rownum=None):
        if lensID is not None and rownum is not None:
            raise ValueError("Need to define either lensID or rownum, not both.")
        
        if lensID is not None:
            return self._catalog.get_lens(lensID)
        elif rownum is not None:
            return self._catalog.get_lens(self._catalog.sample[i])
    
    def draw_lens_random_date(self, lensID=None, rownum=None, save_dir=None):
        """                                                                                                                   
        Draws a lens with the given lensID under randomly chosen observation conditions

        Keyword arguments:
        - lensID: an integer specifying the ID of the lens system in the OM10 catalog

        Returns:
        - A dictionary of the properties of the observed lens
        # TODO say which properties

        """
        # Randomly select observation ID
        obsID = random.randint(0, self.num_obs)
        # Render an image of the system using GalSim
        img, obj = nd.draw_system(self._get_lens(lensID=lensID, rownum=rownum),\
                                  self._get_observation(obsID),\
                                  save_dir)
        # Draw the lens system at the epoch defined by obsID
        fig, axes = plt.subplots(2, figsize=(5, 10))
        axes[0].imshow(img.array, interpolation='none', aspect='auto')
        axes[0].set_title("TRUE MODEL IMAGE, epochID: %d" %obsID)
        params = self.get_lens_params(self.catalog.get_lens(lensID),\
                                       self.observation.loc[obsID])
        # Define a GalSim object using the properties derived from the aggregate system 
        galaxy = galsim.Gaussian(flux=params['flux'], sigma=params['size'])\
                       .shift(float(params['x']), float(params['y']))\
                       .shear(e1=params['e1'], e2=params['e2'])
        # Draw the emulated image
        img = galaxy.drawImage(scale=utils.get_pixel_scale(), method='no_pixel')
        axes[1].imshow(img.array, interpolation='none', aspect='auto')
        axes[1].set_title("EMULATED IMAGE, epochID %d" %obsID)
        if save_dir is not None:
            plt.savefig(save_dir + 'deblending.png')
        return params
    
    