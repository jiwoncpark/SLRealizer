#from __future__ import absolute_import
#from __future__ import division
#from __future__ import print_function

from realize_sl import SLRealizer
from utils.constants import *
from utils.utils import *
import pandas
import numpy as np
import galsim

class OM10Realizer(SLRealizer):
    
    def __init__(self, observation, catalog, debug=False):
        #super(OM10Realizer, self).__init__(observation) # Didn't work for some reason
        self.as_super = super(OM10Realizer, self)
        self.as_super.__init__(observation)
        self.catalog = catalog
        self.num_systems = len(self.catalog.sample)
        self.DEBUG = debug
        
    def get_lensInfo(self, lensID=None, rownum=None):
        if lensID is not None and rownum is not None:
            raise ValueError("Need to define either lensID or rownum, not both.")
        
        if lensID is not None:
            return self.catalog.get_lens(lensID) # Don't use this (?)
        elif rownum is not None:
            return self.catalog.sample[rownum]
       
    def _from_om10_to_galsim(self, lensInfo, band):
        """
        Converts OM10's column values into GalSim terms
        
        Keyword arguments:
        lensInfo -- a row of the OM10 DB
        band -- the filter used to observe

        Returns:
        A dictionary (named galsimInput) containing properties that 
        can be passed into GalSim (See code for which properties)
        """
        
        # We will input flux in units of nMgy
        mag   = lensInfo[band + '_SDSS_lens'] 
        flux  = from_mag_to_flux(mag, to_unit='nMgy') # in nMgy
        hlr   = lensInfo['REFF_T'] # REFF_T is in arcsec
        e     = lensInfo['ELLIP']
        beta  = lensInfo['PHIE'] * galsim.degrees # PHIE is in degrees
        
        galsimInput={'flux': flux,
                       'half_light_radius': hlr,
                       'e': e,
                       'beta': beta,
                       'num_objects': lensInfo['NIMG']}
        
        for obj in xrange(lensInfo['NIMG']):
            obj_mag = lensInfo[band + '_SDSS_quasar']\
                      + from_flux_to_mag(abs(lensInfo['MAG'][obj]))
                      #+ from_flux_to_mag(lensInfo['MAG'][obj] + get_filter_AB_offset())
            galsimInput['flux_'+str(obj)] = from_mag_to_flux(obj_mag, to_unit='nMgy') # don't subtract AB offset?
            galsimInput['xy_'+str(obj)] = lensInfo['XIMG'][obj], lensInfo['YIMG'][obj]
            
        return galsimInput        
            
    def draw_system(self, obsInfo, lensInfo, save_dir=None):
        galsimInput = self._from_om10_to_galsim(lensInfo, obsInfo['filter'])
        return self.as_super.draw_system(galsimInput=galsimInput, obsInfo=obsInfo, save_dir=save_dir)
    
    def estimate_hsm(self, obsInfo, lensInfo):
        galsim_img = self.draw_system(lensInfo=lensInfo, obsInfo=obsInfo, save_dir=None)
        return self.as_super.estimate_hsm(galsim_img=galsim_img, obsInfo=obsInfo)
    
    def draw_emulated_system(self, obsInfo, lensInfo):
        hsmOutput = self.estimate_hsm(obsInfo, lensInfo)
        return self.as_super.draw_emulated_system(hsmOutput)
    
    def create_source_row(self, obsInfo, lensInfo):
        hsmOutput = self.estimate_hsm(obsInfo, lensInfo)
        if hsmOutput is None:
            return None
        objectId = lensInfo['LENSID']
        return self.as_super.create_source_row(hsmOutput=hsmOutput, objectId=objectId, obsInfo=obsInfo)
    
    #def make_source_table(self, save_file)