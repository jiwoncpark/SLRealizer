#from __future__ import absolute_import
#from __future__ import division
#from __future__ import print_function

from realize_sl import SLRealizer
from utils.constants import *
from utils.utils import *
import pandas
import numpy as np
import galsim

class SDSSRealizer(SLRealizer):
    
    def __init__(self, observation, catalog, debug=False):
        #super(SDSSRealizer, self).__init__(observation) # Didn't work for some reason
        self.as_super = super(SDSSRealizer, self)
        self.as_super.__init__(observation)
        self.catalog = catalog
        self.num_systems = len(self.catalog)
        self.DEBUG = debug
        
    def get_lensInfo(self, objID=None, rownum=None):
        if objID is not None and rownum is not None:
            raise ValueError("Need to define either objID or rownum, not both.")
        
        if objID is not None:
            return self.catalog.get_lens(objID) # Don't use this (?)
        elif rownum is not None:
            return self.catalog.loc[rownum]
    
    def create_source_row(self, obsInfo, lensInfo):
        
        histID, MJD, band, PSF_FWHM, sky_mag = obsInfo
        row = {}
        
        row['filter'] = band
        row['MJD'] = MJD
        row['psf_fwhm'] = PSF_FWHM
        row['skyErr'] = from_mag_to_flux(sky_mag-22.5)/5.0
        
        row['objectId'] = lensInfo['objectId']
        row['appFlux'] = lensInfo['modelFlux_' + band]
        row['x'] = np.cos(np.deg2rad(lensInfo['offsetDec_' + band]*3600.0)) * lensInfo['offsetRa_' + band]
        row['y'] = lensInfo['offsetDec_' + band]
        row['trace'] = lensInfo['mRrCc_' + band] + 2.0*fwhm_to_sigma(PSF_FWHM)**2.0
        row['e1'] = lensInfo['mE1_' + band]
        row['e2'] = lensInfo['mE2_' + band]
        
        return row
    
    #def make_source_table INHERITED
