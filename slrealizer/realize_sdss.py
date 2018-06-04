#from __future__ import absolute_import
#from __future__ import division
#from __future__ import print_function

from realize_sl import SLRealizer
from utils.constants import *
from utils.utils import *
import pandas as pd
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
        self.sdss_pixel_scale = 0.396 #arcsec/pixel
        
    def get_lensInfo(self, objID=None, rownum=None):
        if objID is not None and rownum is not None:
            raise ValueError("Need to define either objID or rownum, not both.")
        
        if objID is not None:
            return self.catalog.get_lens(objID) # Don't use this (?)
        elif rownum is not None:
            return self.catalog.loc[rownum]
    
    def _sdss_to_galsim(self, lensInfo, band):
        raise NotImplementedError
    
    def _sdss_to_lsst(self, obsInfo, lensInfo):
        raise NotImplementedError
    
    def draw_system(self, obsInfo, lensInfo, save_dir=None):
        raise NotImplementedError
        
    def estimate_hsm(self, obsInfo, lensInfo):
        raise NotImplementedError
    
    def draw_emulated_system(self, obsInfo, lensInfo):
        raise NotImplementedError
    
    def create_source_row(self, obsInfo, lensInfo, use_hsm=False):
        
        histID, MJD, band, PSF_FWHM, sky_mag = obsInfo
        row = {}
        
        row['filter'] = band
        row['ccdVisitId'] = histID
        row['MJD'] = MJD
        row['psf_fwhm'] = PSF_FWHM
        row['apFluxErr'] = from_mag_to_flux(sky_mag-22.5)/5.0
        
        row['objectId'] = lensInfo['objectId']
        row['apFlux'] = lensInfo['modelFlux_' + band]
        row['x'] = np.cos(np.deg2rad(lensInfo['offsetDec_' + band]*3600.0)) * lensInfo['offsetRa_' + band]
        row['y'] = lensInfo['offsetDec_' + band]
        row['trace'] = lensInfo['mRrCc_' + band]*(self.sdss_pixel_scale)**2.0 + 2.0*fwhm_to_sigma(PSF_FWHM)**2.0
        row['e1'] = lensInfo['mE1_' + band]
        row['e2'] = lensInfo['mE2_' + band]
        
        return row
    
    def make_source_table_vectorized(self, save_path):
        import gc # need this to optimize memory usage
        import time
        
        start = time.time()

        ####################################
        # Merging catalog with observation #
        ####################################
        catalog = self.catalog.copy()
        observation = self.observation.copy()
        catalog['key'] = 0
        observation['key'] = 0
        src = catalog.merge(observation, how='left', on='key')
        src.drop('key', 1, inplace=True)
        gc.collect()
        
        ####################################
        # Collapsing multi-band properties #
        # into one of observed band        #
        ####################################
        setZeroDict = {}
        propsToCollapse = ['modelFlux', 'offsetRa', 'offsetDec', 'mRrCc', 'mE1', 'mE2', ]
        # Initialize dictionary of columns we want to collapse
        for p in propsToCollapse:
            setZeroDict[p] = [p + '_' + b for b in 'ugriz']
        # Set unused column values to zero
        for b in 'ugriz': # b = observed filter
            for p in propsToCollapse: # multi-band columns to collapse
                setZeroCols = setZeroDict[p][:]
                setZeroCols.remove(p + '_' + b)
                src.loc[src['filter'] == b, setZeroCols] = 0.0
        # Collapse
        for p in propsToCollapse: # multi-band columns to collapse
            src[p] = src[[p + '_' + b for b in 'ugriz']].sum(axis=1)
            src.drop([p + '_' + b for b in 'ugriz'], axis=1, inplace=True)
        gc.collect()
        
        ###################
        # Unit conversion #
        ###################
        src['apFluxErr'] = from_mag_to_flux(src['fiveSigmaDepth'] - 22.5)/5.0
        src['x'] = np.cos(np.deg2rad(src['offsetDec']*3600.0))*src['offsetRa']
        src['y'] = src['offsetDec']
        src['trace'] = src['mRrCc']*(self.sdss_pixel_scale**2.0) + 2.0*np.power(fwhm_to_sigma(src['FWHMeff']), 2.0)
        
        #####################################################
        # Final column renaming/reordering & saving to file #
        #####################################################
        src.rename(columns={'obsHistID': 'ccdVisitId',
                            'expMJD': 'MJD',
                            'FWHMeff': 'psf_fwhm',
                            'modelFlux': 'apFlux',
                            'mE1': 'e1',
                            'mE2': 'e2'}, inplace=True)
        src.drop(['mRrCc', 'offsetRa', 'offsetDec', 'fiveSigmaDepth'], axis=1, inplace=True)
        gc.collect()
        
        src = src[self.sourceCols]
        src.set_index('objectId', inplace=True)
        src.to_csv(save_path)
        end = time.time()
        
        print("Done making the source table with %d row(s) in %0.2f seconds using vectorization." %(len(src), end-start))
        
        self.sourceTable = src
        if self.DEBUG:
            return src
        
    
    #def make_source_table INHERITED
