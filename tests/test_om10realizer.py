# *-* encoding: utf-8 *-*
# Unit tests for OM10Realizer functions

# ======================================================================
from __future__ import print_function
import os, unittest
import galsim
from om10 import DB
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys
realizer_path = os.path.join(os.environ['SLREALIZERDIR'], 'slrealizer')
sys.path.insert(0, realizer_path)
from realize_om10 import OM10Realizer
from utils.utils import *
#from realize_sl import SLRealizer
# ======================================================================

class OM10RealizerTest(unittest.TestCase):
    
    """
    Tests the OM10 Realizer subclass.
    
    NOTE
    Execute these tests with:
        nosetests
    from anywhere in the module, provided you have run
        pip install nose
    """
    
    # Thanks to http://stezz.blogspot.com/2011/04/calling-only-once-setup-in-unittest-in.html
    # for the single setUp example!
    classIsSetup = False
    
    def setUp(self):
        # If it was not setup yet, do it
        if not self.classIsSetup:
            print("Initializing testing environment...")
            # run the real setup
            self.setupClass()
            # remember that it was setup already
            self.__class__.classIsSetup = True
                               
    def setupClass(self):
        # Do the real setup
        unittest.TestCase.setUp(self)
        # you want to have persistent things to test
        
        data_path = os.path.join(os.environ['SLREALIZERDIR'], 'data')
        test_path = os.path.join(os.environ['SLREALIZERDIR'], 'tests')
        test_catalog_f = os.path.join(data_path, 'test_catalog.fits')
        observation_f = os.path.join(data_path, 'twinkles_observation_history.csv')
    
        test_db = DB(catalog=test_catalog_f)
        test_db.paint(synthetic=True)
        test_obs = pd.read_csv(observation_f).sample(1, random_state=123).reset_index(drop=True)
        self.__class__.realizer = OM10Realizer(observation=test_obs, catalog=test_db, debug=True) 
        # We can access lens and observation rows through these variables.
        self.__class__.test_lensInfo = test_db.sample[0]
        self.__class__.test_obsInfo = test_obs.loc[0]
        # Path where we will save our test source table
        self.__class__.test_analytical_savepath = os.path.join(test_path, 'test_ana_source_table.csv')
        self.__class__.test_numerical_savepath = os.path.join(test_path, 'test_num_source_table.csv')
        self.__class__.test_vectorized_savepath = os.path.join(test_path, 'test_vectorized_source_table.csv')
        self.__class__.test_object_savepath = os.path.join(test_path, 'test_object_table.csv')
    
    def test_om10_to_galsim(self):
        self.realizer._om10_to_galsim(lensInfo=self.test_lensInfo, band=self.test_obsInfo['filter'])
    
    def test_draw_system(self):
        self.realizer.draw_system(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo) 
    
    def test_estimate_hsm(self):
        self.realizer.estimate_hsm(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo)

    def test_om10_to_lsst(self):
        self.realizer._om10_to_lsst(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo)
        
    def test_create_source_row_numerical(self):
        print("numerical: ", self.realizer.create_source_row(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo, use_hsm=True))
        
    def test_create_source_row_analytical(self):
        print("analytical: ", self.realizer.create_source_row(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo, use_hsm=False))
        
    def test_make_source_table_numerical(self):
        self.realizer.make_source_table(save_file=self.test_numerical_savepath, use_hsm=True)

    #def test_make_source_table_analytical(self):
    #    self.realizer.make_source_table(save_file=self.test_analytical_savepath, use_hsm=False)

    def test_vectorized_source_table(self):

        rowbyrow = self.realizer.make_source_table(save_file=self.test_analytical_savepath, use_hsm=False)[['e1', 'e2', 'trace', 'x', 'y', 'apFlux']].values.astype(float)
        #print("Row by row: ", rowbyrow)
        
        vectorized = self.realizer.make_source_table_vectorized(save_file=self.test_vectorized_savepath)[['e1', 'e2', 'trace', 'x', 'y', 'apFlux']].values.astype(float)
        #print("Vectorized: ", vectorized)
        
        self.assertTrue(np.allclose(rowbyrow, vectorized, rtol=1e-04, atol=1e-04))

    def test_make_object_table(self):
        self.realizer.make_object_table(sourceTablePath=self.test_vectorized_savepath, objectTablePath=self.test_object_savepath)

if __name__ == '__main__':
    unittest.main()
