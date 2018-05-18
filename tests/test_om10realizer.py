# *-* encoding: utf-8 *-*
# Unit tests for OM10Realizer functions

# ======================================================================
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
            print "Initializing testing environment"
            # run the real setup
            self.setupClass()
            # remember that it was setup already
            self.__class__.classIsSetup = True
                               
    def setupClass(self):
        # Do the real setup
        unittest.TestCase.setUp(self)
        # you want to have persistent things to test
        
        data_path = os.path.join(os.environ['SLREALIZERDIR'], 'data')
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
        self.__class__.test_datapath = os.path.join(os.environ['SLREALIZERDIR'], 'tests', 'test_source_table.csv')
    
    def test_from_om10_to_galsim(self):
        self.realizer._from_om10_to_galsim(lensInfo=self.test_lensInfo, band=self.test_obsInfo['filter'])
    
    def test_draw_system(self):
        self.realizer.draw_system(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo) 
    
    def test_estimate_hsm(self):
        self.realizer.estimate_hsm(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo)
    
    def test_create_source_row(self):
        print(self.realizer.create_source_row(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo))
        
    def test_make_source_table(self):
        self.realizer.make_source_table(save_dir=self.test_datapath)

if __name__ == '__main__':
    unittest.main()