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
from realize_sdss import SDSSRealizer
from utils.utils import *
#from realize_sl import SLRealizer
# ======================================================================

class SDSSRealizerTest(unittest.TestCase):
    
    """
    Tests the SDSS Realizer subclass.
    
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
        catalog_f = os.path.join(data_path, 'sdss_test_processed.csv')
        observation_f = os.path.join(data_path, 'twinkles_observation_history.csv')
    
        test_db = pd.read_csv(catalog_f).sample(1, random_state=123).reset_index(drop=True)
        test_obs = pd.read_csv(observation_f).sample(1, random_state=123).reset_index(drop=True)
        
        self.__class__.realizer = SDSSRealizer(observation=test_obs, catalog=test_db, debug=False) 
        # We can access lens and observation rows through these variables.
        self.__class__.test_lensInfo = test_db.loc[0]
        self.__class__.test_obsInfo = test_obs.loc[0]
        # Path where we will save our test source table
        self.__class__.test_datapath = os.path.join(os.environ['SLREALIZERDIR'], 'tests', 'test_nonlens_source_table.csv')
    
    def test_create_source_row(self):
        print(self.realizer.create_source_row(lensInfo=self.test_lensInfo, obsInfo=self.test_obsInfo))
        
    def test_make_source_table(self):
        self.realizer.make_source_table(save_file=self.test_datapath)

if __name__ == '__main__':
    unittest.main()
