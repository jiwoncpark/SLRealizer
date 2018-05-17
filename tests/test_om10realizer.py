# *-* encoding: utf-8 *-*
# Unit tests for OM10Realizer functions

# ======================================================================
import os, unittest
from desc.slrealizer.realize_om10 import OM10Realizer
import om10
import numpy as np
import pandas as pd
import galsim
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
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    
    catalog_f = os.path.join(curr_dir, '..', '..', '..', '..', 'data', 'qso_mock.fits')
    observation_f = os.path.join(curr_dir, '..', '..', '..', '..', 'data', 'twinkles_observation_history.csv')
    
    # Query the OM10 catalog to yield one lens system.
    test_db = om10.DB(catalog=catalog_f)
    test_db.select_random(maglim=23.3, area=10.0, IQ=0.75) 
    test_db.paint(synthetic=True)

    # Query the Twinkles catalog to yield one observation condition.
    test_obs = pd.read_csv(observation_f).sample(1, random_state=123).reset_index(drop=True)
    
    # Create OM10Realizer instance to be tested.
    realizer = OM10Realizer(observation=test_obs, catalog=test_db)
    
    # We can access lens and observation databases in dictionary form.
    db_dict = test_db.sample[0]
    obs_dict = test_obs.loc[0].to_dict()
        
    def test_hsm(self):
        

if __name__ == '__main__':
    unittest.main()