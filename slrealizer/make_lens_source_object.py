from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
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

if __name__=='__main__':
    data_path = os.path.join(os.environ['SLREALIZERDIR'], 'data')
    catalog_f = os.path.join(data_path, 'qso_mock.fits')
    observation_f = os.path.join(data_path, 'twinkles_observation_history.csv')

    output_lens_source_path = os.path.join(data_path, 'lens_source_table.csv')
    output_lens_object_path = os.path.join(data_path, 'lens_object_table.csv')
    output_lens_tvar_source_path = os.path.join(data_path, 'lens_tvar_source_table.csv')

    db = DB(catalog=catalog_f)
    db.select_random(maglim=23.3, area=1000, IQ=0.75)
    db.paint(synthetic=True)
    
    obs = pd.read_csv(observation_f)\
            .query("(expMJD < 61000) & (filter != 'y')")\
            .reset_index(drop=True)
    realizer = OM10Realizer(observation=obs, catalog=db, debug=True)

    realizer.make_source_table_vectorized(save_file=output_lens_source_path)
    realizer.make_object_table(sourceTablePath=output_lens_source_path,
                               objectTablePath=output_lens_object_path)
    
    # Optionally add time variability
    realizer.add_time_variability(input_source_path=output_lens_source_path,
                                  output_source_path=output_lens_tvar_source_path)
    
