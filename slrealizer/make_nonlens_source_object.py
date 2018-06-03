import os, sys
import galsim
import pandas as pd
import numpy as np

realizer_path = os.path.join(os.environ['SLREALIZERDIR'], 'slrealizer')
sys.path.insert(0, realizer_path)
from realize_sdss import SDSSRealizer
from utils.utils import *

if __name__=='__main__':
    data_path = os.path.join(os.environ['SLREALIZERDIR'], 'data')
    catalog_f = os.path.join(data_path, 'sdss_toy_processed.csv')
    observation_f = os.path.join(data_path, 'twinkles_observation_history.csv')

    db = pd.read_csv(catalog_f).reset_index(drop=True)
    obs = pd.read_csv(observation_f)\
            .query("(expMJD < 60919) & (filter != 'y')")\
            .reset_index(drop=True)
    realizer = SDSSRealizer(observation=obs, catalog=db, debug=False)
    table_path = os.path.join(data_path, 'nonlens_source_table.csv')

    realizer.make_source_table(save_file=table_path, use_hsm=False)
    realizer.make_object_table(sourceTablePath=table_path, objectTablePath='nonlens_object_table.csv')

    
