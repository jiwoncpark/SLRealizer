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
    catalog_f = os.path.join(data_path, 'sdss_processed.csv')
    observation_f = os.path.join(data_path, 'twinkles_observation_history.csv')

    output_nonlens_source_path = os.path.join(data_path, 'nonlens_source_table.csv')
    output_nonlens_object_path = os.path.join(data_path, 'nonlens_object_table.csv')
    
    db = pd.read_csv(catalog_f).sample(20000, random_state=123).query('(u_trace < 5.12)').reset_index(drop=True)
    obs = pd.read_csv(observation_f)\
            .query("(expMJD < 62450) & (filter != 'y')")\
            .reset_index(drop=True)
    realizer = SDSSRealizer(observation=obs, catalog=db, debug=False)

    realizer.make_source_table_vectorized(save_path=output_nonlens_source_path)
    realizer.make_object_table(sourceTablePath=output_nonlens_source_path,
                               objectTablePath=output_nonlens_object_path)

    
