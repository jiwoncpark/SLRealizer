"""
script that runs to generate toy catalogs
"""

import lenspop
import om10
import desc.slrealizer


db = om10.DB(vb=False, catalog='../../../data/qso_mock.fits')
db.select_random(maglim=23.3,area=20000.0,IQ=0.75)
db.paint(synthetic=True)
realizer = desc.slrealizer.SLRealizer(catalog=db, observation="../../../data/twinkles_observation_history.csv")
realizer.make_source_catalog(dir='../../../data/source_catalog_galsim_noise.csv')
realizer.make_object_catalog(source_table_dir='../../../data/source_catalog_galsim_noise.csv', save_dir='../../../data/object_catalog_galsim_noise.csv')
#realizer.overlap()
#realizer.make_source_catalog(galsim=False, dir='../../../data/source_catalog_c.csv')
#realizer.make_object_catalog(source_table_dir='../../../data/source_catalog_c.csv', save_dir='../../../data/object_catalog_c.csv')
