"""
script that runs to generate toy catalogs
"""

import lenspop
import om10
import desc.slrealizer


db = om10.DB(vb=False, catalog='../../../data/qso_mock.fits')
db.paint(synthetic=True)
realizer = desc.slrealizer.SLRealizer(catalog=db, observation="../../../data/twinkles_observation_history.csv")
realizer.make_source_catalog()
realizer.make_object_catalog()
