import lenspop
import om10
import desc.slrealizer


db = om10.DB()
db.paint(synthetic=True)
realizer = desc.slrealizer.SLRealizer(catalog=db, observation="../../../data/twinkles_observation_history.csv")
realizer.make_source_catalog()
realizer.make_object_catalog()
