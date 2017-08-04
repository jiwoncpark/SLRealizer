import lenspop
import om10
import desc.slrealizer

def generate_catalog():
    db = om10.DB()
    db.paint(synthetic=True)
    realizer = desc.slrealizer.SLRealizer(catalog=db, observation="../../../data/twinkles_observation_history.csv")
    realizer.make_catalog(num_system = 3)
