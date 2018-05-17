from desc.slrealizer.realize_sl import SLRealizer
from desc.slrealizer.realize_om10 import OM10Realizer

if __name__=='__main__':
    r = OM10Realizer(1, 2)
    print r.observation
