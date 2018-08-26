from __future__ import absolute_import

# Assumes sl-realizer is package root
import slrealizer
print(dir(slrealizer))
yo = slrealizer.ThisDoer(3)
#from slrealizer.do_this import ThisDoer
# Assumes slrealizer is package root
#from do_this import ThisDoer
yo = slrealizer.do_this.ThisDoer(3)
yo.get_this()

import os.path
# Assumes sl-realizer is package root
print('SLRealizer/slrealizer/data/data.py', os.path.exists('SLRealizer/slrealizer/data/data.py'))
# Assumes slrealizer is package root
print('data/data.py', os.path.exists('data/data.py'))
