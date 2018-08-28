import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from om10 import DB

d = DB()
print(dir(om10))
d.paint(synthetic=True)
