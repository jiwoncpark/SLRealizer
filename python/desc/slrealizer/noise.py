#hello

"""An exposure time calculator for LSST.  Uses GalSim to draw a galaxy with specified magnitude,
shape, etc, and then uses the same image as the optimal weight function.  Derived from D. Kirkby's
notes on deblending.
"""
from __future__ import print_function

import numpy as np

import galsim

# Some constants
# --------------
#
# LSST effective area in meters^2
A = 319/9.6  # etendue / FoV.  I *think* this includes vignetting

# zeropoints from DK notes in photons per second per pixel
# should eventually compute these on the fly from filter throughput functions.
s0 = {'u': A*0.732,
      'g': A*2.124,
      'r': A*1.681,
      'i': A*1.249,
      'z': A*0.862,
      'Y': A*0.452}

# exposure time per visit
# For LSST it is 15.0 sec
visit_time = 30
# Sky brightness per arcsec^2 per second

def return_noise(sky_mag, band):
    pixel_scale = 0.2
    sky = s0[band] * 10**(-0.4*(sky_mag-24.0)) * visit_time * pixel_scale**2
    return np.sqrt(sky)
