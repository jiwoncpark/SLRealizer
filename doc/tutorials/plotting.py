import pylab
import matplotlib
import math
import matplotlib.pyplot as plt

def plot_lens(lens, obsHist, convolve=False):
	fig = plt.figure()
	scale_factor = 2
	print obsHist
	sourceX = lens['XSRC'][0]
	sourceY = lens['YSRC'][0]
	source = plt.Circle((sourceX, sourceY), radius=obsHist[2]/scale_factor, alpha=0.1, fc='r')
	lensX = lens['XIMG'][0]
	lensY = lens['YIMG'][0]
	plotX = 0
	plotY = 0
	# REALLY BAD STYLE JENNY YOU WOULD HAVE GOTTEN C FOR THIS
	print abs(lensX)
	print abs(sourceX)
	if (abs(lensX[0]))>(abs(sourceX)):
		plotX = abs(lensX[0])
	else:
		plotX = abs(sourceX)
	if (abs(lensY[0]))>(abs(sourceY)):
		plotY = abs(lensY[0])
	else:
		plotY = abs(sourceY)
	seeing = plt.Circle((-plotY-obsHist[2]/scale_factor, -plotX-obsHist[2]/scale_factor), radius=obsHist[2]/scale_factor, alpha=0.1, fc='g')	
	plt.ylim(-plotY-2*obsHist[2]/scale_factor, plotY+2*obsHist[2]/scale_factor)
	plt.xlim(-plotX-2*obsHist[2]/scale_factor, plotX+2*obsHist[2]/scale_factor)
	lens = plt.Circle((lensX[0], lensY[0]), radius=obsHist[2]/scale_factor, alpha=0.1, fc='b')	
	plt.gca().add_patch(seeing)
	plt.gca().add_patch(lens)
	plt.gca().add_patch(source)

