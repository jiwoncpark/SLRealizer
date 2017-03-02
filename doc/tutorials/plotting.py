import pylab
import matplotlib
import math
import matplotlib.pyplot as plt

def plot_lens(lens, obsHist, convolve=False):
    #obsHist has MJD Filter FWHM 5sigmag
    filterQuasar = obsHist[1] + '_SDSS_quasar'
    filterLens = obsHist[1] + '_SDSS_Lens'
    scale_factor = 2
    sourceX = lens['XSRC'][0]
    sourceY = lens['YSRC'][0]
    lensX = lens['XIMG'][0]
    lensY = lens['YIMG'][0]
    if (abs(lensX[0]+sourceX))>(abs(sourceX)):
        plotX = abs(lensX[0]+sourceX)
    else:
        plotX = abs(sourceX)
    if (abs(lensY[0]+sourceY))>(abs(sourceY)):
        plotY = abs(lensY[0]+sourceY)
    else:
        plotY = abs(sourceY)
    fig = plt.figure(figsize=(4, 16))

    if(True):
        for i in range(4):
            plotNumStr = "4"+"1"+str(i+1)
            plotNum = int(plotNumStr)
            sub = fig.add_subplot(plotNum)
            sub.set_ylim(-plotY-2*obsHist[2]/scale_factor, plotY+2*obsHist[2]/scale_factor)
            sub.set_xlim(-plotX-2*obsHist[2]/scale_factor, plotX+2*obsHist[2]/scale_factor)
            source = plt.Circle((sourceX, sourceY), radius=obsHist[2]/scale_factor, alpha=0.1, fc='r')
            lens = plt.Circle((lensX[i]+sourceX, lensY[i]+sourceY), radius=obsHist[2]/scale_factor, alpha=0.1, fc='b')
            fig.gca().add_patch(lens)
            fig.gca().add_patch(source)
            seeing = plt.Circle((-plotY+(obsHist[2]*0.7)/scale_factor, -plotX-(obsHist[2]*0.8)/scale_factor), radius=obsHist[2]/scale_factor, alpha=0.1, fc='g')
            sub.set_title('Observation ' + str(i+1) + ' with ' + 'filter')
            sub.legend((source, seeing, lens),('source', 'seeing', 'lens'))
            plt.gca().add_patch(seeing)
