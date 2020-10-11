import numpy as np
import pylab
import scipy.io
import sys, getopt

if __name__ == "__main__":

    #plotData = np.loadtxt('/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/tSNEcoord_314.txt')
    plotData = np.loadtxt('/Users/mrahaman1/Documents/Statelet_V2/fResults/tSNEcoord_HC.txt')
    #print(plotData.shape)
    #lab = scipy.io.loadmat('/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/pDensity_allDOMshapes.mat')
    #labels = lab['p']
    #pylab.scatter(plotData[:, 0], plotData[:, 1], 20, labels)
    pylab.scatter(plotData[:, 0], plotData[:, 1])
    pylab.show()

