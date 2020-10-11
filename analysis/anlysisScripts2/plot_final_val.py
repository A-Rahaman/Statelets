import numpy as np
import pylab
import scipy.io
import sys, getopt
import matplotlib
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import scipy.io
import pylab as pl



if __name__ == "__main__":


    #for i in range(3, 7):

         #plotData = np.loadtxt('/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/tSNEcoord_314.txt')
         data = np.loadtxt('/Users/mrahaman1/Documents/Statelet_V2/fResults/tSNEcoord_314.txt')
         #print(plotData.shape)
         #lab2 = scipy.io.loadmat('/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/pDensity_allDOMshapes.mat')
         lab2 = scipy.io.loadmat('/Users/mrahaman1/Documents/Statelet_V2/fResults/pDensity_allDOMshapes.mat')
         lab = lab2['p']
         #pylab.scatter(plotData[:, 0], plotData[:, 1], 20, labels)
         #pylab.scatter(plotData[:, 0], plotData[:, 1])
         #pylab.show()
         #j = "%03d"%i;
         #data = np.loadtxt('/Users/dsaha/Documents/tsne_python_2/Y_Subject_' + str(j) + '.txt');

         #lab2 = scipy.io.loadmat('/Users/dsaha/Documents/tsne_python_2/label_subject_' + str(i) + '.mat')
         #lab = lab2['labels']

         sizeX,sizeY = data.shape



         (Shared_length_X, Shared_length_Y) = data.shape;
         lab = lab.reshape(Shared_length_X, );
         viri = cm.get_cmap('inferno')
         rgbaValues = viri(lab)
         pl.scatter(data[:, 0], data[:, 1], color=rgbaValues)
         pl.show()
         pl.axis('off')
         pl.legend(loc=0)
         pl.title('tSNE plot')
         pl.savefig('SUB_' + '%03d' % i, format='png', pad_inches=0, bbox_inches='tight')
         pl.close()

