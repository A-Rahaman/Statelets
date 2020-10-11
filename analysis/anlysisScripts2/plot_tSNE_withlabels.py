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


    for i in range(1, 2):

        j = "%03d"%i;
        #data = np.loadtxt('tSNEcoord_314.txt');
        data = np.loadtxt('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/tSNE_allpairshape_HC.txt');

        lab2 = scipy.io.loadmat('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/allPairshapesPD_HC.mat')
        lab = lab2['p']


        (Shared_length_X, Shared_length_Y) = data.shape;
        lab = lab.reshape(Shared_length_X, );
        #lab = (lab- np.min(lab))/(np.max(lab)-np.min(lab))

        idx = np.where(lab == 1)
        hc_color = np.asarray([[0, 0, 1, 0.5]] * len(idx[0]))

        idx1 = np.where(lab == 0)
        sz_color = np.asarray([[1, 0, 0, 0.5]] * len(idx1[0]))

        viri = cm.get_cmap('inferno')
        rgbaValues = viri(lab)
        rgbaValues[idx] = hc_color;
        rgbaValues[idx1] = sz_color;


        pl.scatter(data[:, 0], data[:, 1], color=rgbaValues)
        #pl.colorbar()
        pl.savefig('plot_final', format='png', pad_inches=0, bbox_inches='tight')
        pl.show()
        pl.axis('off')
        pl.legend(loc=0)
        pl.title('tSNE plot')

        #pl.savefig('SUB_' + '%03d' % i, format='png', pad_inches=0, bbox_inches='tight')
        #pl.savefig('savePic', format='png', pad_inches=0, bbox_inches='tight')
        #pl.savefig('tSNE_cord_314',format='png', pad_inches=0, bbox_inches='tight')
        #pl.close()

