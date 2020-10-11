import scipy.io
from scipy import stats



if __name__ == "__main__":
         
         rankSZ = scipy.io.loadmat('/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/PairRanking_sSZ_HC.mat') # SZ sorted
         sortedSZ = rankSZ['ranking_sSZ_HC'];
         tauSZ, p_valueSZ = stats.kendalltau(sortedSZ[:,0], sortedSZ[:,1]);

         rankHC = scipy.io.loadmat('/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/PairRanking_sHC_SZ.mat') # HC sorted
         sortedHC = rankHC['ranking_sHC_SZ'];
         tauHC, p_valueHC = stats.kendalltau(sortedHC[:,0], sortedHC[:,1]);
         
         tau, p_value = stats.kendalltau(sortedHC[:,1], sortedSZ[:,1]);
         
         print("SZ sorted tau: %f , p-value is %f" % (tauSZ, p_valueSZ))
         print("HC sorted tau: %f , p-value is %f" % (tauHC, p_valueHC))
         print("overall tau: %f and p-value is %f" % (tau, p_value))
         
         
