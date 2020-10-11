clear
clc
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/ShapeletFreq_HC.mat
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/ShapeletFreq_SZ.mat
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/ShapeletMeanLength_HC.mat
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/ShapeletMeanLength_SZ.mat

f_pval_tval = zeros(size(HCFreq,1),2);
L_pval_tval = zeros(size(HCavglen,1),2);

for i = 1:size(HCFreq,1)
[~,p,~,tstats] = ttest2(HCFreq(i,:),SZFreq(i,:));
f_pval_tval(i,1) = p;
f_pval_tval(i,2) = tstats.tstat; 
[~,lp,~,ltstats] = ttest2(HCavglen(i,:),SZavglen(i,:));
L_pval_tval(i,1) = lp;
L_pval_tval(i,2) = ltstats.tstat;
end