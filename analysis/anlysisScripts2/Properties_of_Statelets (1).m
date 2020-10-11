clear
clc
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts/d_Statelets.mat
load /export/mialab/users/mrahaman/StateLet/Data/SZsubjectIdx.mat
load /export/mialab/users/mrahaman/StateLet/Data/HCsubjectIdx.mat
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts/Wavelets_314X1081.mat
load '/export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts/AllFreqs.mat'
%%
%{

% Load and combine the frquencies 
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/ShapeletFreq_HC.mat
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/ShapeletFreq_SZ.mat
Freqs = zeros(314,1081);
for hc =1:length(CtIdx)
    Freqs(CtIdx(hc),:)= HCFreq(:,hc);
end
for sz =1:length(SZIdx)
    Freqs(SZIdx(sz),:)= SZFreq(:,sz);
end
%}

for i =1:size(Statelets,1)
    Subjects_i = unique(Statelets{i,1});
    Pairs_i    = Statelets{i,3};
    HCs        = intersect(Subjects_i,CtIdx);
    SZs        = intersect(Subjects_i,SZIdx); 
    mConn_pval_tval = zeros(1081,2);
    len_pval_tval = zeros(1081,2);
    freq_pval_tval = zeros(1081,2);
    
    for j=1:length(Pairs_i)
       SZmean = zeros(length(SZs),1);
       SZlen = zeros(length(SZs),1);
       for k = 1:length(SZs)
           % Have to do something by using the generalized wavelets. whic shapelet of the 5 covers that pair.     
           SZmean(k) = mean(Wavelets{SZs(k),Pairs_i(j)});
           SZlen(k)  = length(Wavelets{SZs(k),Pairs_i(j)});
       end
       HCmean = zeros(length(HCs),1);
       HClen = zeros(length(HCs),1);
       for k = 1:length(HCs)
           % Have to do something by using the generalized wavelets. whic shapelet of the 5 covers that pair.     
           HCmean(k) = mean(Wavelets{HCs(k),Pairs_i(j)});
           HClen(k)  = length(Wavelets{HCs(k),Pairs_i(j)});
       end
       [~,p,~,tstats] = ttest2(HCmean,SZmean);                                   % ttest on mean shapelet
       mConn_pval_tval(Pairs_i(j),1) = p;
       mConn_pval_tval(Pairs_i(j),2) = tstats.tstat; 
       
       [~,p1,~,tstats1] = ttest2(HClen,SZlen);                                   % ttest on length
       len_pval_tval(Pairs_i(j),1) = p1;
       len_pval_tval(Pairs_i(j),2) = tstats1.tstat;
      
       [~,p2,~,tstats2] = ttest2(Freqs(HCs,Pairs_i(j)),Freqs(SZs,Pairs_i(j)));   % ttest on frequency
       freq_pval_tval(Pairs_i(j),1) = p2;
       freq_pval_tval(Pairs_i(j),2) = tstats2.tstat;
      
    end
    pval_tval{i,1} = mConn_pval_tval; % Mean pair value 
    pval_tval{i,2} = len_pval_tval;   % Length of shapelets 
    pval_tval{i,3} = freq_pval_tval;  % Frequency of shapelets 
end

%{
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts/Subject_wise_shIndex.mat
for i =1:size(Statelets,1)
    Subjects_i = unique(Statelets{i,1});
    Pairs_i    = Statelets{i,3};
    HCs        = intersect(Subjects_i,CtIdx);
    SZs        = intersect(Subjects_i,SZIdx);
    comm_pairs_HC = SUB_shLets_idx(HCs(1),:);
    for k = 2:length(HCs)
    comm_pairs_HC = intersect(comm_pairs_HC,SUB_shLets_idx(HCs(k),:)); 
    end
end
%}

