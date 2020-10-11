clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/Statelets.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/SZsubjectIdx.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/HCsubjectIdx.mat

load /Users/mrahaman1/Documents/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts/AllFreqs.mat
load /Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/REP_replicability_of_Leaders.mat
load /Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/REP_SUB_wise_leaderShapelets_uLen.mat
%% get the frequencies and the average lengths 
clear
clc
Freqs  = zeros(314,1081);
avglen = zeros(314,1081);
for i =1:314
    load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/subinfo',['info_SUB_',num2str(i,'%03.f')]));
    for j=1:1081
        Freqs(i,j)  = repWavelet_Details{j,5};
        avglen(i,j) = repWavelet_Details{j,6};
    end
end
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
    len_pval_tval   = zeros(1081,2);
    freq_pval_tval  = zeros(1081,2);
    
    for j=1:length(Pairs_i) % For a given pair, pair_i(j)
       SZmean = zeros(length(SZs),1);
       SZlen = zeros(length(SZs),1);
       for k = 1:length(SZs)
           %%%%%% Which shapelet of the 5 covers that pair.  
           for re = 1:5
           if(length(intersect(Pairs_i(j),reps_COMM{SZs(k),re}))>=1)
               break;
           end
           end
           SZmean(k) = mean(SUB_wise_shapelets_uLen{SZs(k),re});
           SZlen(k)  = length(SUB_wise_shapelets_uLen{SZs(k),re});
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           %SZmean(k) = mean(Wavelets{SZs(k),Pairs_i(j)});
           %SZlen(k)  = length(Wavelets{SZs(k),Pairs_i(j)});
       end
       HCmean = zeros(length(HCs),1);
       HClen = zeros(length(HCs),1);
       for k = 1:length(HCs)
           % Have to do something by using the generalized wavelets. which shapelet of the 5 covers that pair.     
           %%%%%% Which shapelet of the 5 covers that pair.  
           for re = 1:5
           if(length(intersect(Pairs_i(j),reps_COMM{HCs(k),re}))>=1)
               break;
           end
           end
           HCmean(k) = mean(SUB_wise_shapelets_uLen{HCs(k),re});
           HClen(k)  = length(SUB_wise_shapelets_uLen{HCs(k),re});
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %HCmean(k) = mean(Wavelets{HCs(k),Pairs_i(j)});
           %HClen(k)  = length(Wavelets{HCs(k),Pairs_i(j)});
       end
       [~,p,~,tstats] = ttest2(HCmean,SZmean);                                   % ttest on mean shapelet/ has to make an maxshaplet/maxconnectivity
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
%set(gca,'YTickLabel',[]);

%%