clear
clc
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/Statelets.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/Statelets_k10.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/SZsubjectIdx.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/HCsubjectIdx.mat
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
%% ttest on the max conn, mean conn, length, number of occurences of the shapelet for HC-SZ group differences
load /Users/mrahaman1/Documents/Statelet_V2/fResults/avglen.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/Occurences.mat
for i =1:size(Statelets,1)
    
    Subjects_i = Statelets{i,1};
    Pairs_i    = 1:1081;%Statelets{i,3};%all
    HCs        = intersect(Subjects_i,CtIdx);
    SZs        = intersect(Subjects_i,SZIdx);
    
    maxConn_pval_tval = zeros(1081,2);
    len_pval_tval   = zeros(1081,2);
    freq_pval_tval  = zeros(1081,2);
    mConn_pval_tval = zeros(1081,2);
    
    for j=1:length(Pairs_i) % For a given pair, pair_i(j)
       j
       SZmean = zeros(length(SZs),1);
       SZmax = zeros(length(SZs),1);
       for k = 1:length(SZs)
           load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/pairWISELeaders',['sub_',num2str(SZs(k),'%03.f')]));
           SZmax(k)     = max(pairWISELeadershapes{Pairs_i(j)});
           SZmean(k)    = mean(pairWISELeadershapes{Pairs_i(j)});
       end
       HCmean = zeros(length(HCs),1);
       HCmax = zeros(length(HCs),1);
       for k = 1:length(HCs)
           load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/pairWISELeaders',['sub_',num2str(HCs(k),'%03.f')]));
           HCmax(k)     = max(pairWISELeadershapes{Pairs_i(j)});
           HCmean(k)    = mean(pairWISELeadershapes{Pairs_i(j)});
       end
       
       % ttests 1: max amplitude of the shapelet
       [~,p,~,tstats] = ttest2(HCmean,SZmean);                                  % ttest on maxshaplet/maxconnectivity
       maxConn_pval_tval(Pairs_i(j),1) = p;
       maxConn_pval_tval(Pairs_i(j),2) = tstats.tstat; 
       
       [~,p,~,tstats] = ttest2(HCmax,SZmax);                                  % ttest on mean shapelet/mean connectivity
       mConn_pval_tval(Pairs_i(j),1) = p;
       mConn_pval_tval(Pairs_i(j),2) = tstats.tstat; 
       
       [~,p1,~,tstats1] = ttest2(avglen(HCs,Pairs_i(j)),avglen(SZs,Pairs_i(j))); % ttest on length
       len_pval_tval(Pairs_i(j),1) = p1;
       len_pval_tval(Pairs_i(j),2) = tstats1.tstat;
      
       [~,p2,~,tstats2] = ttest2(Freqs(HCs,Pairs_i(j)),Freqs(SZs,Pairs_i(j)));   % ttest on frequency
       freq_pval_tval(Pairs_i(j),1) = p2;
       freq_pval_tval(Pairs_i(j),2) = tstats2.tstat;
      
    end
    
    pval_tval{i,1} = mConn_pval_tval;   % Mean pair value 
    pval_tval{i,2} = maxConn_pval_tval; % Max pair value 
    pval_tval{i,3} = len_pval_tval;     % Length of shapelets 
    pval_tval{i,4} = freq_pval_tval;    % Frequency of shapelets 
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

%% based on the real shapes included in each statelet
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/Statelets.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes_real.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/shapestosubjectclassMapping.mat % shapes represnts which group SZ or HC
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allshapestoPairmapping.mat % shapes coming from which pair 

for i =1:size(Statelets,1) 
    shapes = Statelets{i,2}; 
    pairs  = shape_to_pair(shapes);
    subcls = shape_to_subcl(shapes);

    maxConn_pval_tval = zeros(1081,2);
    len_pval_tval   = zeros(1081,2);
    %freq_pval_tval  = zeros(1081,2);
    mConn_pval_tval = zeros(1081,2);
    Pairs_i = unique(pairs);
    
    for j=1:length(Pairs_i) % For a given pair, pair_i(j)
        
       shape_idxPJ = find(pairs==Pairs_i(j));
       
       %SZ
       SZmax = [];
       SZmean = [];
       SZlen =[];
       szIdx = find(subcls(shape_idxPJ)==0); 
       SZshapes_PJ = shape_idxPJ(szIdx);
       for k =1:length(SZshapes_PJ)
           SZmax(k)     = max(allSUBshapes_real{SZshapes_PJ(k)});
           SZmean(k)    = mean(allSUBshapes_real{SZshapes_PJ(k)});
           SZlen(k)     = length(allSUBshapes_real{SZshapes_PJ(k)});
       end
       % HC
       HCmax = [];
       HCmean = [];
       HClen =[];
       hcIdx = find(subcls(shape_idxPJ)==1); 
       HCshapes_PJ = shape_idxPJ(hcIdx);
       for k =1:length(HCshapes_PJ)
           HCmax(k)     = max(allSUBshapes_real{HCshapes_PJ(k)});
           HCmean(k)    = mean(allSUBshapes_real{HCshapes_PJ(k)});
           HClen(k)     = length(allSUBshapes_real{HCshapes_PJ(k)});
       end
       if(isempty(HCmean))
           HCmean=zeros(length(SZmean),1);
       elseif(isempty(SZmean))
           SZmean=zeros(length(HCmean),1);
       end
       if(isempty(HCmax))
           HCmax=zeros(length(SZmax),1);
       elseif(isempty(SZmax))
           SZmax=zeros(length(HCmax),1);
       end
       if(isempty(HClen))
           HClen=zeros(length(SZlen),1);
       elseif(isempty(SZlen))
           SZlen=zeros(length(HClen),1);
       end
       % ttests 1: max amplitude of the shapelet
       [~,p,~,tstats] = ttest2(HCmax,SZmax);                                  % ttest on maxshaplet/maxconnectivity
       maxConn_pval_tval(Pairs_i(j),1) = p;
       maxConn_pval_tval(Pairs_i(j),2) = tstats.tstat; 
       
       [~,p,~,tstats] = ttest2(HCmean,SZmean);                                  % ttest on mean shapelet/mean connectivity
       mConn_pval_tval(Pairs_i(j),1) = p;
       mConn_pval_tval(Pairs_i(j),2) = tstats.tstat; 
       
       %[~,p1,~,tstats1] = ttest2(avglen(HCs,Pairs_i(j)),avglen(SZs,Pairs_i(j))); % ttest on length
       [~,p1,~,tstats1] = ttest2(HClen,SZlen);
       len_pval_tval(Pairs_i(j),1) = p1;
       len_pval_tval(Pairs_i(j),2) = tstats1.tstat;
      %{
       [~,p2,~,tstats2] = ttest2(Freqs(HCs,Pairs_i(j)),Freqs(SZs,Pairs_i(j)));   % ttest on frequency
       freq_pval_tval(Pairs_i(j),1) = p2;
       freq_pval_tval(Pairs_i(j),2) = tstats2.tstat;
       %}
      
    end
    
    pval_tval{i,1} = mConn_pval_tval;   % Mean pair value 
    pval_tval{i,2} = maxConn_pval_tval; % Max pair value 
    pval_tval{i,3} = len_pval_tval;     % Length of shapelets 
    %pval_tval{i,4} = freq_pval_tval;    % Frequency of shapelets 
end
%%








