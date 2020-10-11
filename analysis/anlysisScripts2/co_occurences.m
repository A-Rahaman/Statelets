%% co-occurences of the shapelets. Lets hammer the temporal information of extracted shapelets 
% For the cooccurrences, for a given subject, a given pair, get each occurrence of shapelet and go down across all other pairs 
% using that time information. Count how many other matched - I mean occurred as well probably with a different shape.  
% Fr = #time shapelet of pair j (1 to 1080) occurred/ #time shapelet of pair i occurred  
% Get 1081 X 1080 matrix for each subject  
% Analyze it for HC and SZ separately. Get one mean occurrences across all the HCs and SZs  
% For each group you get a mean occurrence matrix of 1081 X 1080 
% Print that as you print tvals.  
% If you get higher value there that means it is happening across the subjects. So consistent  
% Hopefully, get some modularity in it.      
% t-test on occurrences for group differences     
% by //munna//
% Dated 11 July 2019
%% Load the data
clear
clc
load /Users/mrahaman1/Documents/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/Info_HCs.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/Info_SZs.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/Candidates_10-50.mat
co_occurences_SZ = zeros(size(Info_SZs,1),1081,1080);
co_occurences_HC = zeros(size(Info_HCs,1),1081,1080);

 %% SZ
 
 for i = 1:size(Info_SZs,1)
     Info_i  = Info_SZs{i,1}.repWavelet_Details;
     for j = 1:size(Info_i,1)
        cands_ij = candidates{1,Info_i{j,1}};                              % Candidates for pair j and subject i
        occur_ij = Info_i{j,3};                                            % Occurences for pair j of subject i 
        pairs = 1:1081;
        k_pairs  = pairs(pairs~=j); 
        for k = 1:size(k_pairs,2)
            cands_ik = candidates{1,Info_i{k_pairs(k),1}};                          % Candidates for pair k and subject i
            occur_ik = Info_i{k_pairs(k),3};                               % Occurences for pair k of subject i 
            cooc_count = 0;                                                % Co-occurence count init
            for oc_j = 1:length(occur_ij)
                for oc_k =1:length(occur_ik)
                    if(intersect(cands_ij{occur_ij(oc_j)}, cands_ik{occur_ik(oc_k)}))
                        cooc_count = cooc_count+1;                         % Co-occurence count incremented
                        break;
                    end
                end
            end
            co_occurences_SZ(i,j,k) = cooc_count/length(occur_ij);         % Co-occurences betn pair j and k of sub i
        end
     end
 end
 %% HC
 
 for i = 1:size(Info_HCs,1)
     Info_i  = Info_HCs{i,1}.repWavelet_Details;
     for j = 1:size(Info_i,1)
        cands_ij = candidates{1,Info_i{j,1}};                              % Candidates for pair j and subject i
        occur_ij = Info_i{j,3};                                            % Occurences for pair j of subject i 
        pairs = 1:1081;
        k_pairs  = pairs(pairs~=j); 
        for k = 1:size(k_pairs,2)
            cands_ik = candidates{1,Info_i{k_pairs(k),1}};                 % Candidates for pair k and subject i
            occur_ik = Info_i{k_pairs(k),3};                               % Occurences for pair k of subject i 
            cooc_count = 0;                                                % Co-occurence count init
            for oc_j = 1:length(occur_ij)
                for oc_k =1:length(occur_ik)
                    if(intersect(cands_ij{occur_ij(oc_j)}, cands_ik{occur_ik(oc_k)}))
                        cooc_count = cooc_count+1;                         % Co-occurence count incremented
                        break;
                    end
                end
            end
            co_occurences_HC(i,j,k) = cooc_count/length(occur_ij);         % Co-occurences betn pair j and k of sub i
        end
     end
 end
 %% Analysis
 clear
 clc
 load /Users/mrahaman1/Documents/StateLets/metaData/co_occurences_HC
 load /Users/mrahaman1/Documents/StateLets/metaData/co_occurences_SZ
 
 meanCooccurenceSZ = squeeze(mean(co_occurences_SZ)); % across all SZ
 meanCooccurenceHC = squeeze(mean(co_occurences_HC)); % across all HC
 pairwiseMeanSZ = mean(meanCooccurenceSZ,2); % across all the SZ and 1 pair all others
 pairwiseMeanHC = mean(meanCooccurenceHC,2); % across all the SZ and 1 pair all others
 subjectwise_pairwiseMeanSZ = mean(co_occurences_SZ,3);
 subjectwise_pairwiseMeanHC = mean(co_occurences_HC,3);
%%
figure()
subplot(1,2,1);
im1 = imagesc(meanCooccurenceHC);
 %CLIM = [-max(abs(meanCooccurenceHC(:))) max(abs(meanCooccurenceHC(:)))];
 %set(gca, 'clim', CLIM); %axis ij;
% colormap('jet');
 %set(im1,'Title','HC subject Co-occurence');
subplot(1,2,2);
im2 = imagesc(meanCooccurenceSZ); %colormap('jet');
 %CLIM = [-max(abs(meanCooccurenceSZ(:))) max(abs(meanCooccurenceSZ(:)))];
 %set(gca, 'clim', CLIM); %axis ij;
%% mat plot
clear
clc
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'jet';
% co-occurences
load /Users/mrahaman1/Documents/StateLets/metaData/pairwiseMeanHC_co.mat
load /Users/mrahaman1/Documents/StateLets/metaData/pairwiseMeanSZ_co.mat

%% Plot
% extra code segment for aligning the data dimension 
pval_tval(:,1) = pairwiseMeanHC';
pval_tval(:,2) = pairwiseMeanSZ';
%------------------------------------
for i=1:size(pval_tval,2) 
    
Ltvals = pval_tval(:,i);
Lpvals = zeros(1081,1);
%Correctedpval = myFDR(Lpvals);
tempid = find(Ltvals>=0.93);
temptvals(tempid) = Ltvals(tempid);
FDRtvals = temptvals;
temptvals = zeros(length(Ltvals),1);
datt = zeros(47);    
datt(idx) = Ltvals;  % tval
datt2 = rot90(fliplr(datt));
Fdatt = zeros(47);
Fdatt(idx) = FDRtvals;
FNC = datt2+Fdatt;
CLIM = [-max(abs(datt2(:))) max(abs(datt2(:)))];
T = 't-values';
[F,A,C,I] = plot_FNC(datt2, CLIM, LABEL, RSN_I,'',T,MOD,cm,Lpvals,0); 
%fName = ['t-val_K',num2str(size(pval_tval,1)),'_',num2str(iz(i)),'-',num2str(v),'.fig'];
%saveas(F,fullfile(figDir,fName));
%close
end
%% tvals
clear
clc
load /Users/mrahaman1/Documents/StateLets/metaData/subpair_HC_co_ttest.mat
load /Users/mrahaman1/Documents/StateLets/metaData/subpair_SZ_co_ttest.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'jet';
CoOccur_pval_tval = zeros(1081,2);
for i =1:size(subjectwise_pairwiseMeanHC,2)
    [~,p,~,tstats] = ttest2(subjectwise_pairwiseMeanHC(:,i),subjectwise_pairwiseMeanSZ(:,i));                                  
    CoOccur_pval_tval(i,1) = p;
    CoOccur_pval_tval(i,2) = tstats.tstat; 
end
%% Plot tvals
PT =   CoOccur_pval_tval;
Lpvals = PT(:,1);
Ltvals = PT(:,2);

temptvals = zeros(length(Ltvals),1);
Correctedpval = myFDR(Lpvals);
tempid = find(Correctedpval<0.05);
temptvals(tempid) = Ltvals(tempid);
FDRtvals = temptvals;

% tval matrix printing 
datt = zeros(47);    
datt(idx) = Ltvals;  % tval
datt2 = rot90(fliplr(datt));
Fdatt = zeros(47);
Fdatt(idx) = FDRtvals;
FNC = datt2+Fdatt; % Adding two matrices 47X47. One have lower traingle all incorrected values and one have upper traingle with FDR corrected.  
%pvals = tvalMat{db3,1}; % pval
CLIM = [-max(abs(Ltvals(:))) max(abs(Ltvals(:)))];
T = 't-values';
%fName = ['t-val_K',num2str(size(CoOccur_pval_tval,1)),'_',num2str(iz(i)),'-',num2str(v),'.fig'];
%[F,A,C,I] = plot_FNC(rot90(fliplr(datt)), CLIM, LABEL, RSN_I,'',T,MOD,cm,pvals,0); % db3 as last argument for normal plot 
[F,A,C,I] = plot_FNC(FNC, CLIM, LABEL, RSN_I,'',T,MOD,cm,Lpvals,0); 
%saveas(F,fullfile(figDir,fName));
%close
%%
 
 
 
 
 
 
 
 
 
 
 
 
 