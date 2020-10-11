clear
clc
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
%load /Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/REP_pval-tval.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/SM_pval-tval.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/pval-tval_generalized.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/CORR/pval-tval_CoRR_K5.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/CORR/pval-tval_CoRR_K7.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/eucledian/pval-tval_eu_K5.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/eucledian/pval-tval_eu_K7.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/pval-tval_SM_k7.mat
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'jet';
% co-occurences
load /Users/mrahaman1/Documents/StateLets/metaData/pairwiseMeanHC_co.mat
load /Users/mrahaman1/Documents/StateLets/metaData/pairwiseMeanSZ_co.mat
% FDR correction of tvals then keeping only the significant tvals    
FDRtvals = [];
%% Plotting tvals for meanLen and Frequency  
figDir = '/Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/K7';

iz = 1:size(pval_tval,1);

for i =1:length(iz)
    
% Average length of shapelets 
%Lpvals = L_pval_tval(:,1);
%Ltvals = L_pval_tval(:,2);

% For Frequencey
for v = 1:1size(pval_tval,2)
PT = pval_tval{iz(i),v};
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
fName = ['t-val_K',num2str(size(pval_tval,1)),'_',num2str(iz(i)),'-',num2str(v),'.fig'];
%[F,A,C,I] = plot_FNC(rot90(fliplr(datt)), CLIM, LABEL, RSN_I,'',T,MOD,cm,pvals,0); % db3 as last argument for normal plot 
[F,A,C,I] = plot_FNC(FNC, CLIM, LABEL, RSN_I,'',T,MOD,cm,Lpvals,0); 
%saveas(F,fullfile(figDir,fName));
%close
end
end

%%