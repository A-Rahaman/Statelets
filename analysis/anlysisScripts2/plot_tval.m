clear
clc
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/mrahaman/d_NBiC/Data/idxUpper.mat % Upper 1081 idx (IND1)
%load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/pval-tval_Freq.mat % Frequency 
%load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/pval-tval_AvgLen.mat % Mean shapelet length
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/pval_tval.mat % all pval-tval meanconn-maxconn-len-freq 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/pval-tval_shapesonlyk5.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/pval_tval_k10.mat % all pval-tval meanconn-maxconn-len-freq 
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'jet';
% FDR correction of tvals then keeping only the significant tvals    
FDRtvals = [];
figDir ='/Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/ttestFigureshapesonly';
%% Plotting tvals for meanLen and Frequency  
%figDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/ttestFigures';
%figDir ='/Users/mrahaman1/Documents/Statelet_V2/fResults/KMeans/ttestFigures/pngs';

%N_inits = ['meanConn','maxConn','len','freq'];
for i =1:size(pval_tval,2) % iterate over properties
for j=1:size(pval_tval,1)  % iterate over statelets
PT = pval_tval{j,i};
Lpvals = PT(:,1);
Ltvals = PT(:,2);
px = find(isnan(Ltvals)==1);
Ltvals(px) = 0;
fName = ['tval_',num2str(i),'_','statelet_',num2str(j),'.png']; 

%temptvals = zeros(length(Ltvals),1);
%Correctedpval = myFDR(Lpvals);
%tempid = find(Correctedpval<0.05);
%tempid = find(Lpvals<0.05);        % thresholding 
%temptvals(tempid) = Ltvals(tempid);
%FDRtvals = temptvals;

% tval matrix printing 
datt = zeros(47);    
datt(idx) = Ltvals;  % tval
datt2 = rot90(fliplr(datt));
Fdatt = zeros(47);
Fdatt(idx) = Ltvals;%FDRtvals;

FNC = datt2+Fdatt; % Adding two matrices 47X47. One have lower traingle all incorrected values and one have upper traingle with FDR corrected.  
%FNC = datt2+datt;
%pvals = tvalMat{db3,1}; % pval
%CLIM = [-max(abs(Ltvals(:))) max(abs(Ltvals(:)))];

CLIM = [min((Ltvals(:))) max((Ltvals(:)))];

T = 't-values';
%fName = ['t-val',num2str(db3) '.fig'];
%[F,A,C,I] = plot_FNC(rot90(fliplr(datt)), CLIM, LABEL, RSN_I,'',T,MOD,cm,pvals,0); % db3 as last argument for normal plot 
[F,A,C,I] = plot_FNC(FNC, CLIM, LABEL, RSN_I,'',T,MOD,cm,Lpvals,0); 
set(gcf, 'Position',  [100, 100, 400, 300])
saveas(F,fullfile(figDir,fName));
close
end
end
%%