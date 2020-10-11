%% run fft on the pairwise leader shape for 9 different frequency bands 
clear
clc
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/pairWISELeaders';
%load /Users/mrahaman1/Documents/StateLets/metaData/Wavelets_314X1081.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/SM_SUB_wise_leaderShapelets_uLen.mat
%Wavelets = SUB_wise_shapelets_uLen;
for i = 1:314
    load(fullfile(dataDir,['sub_',num2str(i,'%03.f'),'.mat']));
    for j=1:size(pairWISELeadershapes,1)
        w = pairWISELeadershapes{j} - mean(pairWISELeadershapes{j}); % removing the mean 
        Freqs(i,j,:) = fft(w,9);  % fft on the leader shape of ith subject and jth pair 
    end
end
save('/Users/mrahaman1/Documents/Statelet_V2/fResults/FFT/fft_frequencies.mat','Freqs');
%%
clear
clc
%load /Users/mrahaman1/Documents/StateLets/metaData/Statelets_SM_k5.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/SZsubjectIdx.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/HCsubjectIdx.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/frequencies_fft_rmovingmean.mat  
load /Users/mrahaman1/Documents/Statelet_V2/fResults/FFT/fft_frequencies.mat % frequencies 

HCs = CtIdx;
SZs = SZIdx; 
for i = 1:size(Freqs,3)
    freq_pval_tval = zeros(1081,2);
    for j = 1:size(Freqs,2) 
       %[~,p,~,tstats] = ttest2(abs(Freqs(HCs,j,i)),abs(Freqs(SZs,j,i)));   
       [~,p,~,tstats] = ttest2(abs(Freqs(HCs,j,i)),abs(Freqs(SZs,j,i)));    % ttest on fft freqs
       freq_pval_tval(j,1) = p;
       freq_pval_tval(j,2) = tstats.tstat; 
    end 
    pval_tval{i,1} = freq_pval_tval;  % Frequency of shapelets 
end
save('/Users/mrahaman1/Documents/Statelet_V2/fResults/FFT/pval-tval_fft.mat','pval_tval');
%% Plotting ttest results 
clear
clc
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/StateLets/metaData/pval-tval_fft_rm.mat
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'jet';
figDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/FFT/Figures';
%% FDR correction of tvals then keeping only the significant tvals    
FDRtvals = [];

%figure(60);
for i = 1:size(pval_tval,1)
PT = pval_tval{i,1};
Lpvals = PT(:,1);
Ltvals = PT(:,2);
%Ltvals = real(PT(:,2));

temptvals = zeros(length(Ltvals),1);
Correctedpval = myFDR(Lpvals);
tempid = find(Correctedpval<0.05);   % FDR correction 
%tempid = find(Lpvals<0.002);        % thresholding 
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
%subplot(1,9,i)
[F,A,C,I] = plot_FNC(FNC, CLIM, LABEL, RSN_I,'',T,MOD,cm,Lpvals,0); 
fName = ['t-val_fft_',num2str(i),'.fig'];
set(gcf, 'Position',  [100, 100, 400, 300])% size of 400 X 300 pixels 
saveas(F,fullfile(figDir,fName));

close
end
%%