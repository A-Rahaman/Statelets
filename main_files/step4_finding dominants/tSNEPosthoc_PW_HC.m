% to get the dominant shapes in a given subject
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
dataDir ='/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/HC_PW_tSNE'; 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/probabilityDensity_HC.mat
%% Map the coordinates a 2D matrix and weight each of the cell with probability density
pMats = cell(314,1);
for i = 1:1081
    i
    Ps = pd{1,i};
    Ysub = load(fullfile(dataDir,['pair_',num2str(i),'.txt']));
    X = Ysub(:,1);
    Y = Ysub(:,2);
    pMat = zeros(1+ceil((max(X)-min(X))),1+ceil((max(Y)-min(Y))));
    %pMatfreq = zeros(1+ceil((max(X)-min(X))),1+ceil((max(Y)-min(Y))));
    for j = 1:length(X)
        x = 1+round(X(j)+(-1)*min(X));
        y = 1+round(Y(j)+(-1)*min(Y));
        pMat(x,y) = pMat(x,y)+Ps(j);
        %pMatfreq(x,y)
    end
    pMats{i} = pMat;
end
save(fullfile(outdir,'pMats_HC.mat'),'pMats');
%% Find the peaks
%[pks,locs] = findpeaks(bpMat,'MinPeakDistance',6);
%[pks,locs] = findpeaks(bpMat);
%[cent, cm] = FastPeakFind(pMat);
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
load(fullfile(outdir,'pMats_HC.mat'))
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/anlysisScripts/FastPeakFind'));
peaklocs = cell(size(pMats,1),1);

for i = 1:length(pMats)
pMat = pMats{i,1};
bpMat = imgaussfilt(pMat,4.5);
%bpMat = bpMat/max(bpMat(:));
[cent, cm] = FastPeakFind(bpMat);
%ii = imregionalmax(bpMat);
%[xloc,yloc,~] = find(ii);        % coordinates
[xloc,yloc] = find(cm==1);
locs = [xloc yloc];
peaklocs{i} = locs;
end
save(fullfile(outdir,'Peaklocs_HC.mat'),'peaklocs');
% imagesc(bpMat)
% hold on
% for fn=1:length(x)
%     text(y(fn),x(fn),'*','FontSize',25)
% end
% hold off
%mesh(0.01*cm+bpMat)
%mesh(ii+bpMat)
%% Take 10 closest neigbours and get the max of them
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/HC_PW_tSNE';
load(fullfile(outdir,'Peaklocs_HC.mat'))
load(fullfile(outdir,'probabilityDensity_HC.mat'))
for i = 1:length(peaklocs)
    i
    coord = peaklocs{i};
    Ysub = load(fullfile(dataDir,['pair_',num2str(i),'.txt']));
    P = pd{i};
    
    peaks_i_j = []; % Crucial   
    for j=1:length(coord)
        cx = coord(j,:);
        cx(1) = cx(1)+min(Ysub(:,1)); % map to tSNE plot coordiante x
        cx(2) = cx(2)+min(Ysub(:,2)); % map to tSNE plot coordiante y
        
        % Find the closest match ---------------------------------
        d = zeros(length(Ysub),1); % Crucial
        for k=1:length(Ysub)
            cy = Ysub(k,:);
            d(k) = norm(cx-cy);
        end
        [~,clId] = min(d);
        clCoord = Ysub(clId,:);  
        realPeaks(j) = clId;
        
        % --------------------------------------------------------
        % distance from that real point in tSNE to get 10 nearest
        dd = zeros(length(Ysub),1); % Crucial
        for k=1:length(Ysub)    
            cy = Ysub(k,:);
            dd(k) = norm(clCoord-cy);
        end
        [sortedDD,idx] = sort(dd);
        I =  idx(1:10);           
        PI = P(I); 
        % ---------------------------------------------------------
        % Find one with the optimal probability within that 10 
        [~,mPI] = max(PI);
        peak_j = I(mPI);
        peaks_i_j(j) = peak_j; 
        % ---------------------------------------------------------
    end
    peaksIJ{i} = peaks_i_j; 
    % Filter some of the peaks which are very close to each other. keep maxp always
%     peaks_i = peaksIJ{i};
%     for l = 1:length(peaks_i)
%         
%     end
    % --------------------------------------------------------------
end
save(fullfile(outdir,'pairWISEPeaks_HC.mat'),'peaksIJ');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/subWISEPeaks.mat','peaksIJ');
%% Getting the most dominant shapes from each pair. 
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
load(fullfile(outdir,'pairWISEPeaks_HC.mat'))
count = 1;
for h =1:size(peaksIJ,2)
    h
  load (fullfile(outdir,'HC_pairwise_shapelet',['pair_',num2str(h,'%04.f'),'.mat']));
  kk = peaksIJ{h};
  for pp = 1:length(kk)
      allPairshapes{count,1} = shapes(kk(pp)).real_length;%allshapes_s{kk(pp),1}; % for grabing the real length 
      allPairshapes{count,2} = shapes(kk(pp)).extrapolated;%allshapes_s{kk(pp),2}; % extrapolated
      allPairshapes{count,3} = shapes(kk(pp)).whichsubject;%allshapes_s{kk(pp),3}; % subject
      allPairshapes{count,4} = h;
      count = count+1;
  end
end
%save(fullfile(dirr,['pairWISEshapes_SZ' '.mat']),'pairWISEshapes');
save(fullfile(outdir,['allPairshapes_HC' '.mat']),'allPairshapes');

%% get EMD for all the shapes from all the pairs
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
load(fullfile(outdir,'allPairshapes_HC.mat'));
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/shapeletAlgorithmScripts'));
simscore = zeros(size(allPairshapes,1));
for l=1:size(allPairshapes,1)
allPairshapes_s{l} = allPairshapes{l,2};
end

for j = 1:size(allPairshapes_s,2)
  j
  x = allPairshapes_s{j};  
  for k = 1:size(allPairshapes_s,2)
    y = allPairshapes_s{k};
    x = x-min(x);
    x = x/sum(x);
    y = y-min(y);
    y = y/sum(y);
    simscore(j,k) = dEMD(x,y);
  end
end
save(fullfile(outdir,['allPairEMD_HC','.mat']),'simscore');
%{
%  ---------------This portion is for other analysis like frequency ---------------------------  
% get pairwise leader shape. One form each pair  
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load pairWISEPeaks_HC.mat
load probabilityDensity_HC.mat
cd('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/HC_pairwise_shapelet')
for i = 1:1081
    i
    load(['pair_' num2str(i,'%04.f') '.mat']);
    peaks = peaksIJ{i};
    pds = pd{i};
    highpd = pds(peaks);
    idx = peaks(highpd==max(highpd));
    leaderFreq_HC(i) = max(highpd);
    pairwiseleader_HC{i} = shapes(idx).real_length;
end

%% Peaks count group Wise 
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load pairWISEPeaks_HC.mat
counts_HC = zeros(1081,1);
for i=1:length(peaksIJ)
    counts_HC(i) = length(peaksIJ{i});
end
%% plotting counts
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/mrahaman/d_NBiC/Data/idxUpper.mat % Upper 1081 idx (IND1)
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'hot';
%cm = 'jet';
% FDR correction of tvals then keeping only the significant tvals    
%% HC
data_HC = zeros(47);    
%datt1(idx) = log(HC_i+eps); 
data_HC(idx) = (counts_HC);  
datt21 = rot90(fliplr(data_HC));
Fdatt_HC = zeros(47);
%Fdatt1(idx) = log(HC_i+eps);
Fdatt_HC(idx) = (counts_HC);

FNC1 = datt21+Fdatt_HC; 
%FNC = datt2+datt;
%pvals = tvalMat{db3,1}; % pval
%CLIM = [-max(abs(HC_i(:))) max(abs(HC_i(:)))];
CLIM = [0 max((Fdatt_HC(:)))];
T = 'HC Counts';
%fName = ['t-val',num2str(db3) '.fig'];
%[F,A,C,I] = plot_FNC(rot90(fliplr(datt)), CLIM, LABEL, RSN_I,'',T,MOD,cm,pvals,0); % db3 as last argument for normal plot 
[F1,A1,C1,I1] = plot_FNC(FNC1, CLIM, LABEL, RSN_I,'',T,MOD,cm,Fdatt_HC,0); 

%}
%%







