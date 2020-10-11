% call KDE
% to get the dominant shapes in a given pair
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/SZ_PW_tSNE';
load /Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/probabilityDensity_SZ.mat
%% Map the coordinates a 2D matrix and weight each of the cell with probability density
pMats = cell(314,1);
for i = 1:1081
    i
    Ps = pd{1,i};
    %Ysub = load(fullfile(dataDir,['pair_',num2str(i),'.txt']));
    Ysub = load('tSNE_allpairshape_SZ.txt');
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
save(fullfile(outdir,'pMats_SZ.mat'),'pMats');
%% Find the peaks
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
load(fullfile(outdir,'pMats_SZ.mat'))
addpath(genpath('/data/mialab/users/mrahaman/statelets/scripts/shapeletAlgorithmScripts/FastPeakFind'));
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
save(fullfile(outdir,'Peaklocs_SZ.mat'),'peaklocs');
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
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/SZ_PW_tSNE';
load(fullfile(outdir,'Peaklocs_SZ.mat'))
load(fullfile(outdir,'probabilityDensity_SZ.mat'))
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
save(fullfile(outdir,'pairWISEPeaks_SZ.mat'),'peaksIJ');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/subWISEPeaks.mat','peaksIJ');
%% Getting the most dominant shapes 
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
load(fullfile(outdir,'pairWISEPeaks_SZ.mat'))
count = 1;
for h =1:size(peaksIJ,2)
    h
  load (fullfile(outdir,'SZ_pairwise_shapelet',['pair_',num2str(h,'%04.f'),'.mat']));
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
save(fullfile(outdir,['allPairshapes_SZ' '.mat']),'allPairshapes');
%% Leaders
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load pairWISEPeaks_SZ.mat
load probabilityDensity_SZ.mat
cd('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/SZ_pairwise_shapelet')
for i = 1:1081
    i
    load(['pair_' num2str(i,'%04.f') '.mat']);
    peaks = peaksIJ{i};
    pds = pd{i};
    highpd = pds(peaks);
    idx = peaks(highpd==max(highpd));
    leaderFreq_SZ(i) = max(highpd);
    pairwiseleader_SZ{i} = shapes(idx).real_length;
end
%% get EMD for all the shapes from all the pairs
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
load(fullfile(outdir,'allPairshapes_SZ.mat'));
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
    % already normalized
    %x = x-min(x);
    %x = x/sum(x);
    %y = y-min(y);
    %y = y/sum(y);
    simscore(j,k) = dEMD(x,y);
  end
end
save(fullfile(outdir,['allPairEMD_SZ','.mat']),'simscore');

%% Peak counts SZ 
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load pairWISEPeaks_SZ.mat
counts_SZ = zeros(1081,1);
for i=1:length(peaksIJ)
    counts_SZ(i) = length(peaksIJ{i});
end
%% Combined HC&SZ pairwise leaders plotting 
clear
clc
cd('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load pairwiseleaders_HC.mat
load pairwiseleaders_SZ.mat
figure()
for k=1:1081
    subplot(33,33,k)
    tempSZ = pairwiseleader_SZ{k};
    plot(tempSZ,'r');
    hold on
    tempHC = pairwiseleader_HC{k};
    plot(tempHC,'b');
    maxl = max(length(tempSZ),length(tempHC));
    plot(zeros(maxl,1),'k')
    box off
    axis tight
    axis off 
end

%% SZ&HC pair coutns plotting
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis'))
clear
clc
%load Peakcounts_HC.mat
%load Peakcounts_SZ.mat
load HC_pairwise_probability.mat
load SZ_pairwise_probability.mat

[SZs,sortedSZidx] = sort(SZ_i,'descend');
HC_as_SZ = HC_i(sortedSZidx);
[HCs,sortedHCidx] = sort(HC_i,'descend');
SZ_as_HC = SZ_i(sortedHCidx);

SZ_HC = [SZs HC_as_SZ];
HC_SZ = [HCs SZ_as_HC];
figure()
subplot(2,1,1)
bar(SZ_HC)
axis tight
subplot(2,1,2)
bar(HC_SZ)
axis tight
%% SZ-HC pairwise probability %% Sum the probabiltity of all shapes coming from a given pair p
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load allPairshapesPD_HC.mat
load allPairshapes_HC.mat
% HC
for i =1:size(allPairshapes,1)
HC_shape_to_pair_mapping(i,1) = allPairshapes{i,4}; 
end
HC_pairwise_maxPD = zeros(1081,1);
HC_pairwise_PD = zeros(1081,1);
HC_pairwise_dom = zeros(1081,1);
for i=1:length(HC_pairwise_PD)
    shapes_i = find(HC_shape_to_pair_mapping==i);
    [HC_pairwise_maxPD(i),midx] = max(p(shapes_i));
    HC_pairwise_dom(i) = max(shapes_i(midx));
    HC_pairwise_PD(i) = sum(p(shapes_i));
end

% SZ
clear
clc
load allPairshapesPD_SZ.mat
load allPairshapes_SZ.mat

for i =1:size(allPairshapes,1)
SZ_shape_to_pair_mapping(i,1) = allPairshapes{i,4}; 
end
SZ_pairwise_maxPD = zeros(1081,1);
SZ_pairwise_PD = zeros(1081,1);
SZ_pairwise_dom = zeros(1081,1);
for i=1:length(SZ_pairwise_PD)
    shapes_i = find(SZ_shape_to_pair_mapping==i);
    [SZ_pairwise_maxPD(i),midx] = max(p(shapes_i));
    SZ_pairwise_dom(i) = max(shapes_i(midx));
    SZ_pairwise_PD(i) = sum(p(shapes_i));
end

%% Sorting SZ and HC probability density values w.r.t. others and plot
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load HC_pairwise_probability.mat
load SZ_pairwise_probability.mat

SZ_i = SZ_pairwise_PD;
HC_i = HC_pairwise_PD;
% sort SZ
[~,sortedSZ_idx] = sort(SZ_i,'descend');
sortedSZs = SZ_i(sortedSZ_idx);        % SZ sorted hight to low probability density   
sortedHCs = HC_i(sortedSZ_idx);        % HC sorted accrodingly
sortedSZ_HC = [sortedSZs sortedHCs];

% sort HC
[~,sortedHC_idx] = sort(HC_i,'descend');
sortedHCs1 = HC_i(sortedHC_idx);        % HC sorted hight to low probability density   
sortedSZs1 = SZ_i(sortedHC_idx);        % SZ sorted accrodingly 
sortedHC_SZ = [sortedSZs1 sortedHCs1];  

figure()
subplot(1,2,1)
h = barh(sortedSZ_HC);
axis tight
h(1).FaceColor = 'r'; % color
h(2).FaceColor = 'b';%[0.3 0.78 0.5];% rgb
subplot(1,2,2)
h1 = barh(sortedHC_SZ);
h1(1).FaceColor = 'r'; % color
h1(2).FaceColor = 'b';%[0.3 0.78 0.5];% rgb
axis tight

% Plotting top 10 for both combinations HC sorted and SZ sorted 
%% Create the colormap
%% gradient map
clear
clc
Bar_range = [-14 14];
colmap = colormap(jet);
ini_idx = 20; % small:edge close to black;
zero_modified_idx = 2;
ratio = (Bar_range(2) - 0)/(0 - Bar_range(1));
num_hot = 100;
num_col = round((num_hot-zero_modified_idx)/ratio);
% hot map
colormap(rand(num_hot,3))
hot_map = colormap(hot);
hot_map_use = hot_map(end-zero_modified_idx:-1:(1+ini_idx),:);
% col map
colormap(rand(num_col,3))
hot_map = colormap(hot);
col_map_use = zeros(size(hot_map,1)-round(ini_idx/ratio),size(hot_map,2));
col_map_use(:,1) = hot_map((1+round(ini_idx/ratio):end),3);
col_map_use(:,2) = hot_map((1+round(ini_idx/ratio):end),2);
col_map_use(:,3) = hot_map((1+round(ini_idx/ratio):end),1);

colmap = [col_map_use;hot_map_use];
colormap(colmap);

%% heat map; just plot the pairwise probability density 

% The probability density by a pair is the sum of PD by all the shapes coming from that pair
% This has been done for group wise dynamics HC/SZ seperate (HC: 9285 and SZ: 8706 shapes) 
clear
clc
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/mrahaman/d_NBiC/Data/idxUpper.mat % Upper 1081 idx (IND1)
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
%cm = colmap;
cm = 'jet';
%cm = 'parula';
% FDR correction of tvals then keeping only the significant tvals    
%FDRtvals = [];
load HC_pairwise_probability.mat
load SZ_pairwise_probability.mat
%load Peakcounts_HC.mat
%load Peakcounts_SZ.mat
%% HC
HC_i = HC_pairwise_PD;
%HC_i = counts_HC;
datt1 = zeros(47);    
%datt1(idx) = log(HC_i+eps); 
datt1(idx) = (HC_i);  
datt21 = rot90(fliplr(datt1));
Fdatt1 = zeros(47);
%Fdatt1(idx) = log(HC_i+eps);
Fdatt1(idx) = (HC_i);
FNC1 = datt21+Fdatt1; 
%CLIM = [-max(abs(HC_i(:))) max(abs(HC_i(:)))];
%CLIM = [0 0.5656]; % generalized for both SZ and HC
%CLIM = [-14 14]; % for the counts 
%CLIM = [min((Fdatt1(:))) max((Fdatt1(:)))];
T1 = 'HC';
%[F1,A1,C1,I1] = plot_FNC(FNC1, CLIM, LABEL, RSN_I,'',T1,MOD,cm,Fdatt1,0); 
%% SZ
SZ_i = SZ_pairwise_PD;
%SZ_i = counts_SZ;
datt = zeros(47);    
%datt(idx) = log(SZ_i+eps);  % tval
datt(idx) = (SZ_i);  % tval
datt2 = rot90(fliplr(datt));
Fdatt = zeros(47);
%Fdatt(idx) = log(SZ_i+eps);
Fdatt(idx) = (SZ_i);
FNC2 = datt2+Fdatt; 
%FNC = datt2+datt;
%CLIM = [-max(abs(HC_i(:))) max(abs(HC_i(:)))];
%CLIM = [min((Fdatt(:))) max((Fdatt(:)))];
%CLIM = [-14 14]; % for the counts 

%CLIM = [0 0.5656];
T2 = 'SZ';
%fName = ['t-val',num2str(db3) '.fig'];
%[F2,A2,C2,I2] = plot_FNC(FNC2, CLIM, LABEL, RSN_I,'',T2,MOD,cm,log(SZ_i+eps),0);

%% differeces 
FF=zeros(47);
%FFF =zeros(47);
FF(idx) = HC_i - SZ_i;
FF = rot90(fliplr(FF));
pal1 = zeros(47);
pal1(idx) = HC_i - SZ_i;
FNCC = FF+pal1;

%FFF(idx) = SZ_i;
%pal2 = rot90(fliplr(FFF));
%FFF = FFF +pal2;
%F1_F2   = FF-FFF; 
%clear FNCC
%FNCC = F1_F2;
%FNCC = FF;
%FNCC = log(F1_F2+eps);
%CLIM = [-max(abs(FNCC(:))) max(abs(FNCC(:)))];

%CLIM = [-14 14]; % for the counts 
M1 = max(max(abs(FNC1(:))),max(abs(FNC2(:))));
M2 = max(M1,max(abs(FNCC(:))));
CLIM = [-M2 M2]; % for the counts 
T3 = 'HC-SZ';
%fName = ['t-val',num2str(db3) '.fig'];
%[F,A,C,I] = plot_FNC(rot90(fliplr(datt)), CLIM, LABEL, RSN_I,'',T,MOD,cm,pvals,0); % db3 as last argument for normal plot 
[F,A,C,I] = plot_FNC(FNCC, CLIM, LABEL, RSN_I,'',T3,MOD,cm,FNCC,0);
[F1,A1,C1,I1] = plot_FNC(FNC1, CLIM, LABEL, RSN_I,'',T1,MOD,cm,Fdatt1,0); 
[F2,A2,C2,I2] = plot_FNC(FNC2, CLIM, LABEL, RSN_I,'',T2,MOD,cm,log(SZ_i+eps),0);
%% clear
clear
clc
cd('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
% sum of Probability across all the shapes of that pair 
load HC_pairwise_probability.mat
load SZ_pairwise_probability.mat
% max Probability across all the shapes of that pair 
%load HC_pairwise_maxPD.mat 
%load SZ_pairwise_maxPD.mat
%load pairwiseleaders_HC.mat
%load pairwiseleaders_SZ.mat
load HC_pairwise_dominants.mat
load SZ_pairwise_dominants.mat

SZs = load('allPairshapes_SZ.mat');
SZshapes = SZs.allPairshapes;
HCs = load('allPairshapes_HC.mat');
HCshapes = HCs.allPairshapes;


[SZpal,SZidx] = sort(SZ_pairwise_PD,'descend');
top_10_SZPairs = SZidx(1:10);


[HCpal,HCidx] = sort(HC_pairwise_PD,'descend');
top_10_HCPairs = HCidx(1:10);



figure()
for k=1:length(top_10_SZPairs)
    subplot(2,5,k)
    tempSZ = SZshapes{SZ_pairwise_dom(top_10_SZPairs(k)),1};
    plot(tempSZ,'r','LineWidth',2);
    hold on
    %tempHC = HCshapes{HC_pairwise_dom(top_10_HCPairs(k)),1};
    tempHC = HCshapes{HC_pairwise_dom(top_10_SZPairs(k)),1}; % just take the same pair from HC 
    plot(tempHC,'b','LineWidth',2);
    maxl = max(length(tempSZ),length(tempHC));
    plot(zeros(maxl,1),'k','LineWidth',2)
    title(['Pair ' num2str(top_10_SZPairs(k))])
    set(gca,'fontsize', 14);
    box off
    axis tight
    axis off 
end

figure()
for k=1:length(top_10_HCPairs)
    subplot(2,5,k)
    
    tempSZ = SZshapes{SZ_pairwise_dom(top_10_HCPairs(k)),1};
    plot(tempSZ,'r','LineWidth',2);
    hold on
    tempHC = HCshapes{HC_pairwise_dom(top_10_HCPairs(k)),1};
    %tempHC = HCshapes{HC_pairwise_dom(top_10_SZPairs(k)),1}; % just take the same pair from HC 
    plot(tempHC,'b','LineWidth',2);
    maxl = max(length(tempSZ),length(tempHC));
    plot(zeros(maxl,1),'k','LineWidth',2)
    title(['Pair ' num2str(top_10_HCPairs(k))])
    set(gca,'fontsize', 14);
    box off
    axis tight
    axis off 
end
%%
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load HC_pairwise_dominants.mat
load SZ_pairwise_dominants.mat
SZs = load('allPairshapes_SZ.mat');
SZshapes = SZs.allPairshapes;
HCs = load('allPairshapes_HC.mat');
HCshapes = HCs.allPairshapes;
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/shapeletAlgorithmScripts'))

pairwise_group_distance = zeros(1081,1);
for i =1:length(HC_pairwise_dom)
    i
    x  = SZshapes{SZ_pairwise_dom(i),2};
    y = HCshapes{HC_pairwise_dom(i),2};
    x = x-min(x);
    x = x/sum(x);
    y = y-min(y);
    y = y/sum(y);
    pairwise_group_distance(i) = dEMD(x,y);
end

load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/mrahaman/d_NBiC/Data/idxUpper.mat % Upper 1081 idx (IND1)
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'hot';  
FF=zeros(47);
FF(idx) = pairwise_group_distance;
FF = rot90(fliplr(FF));
pal1 = zeros(47);
pal1(idx) = pairwise_group_distance;
FNCC = FF+pal1;
CLIM = [-max(abs(FNCC(:))) max(abs(FNCC(:)))];
T = 'HC-SZ connectivity difference';
[F,A,C,I] = plot_FNC(FNCC, CLIM, LABEL, RSN_I,'',T,MOD,cm,FNCC,0);
%% The peaks from SZ and HC
clear
clc
load('SZ_Peaks.mat');
SZP = peaks_i_j;
clear peaks_i_j
load('HC_Peaks.mat');
HCP = peaks_i_j;
load('allPairshapes_SZ.mat');
SZshapes = allPairshapes;
load('allPairshapes_HC.mat');

figure()
for i =1:size(HCP,2)
    subplot(6,8,i)
    
    tempSZ = SZshapes{SZP(i),1};
    plot(tempSZ,'r','LineWidth',2);
    hold on
    %tempHC = HCshapes{HC_pairwise_dom(top_10_HCPairs(k)),1};
    tempHC = allPairshapes{HCP(i),1}; % just take the same pair from HC 
    plot(tempHC,'b','LineWidth',2);
   % maxl = max(length(tempSZ),length(tempHC));
    maxl = length(tempHC);
    hold on
    plot(zeros(maxl,1),'k','LineWidth',2)
    %title(['Pair ' num2str(top_10_HCPairs(k))])
    set(gca,'fontsize', 14);
    box off
    axis tight
    axis off 
end

%% Figure size 
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 4 10]);

%% Pair 417 analysis
clear
clc
load('allPairshapes_SZ.mat');
SZshapes = allPairshapes;
load('allPairshapes_HC.mat');
HCshapes = allPairshapes;
for i =1:size(HCshapes,1)
   HCPairs(i) = HCshapes{i,4};
   HCPshapes{i} = HCshapes{i,2};
end

for i =1:size(SZshapes,1)
   SZPairs(i) = SZshapes{i,4};
   SZPshapes{i} = SZshapes{i,2};
end
SZP417 = find(SZPairs==417);
HCP417 = find(HCPairs==417);
figure()

%SZ

for i=1:length(SZP417)
    subplot(3,4,i)
    plot(SZPshapes{SZP417(i)},'r','LineWidth',2);
    hold on
    plot(zeros(length(SZPshapes{SZP417(i)}),1),'k','LineWidth',2)
    set(gca,'fontsize', 14);
    box off
    axis tight
    axis off 
end

% HC
figure()
for i=1:length(HCP417)
    subplot(4,4,i)
    plot(HCPshapes{HCP417(i)},'b','LineWidth',2);
    hold on
    plot(zeros(length(HCPshapes{HCP417(i)}),1),'k','LineWidth',2)
    set(gca,'fontsize', 14);
    box off
    axis tight
    axis off 
end
%% Pairwise mxCONN testing 
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load HC_pairwise_dominants.mat
load SZ_pairwise_dominants.mat
SZs = load('allPairshapes_SZ.mat');
SZshapes = SZs.allPairshapes;
HCs = load('allPairshapes_HC.mat');
HCshapes = HCs.allPairshapes;
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/shapeletAlgorithmScripts'))
load HC_pairwise_dominants.mat
load SZ_pairwise_dominants.mat
SZCoNN = zeros(1081,1);
HCCoNN =zeros(1081,1);
for i =1:size(SZshapes,1)
    SZrealshapes{i} = SZshapes{i,1};
    SZextshapes{i} = SZshapes{i,2};
    SZshapespair(i) = SZshapes{i,4};
end
for i =1:size(HCshapes,1)
    HCrealshapes{i} = HCshapes{i,1};
    HCextshapes{i} = HCshapes{i,2};
    HCshapespair(i) = HCshapes{i,4};
end

for i =1:1081
    %SZshape_i = SZshapes{SZ_pairwise_dom(i),1};
    SZshapes_i=find(SZshapespair==i);
    allSZ_i = zeros(length(SZshapes_i),50);
    for j=1:length(SZshapes_i)
        allSZ_i(j,:) = SZextshapes{SZshapes_i(j)};
    end
    tempSZ = mean((allSZ_i));
    SZCoNN(i) = max(abs(tempSZ));
    
    
    %HCshape_i = HCshapes{HC_pairwise_dom(i),1};
    HCshapes_i=find(HCshapespair==i);
    allHC_i = zeros(length(HCshapes_i),50);
    for j=1:length(HCshapes_i)
        allHC_i(j,:) = HCextshapes{HCshapes_i(j)};
    end
    tempHC = mean((allHC_i));
    HCCoNN(i) = max(abs(tempHC));
    
end

% Plot things domainwise 
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/mrahaman/d_NBiC/Data/idxUpper.mat % Upper 1081 idx (IND1)
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'hot';  
% SZ
FF=zeros(47);
FF(idx) = SZCoNN;
FF = rot90(fliplr(FF));
pal1 = zeros(47);
pal1(idx) = SZCoNN;
FNCC = FF+pal1;
CLIM = [-max(abs(FNCC(:))) max(abs(FNCC(:)))];
T = 'SZ max connectivity';
[F,A,C,I] = plot_FNC(FNCC, CLIM, LABEL, RSN_I,'',T,MOD,cm,FNCC,0);
% HC

FFHC=zeros(47);
%HCCoNN = HCCoNN-SZCoNN;
FFHC(idx) = HCCoNN;
FFHC = rot90(fliplr(FFHC));
pal2 = zeros(47);
pal2(idx) = HCCoNN;
FNCCHC = FFHC+pal2;
%CLIM = [-max(abs(FNCCHC(:))) max(abs(FNCCHC(:)))];
T = 'HC max connectivity';
[F1,A1,C1,I1] = plot_FNC(FNCCHC, CLIM, LABEL, RSN_I,'',T,MOD,cm,FNCCHC,0);

%% top 10 high frequency shapes in both group 
 clear
 clc
 load('allPairshapesPD_HC.mat')
 HCP =p;
 clear p
 load allPairshapesPD_SZ.mat
 SZP = p;
 clear p
 SZs = load('allPairshapes_SZ.mat');
 SZshapes = SZs.allPairshapes;
 HCs = load('allPairshapes_HC.mat');
 HCshapes = HCs.allPairshapes;

 [top10SZ,top10SZIDX] = sort(SZP,'descend');
 [top10HC,top10HCIDX] = sort(HCP,'descend');
 figure()
 for i =1:10
    subplot(2,5,i)
    plot(SZshapes{top10SZIDX(i),1},'r','LineWidth',2);
    title(['Pair:' num2str(SZshapes{top10SZIDX(i),4}) ' PD:' num2str(top10SZ(i))])
    hold on
    plot(zeros(length(SZshapes{top10SZIDX(i),1}),1),'k','LineWidth',2)
    
    plot(HCshapes{top10HCIDX(i),1},'b','LineWidth',2);
    title(['Pair:' num2str(HCshapes{top10HCIDX(i),4})])
    hold on
    plot(zeros(length(HCshapes{top10HCIDX(i),1}),1),'b','LineWidth',2)
    set(gca,'fontsize', 14);
    box off
    axis tight
    axis off 
 end
 %% Random SZ /HC plotting
 clear
 clc
 load('allPairshapesPD_HC.mat')
 HCP =p_HC;
 load allPairshapesPD_SZ.mat
 SZP = p_SZ;
 load('allPairshapes_SZ.mat');
 SZshapes = SZ_allPairshapes;
 load('allPairshapes_HC.mat');
 HCshapes = HC_allPairshapes;
 [top10SZ,top10SZIDX] = sort(SZP,'descend');
 [top10HC,top10HCIDX] = sort(HCP,'descend');
 
 figure()
 for k =1:10
    subplot(2,5,k)
    %tshape = randi([1 8706],1);
    tshape = top10SZIDX(1000+k);
    plot(SZshapes{tshape,1},'r','LineWidth',2);
    title(['P:' num2str(SZshapes{tshape,4})])
    hold on
    plot(zeros((length(SZshapes{tshape,1})),1),'k','LineWidth',2)
    set(gca,'fontsize', 14);
    box off
    axis tight
    axis off 
 end
 
 figure()
 for k =1:10
    subplot(2,5,k)
    %pshape = randi([1 9285],1);
    pshape = top10HCIDX(k);
    plot(HCshapes{pshape,1},'b','LineWidth',2);
    title(['P' num2str(HCshapes{pshape,4})])
    hold on
    plot(zeros((length(HCshapes{pshape,1})),1),'k','LineWidth',2)
    set(gca,'fontsize', 14);
    box off
    axis tight
    axis off 
 end
 
%% 
 
 
 
 
 
 
 
 
 


