%% Statelts: Motifs summarization for SZ and HC group  
% This script first generate all the information for pairwise high
% frequent motifs from each group e.g., their PD's, shapes, groups etc.
% Then, compute EMD matrix for the groups 
% Compute tSNE and KDE using EMD matrix 
% Find 'k' peaks from each dynamic (HC/SZ) 
%% Getting the most dominant shapes from each pair of SZ subjects. 
clear
clc
%datadir = '/data/mialab/users/mrahaman/StateletFramework';                  % For the server
datadir = '/Users/mrahaman1/Documents/Statelet_V2/StateletFramework';        % Data directory (local) 
addpath(genpath(fullfile(datadir,'utilities')));
%% Get appropriate directories for SZ
PeaksDir = fullfile(datadir,'results','SZ_pairwise_peaks');
ShapesDir = fullfile(datadir,'results','SZ_shapes');
outdir = fullfile(datadir,'results','GroupwiseResults');
%% 1. Collect the motif's information from all pairs of SZ group for computing EMD across
count = 1;
for h =1:1081
  load(fullfile(PeaksDir,['pair_',num2str(h,'%04.f'),'.mat']));
  load (fullfile(ShapesDir,['pair_',num2str(h,'%04.f'),'.mat']));
  for pp = 1:length(peaks_i_j)
      Peak_shapes{count,1} = shapes(peaks_i_j(pp)).real_length;           % for grabing the real length 
      Peak_shapes{count,2} = shapes(peaks_i_j(pp)).extrapolated;          % extrapolated
      Peak_shapes{count,3} = shapes(peaks_i_j(pp)).whichsubject;          % subject
      Peak_shapes{count,4} = h;                                           % pair
      count = count+1;
  end
end

% 2. get EMD across all peak shapes

for l=1:size(Peak_shapes,1)
allshapes{l} = Peak_shapes{l,2};
end
% Compute EMD
simscore = zeros(size(Peak_shapes,1));
for j = 1:size(allshapes,2)
  x = allshapes{j};  
  for k = 1:size(allshapes,2)
    y = allshapes{k};
    x = x-min(x);
    x = x/sum(x);
    y = y-min(y);
    y = y/sum(y);
    simscore(j,k) = dEMD(x,y);
  end
end
save(fullfile(outdir,['group_wise_EMD_SZ','.mat']),'simscore');

% 3. tSNE using EMD distance
cd (fullfile(datadir,'main_files','step2_tSNE'))
exce_str = sprintf('python tSNE_groupwise_SZ.py');
system(exce_str)

% 4. KDE using EMD on all shapes from SZ group 
cd (fullfile(datadir,'main_files','step3_Probability Density(KDE)'))
groupwise_KDE(0)

% 5. Peak finding 
load(fullfile(outdir,'allshapesPD_SZ.mat')); 
    Ps = p; % p is the probability density 
    Ysub = load(fullfile(outdir,['tSNEcoord_SZ','.txt']));
    X = Ysub(:,1);
    Y = Ysub(:,2);
    pMat = zeros(1+ceil((max(X)-min(X))),1+ceil((max(Y)-min(Y))));
    for j = 1:length(X)
        x = 1+round(X(j)+(-1)*min(X));
        y = 1+round(Y(j)+(-1)*min(Y));
        pMat(x,y) = pMat(x,y)+Ps(j);
    end
    
%save(fullfile(outdir,'pMats_SZ.mat'),'pMats');
% Find the peaks
bpMat = imgaussfilt(pMat,4.5);
[cent, cm] = FastPeakFind(bpMat);
[xloc,yloc] = find(cm==1);
locs = [xloc yloc];
peaklocs = locs;
% Take 10 closest neigbours and get the max of them
coord = peaklocs;
%Ysub = load(fullfile(datadir,['pair_',num2str(pr,'%04.f'),'.txt']));
P = p; % PORBABILITY DENSITY
    
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
    
save(fullfile(outdir,['Dominat_Peaks_SZ_' '.mat']),'peaks_i_j')
%% ------------------------------- HC --------------------------------------------------------
clear
clc
%datadir = '/data/mialab/users/mrahaman/StateletFramework';                    % For the server
datadir = '/Users/mrahaman1/Documents/Statelet_V2/StateletFramework';          % For local 
addpath(genpath(fullfile(datadir,'utilities')));
%% Get directories for HC
PeaksDir = fullfile(datadir,'results','HC_pairwise_peaks');
ShapesDir = fullfile(datadir,'results','HC_shapes');
outdir = fullfile(datadir,'results','GroupwiseResults');
%% 1. Collect the motif's information from all pairs of HC group for computing EMD across
count = 1;
for h =1:1081
  load(fullfile(PeaksDir,['pair_',num2str(h,'%04.f'),'.mat']))
  load (fullfile(ShapesDir,['pair_',num2str(h,'%04.f'),'.mat']));
  for pp = 1:length(peaks_i_j)
      Peak_shapes{count,1} = shapes(peaks_i_j(pp)).real_length;           % for grabing the real length 
      Peak_shapes{count,2} = shapes(peaks_i_j(pp)).extrapolated;          % extrapolated
      Peak_shapes{count,3} = shapes(peaks_i_j(pp)).whichsubject;          % subject
      Peak_shapes{count,4} = h;                                           % pair
      count = count+1;
  end
end
%save(fullfile(outdir,['SZ_Peak_shapes' '.mat']),'SZ_peaks');
% get EMD across all peak shapes
simscore = zeros(size(Peak_shapes,1));
for l=1:size(Peak_shapes,1)
allshapes{l} = Peak_shapes{l,2};
end
% Compute EMD
for j = 1:size(allshapes,2)
  x = allshapes{j};  
  for k = 1:size(allshapes,2)
    y = allshapes{k};
    x = x-min(x);
    x = x/sum(x);
    y = y-min(y);
    y = y/sum(y);
    simscore(j,k) = dEMD(x,y);
  end
end
save(fullfile(outdir,['group_wise_EMD_HC','.mat']),'simscore');
% 3. tSNE
cd (fullfile(datadir,'main_files','step2_tSNE'))
exce_str = sprintf('python tSNE_groupwise_HC.py');
system(exce_str)

% 4. KDE using EMD on all shapes from SZ group 
cd (fullfile(datadir,'main_files','step3_Probability Density(KDE)'))
groupwise_KDE(1)

% 5. Peak finding 
load(fullfile(outdir,'allshapesPD_HC.mat')); 
    Ps = p; % p is the probability density 
    Ysub = load(fullfile(outdir,['tSNEcoord_HC','.txt']));
    X = Ysub(:,1);
    Y = Ysub(:,2);
    pMat = zeros(1+ceil((max(X)-min(X))),1+ceil((max(Y)-min(Y))));
    for j = 1:length(X)
        x = 1+round(X(j)+(-1)*min(X));
        y = 1+round(Y(j)+(-1)*min(Y));
        pMat(x,y) = pMat(x,y)+Ps(j);
    end
    
%save(fullfile(outdir,'pMats_SZ.mat'),'pMats');
% Find the peaks
bpMat = imgaussfilt(pMat,4.5);
[cent, cm] = FastPeakFind(bpMat);
[xloc,yloc] = find(cm==1);
locs = [xloc yloc];
peaklocs = locs;
% Take 10 closest neigbours and get the max of them
coord = peaklocs;
%Ysub = load(fullfile(datadir,['pair_',num2str(pr,'%04.f'),'.txt']));
P = p; % PORBABILITY DENSITY
    
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
    
save(fullfile(outdir,['Dominat_Peaks_HC_' '.mat']),'peaks_i_j')
%%


