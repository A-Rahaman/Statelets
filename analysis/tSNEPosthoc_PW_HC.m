% Once we finish all the steps following, we will get EMD for all the dominant shapes across the pairs 
% Steps include:
%       1. coordinate Mapping
%       2. Blurr and 2D peak finding
%       3. Getting 10 dominants to shrink down the set of peaks. Keep 1 per 10
%       4. Final list of peaks per pair
%       5. Get all the peak shapes from all the pairs
%       6. EMD across the shapes 
% Need to do afterwards:
%    - running tSNE using that EMD matrix
%    - rerun the above on one tSNE map
%    - get a final set of peak shapes for a given group (HC)
%%
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

%save(fullfile(outdir,'pMats_HC.mat'),'pMats');
%% Find the peaks
%[pks,locs] = findpeaks(bpMat,'MinPeakDistance',6);
%[pks,locs] = findpeaks(bpMat);
%[cent, cm] = FastPeakFind(pMat);

peaklocs = cell(size(pMats,1),1);
bpMat = imgaussfilt(pMat,4.5);
%bpMat = bpMat/max(bpMat(:));
[cent, cm] = FastPeakFind(bpMat);
%ii = imregionalmax(bpMat);
%[xloc,yloc,~] = find(ii);        % coordinates
[xloc,yloc] = find(cm==1);
locs = [xloc yloc];
peaklocs{i} = locs;

%save(fullfile(outdir,'Peaklocs_HC.mat'),'peaklocs');
% imagesc(bpMat)
% hold on
% for fn=1:length(x)
%     text(y(fn),x(fn),'*','FontSize',25)
% end
% hold off
%mesh(0.01*cm+bpMat)
%mesh(ii+bpMat)
%% Take 10 closest neigbours and get the max of them
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

%% Getting the most dominant shapes across the pairs 
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

%% Leaders
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



%%





