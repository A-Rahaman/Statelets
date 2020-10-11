function fr = pairwise_peak_shapes(pr,gr)
% The script generates dominant motifs in a given pair 'pr' across all
% the subejects of group 'gr'.
% Usually, we run the script on the server using 'arrayjob' properties of 'SLURM' workload manager  
% Steps include:
%       1. coordinate Mapping
%       2. Blurr and 2D peak finding
%       3. Getting R = 10 dominants to shrink down the set of peaksa nd Keeping 1 per group
%       4. Final list of peaks per pair
%       5. Get all the peak shapes from all the pairs
% Returns a list of peak motifs 
% Need to do afterwards:
%    - running tSNE using that EMD matrix
%    - rerun the above on one tSNE map
%    - get a final set of peak shapes for a given group (SZ)
%maindir = '/Users/mrahaman1/Documents/Statelet_V2/StateletFramework/results';  % local 
maindir = '/data/mialab/users/mrahaman/StateletFramework/results';              % Server
if(gr == 1)
    datadir = fullfile(maindir,'HC_tSNE');
    outdir = fullfile(maindir,'HC_pairwise_peaks');
    datadir1 = fullfile(maindir,'HC_pairwise_PD');
else
    datadir = fullfile(maindir,'SZ_tSNE');
    outdir = fullfile(maindir,'SZ_pairwise_peaks');
    datadir1 = fullfile(maindir,'SZ_pairwise_PD');
end

load(fullfile(datadir1,['pair_',num2str(pr,'%04.f'),'.mat']));
%% Map the coordinates a 2D matrix and weight each of the cell with probability density
pMats = cell(314,1);
i = pr;
 
    Ps = p; % p is the probability density
    Ysub = load(fullfile(datadir,['pair_',num2str(pr,'%04.f'),'.txt']));
    X = Ysub(:,1);
    Y = Ysub(:,2);
    pMat = zeros(1+ceil((max(X)-min(X))),1+ceil((max(Y)-min(Y))));
    for j = 1:length(X)
        x = 1+round(X(j)+(-1)*min(X));
        y = 1+round(Y(j)+(-1)*min(Y));
        pMat(x,y) = pMat(x,y)+Ps(j);
    end
    
%save(fullfile(outdir,'pMats_SZ.mat'),'pMats');
%% Find the peaks
peaklocs = cell(size(pMats,1),1);
bpMat = imgaussfilt(pMat,4.5);
%bpMat = bpMat/max(bpMat(:));
[cent, cm] = FastPeakFind(bpMat);
%ii = imregionalmax(bpMat);
%[xloc,yloc,~] = find(ii);        % coordinates
[xloc,yloc] = find(cm==1);
locs = [xloc yloc];
peaklocs{i} = locs;
%% Take 10 closest neigbours and get the max of them
coord = peaklocs{i};
Ysub = load(fullfile(datadir,['pair_',num2str(pr,'%04.f'),'.txt']));
P = p;
    
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
    % peaksIJ{i} = peaks_i_j; 
    % Filter some of the peaks which are very close to each other. keep maxp always
%     peaks_i = peaksIJ{i};
%     for l = 1:length(peaks_i)
%         
%     end
    % --------------------------------------------------------------
save(fullfile(outdir,['pair_' num2str(pr,'%04.f') '.mat']),'peaks_i_j')
end




