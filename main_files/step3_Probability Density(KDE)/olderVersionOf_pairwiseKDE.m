% Kernel Density Estimator for given matrix of distances
clear
clc
outdir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis';
if(cl==1)
    gr = 'HC';
    dataDir ='/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/HC_PW_simscores'; 
else
    gr = 'SZ'; 
    dataDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/SZ_PW_simscores';          % distances
    %sDir    = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/SZ_pairwise_shapelet';        % shapes
end
dirs    = dir(dataDir);
D       = 1;     % dimensions
%% get the generelized bandwidths
for i = 1:size(dirs,1)-2
    i
    load(fullfile(dataDir,['pair_' num2str(i,'%04.f') '.mat']));
    %load(fullfile(sDir,['pair_' num2str(i,'%03.f') '.mat']));
    N = size(simscore,1);
    p = zeros(1,N);  
    for j = 1:N
        dists = squeeze(simscore(j,:));
        std   = sqrt(mean(dists.^2)); 
        h     = 1.06*std*N^-(1/5);               % bandwidth of the kernel, optimal: 1.06*STD*N^(-1/5) for gaussian
        hs(i,j) = h;
    end
end
% Averaged the bandwidths over each pair 
stds = zeros(314,1);
for jj = 1:size(hs,1)
    kl = hs(jj,:);
    hh = find(kl~=0);
    stds(jj) = mean(kl(hh));
end
save(fullfile(outdir,['bandwidths_' gr '.mat']),'stds')
%% For groups. 
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
%load (fullfile(outdir,['bandwidths_' gr '.mat']))

for i = 1:size(dirs,1)-2
    i
    load(fullfile(dataDir,['pair_' num2str(i,'%04.f') '.mat']));
    %load(fullfile(sDir,['pair_' num2str(i,'%03.f') '.mat']));
    N = size(simscore,1);
    p = zeros(1,N);
    h = stds(i);
    
    for j = 1:N
        dists = squeeze(simscore(j,:));
        Gau = zeros(1,N-1);
        for k = 1:N-1
            Gau(k)= ((2*pi*h^2)^-(D/2))*exp(-(dists(1,k)^2)/(2*h^2));  % Gaussian kerenel on point(shapes) 'j' of subejct i  
        end
        p(j)  = mean(Gau);
    end
    pd{i} = p;                                                         % Storing probability density for each of the subject 
end
save(fullfile(outdir,['probabilityDensity_' gr '.mat']),'pd')
%% For all the shapes of 1081 pairs HC
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')  
load allPairEMD_HC.mat
N = size(simscore,1);
p = zeros(1,N);  
 for j = 1:N
        j
        dists = squeeze(simscore(j,:));
        std   = sqrt(mean(dists.^2)); 
        h     = 1.06*std*N^-(1/5);                                   % bandwidth of the kernel, optimal: 1.06*STD*N^(-1/5) for gaussian
        hs(j) = h;
 end
 
 h = 1.8430;
 p = zeros(1,N); 
 D       = 1;     % dimensions
for j = 1:N
    j
    dists = squeeze(simscore(j,:));
    Gau = zeros(1,N-1);
    for k = 1:N-1
            Gau(k)= ((2*pi*h^2)^-(D/2))*exp(-(dists(1,k)^2)/(2*h^2));  % Gaussian kerenel on point(shapes) 'j' of subejct i  
    end
    p(j)  = mean(Gau);
end 
%% SZ
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')  
load allPairEMD_SZ.mat
N = size(simscore,1);
p = zeros(1,N);  
 for j = 1:N
        j
        dists = squeeze(simscore(j,:));
        std   = sqrt(mean(dists.^2)); 
        h     = 1.06*std*N^-(1/5);                           % bandwidth of the kernel, optimal: 1.06*STD*N^(-1/5) for gaussian
        hs(j) = h;
 end
 
 h = mean(hs(:));
 p = zeros(1,N); 
 D       = 1;     % dimensions
for j = 1:N
    j
    dists = squeeze(simscore(j,:));
    Gau = zeros(1,N-1);
    for k = 1:N-1
            Gau(k)= ((2*pi*h^2)^-(D/2))*exp(-(dists(1,k)^2)/(2*h^2));  % Gaussian kerenel on point(shapes) 'j' of subejct i  
    end
    p(j)  = mean(Gau);
end 
%%
    