function fr = getmotifsPairWise(pr)
% The script performs local motifs discovery from the time series of connection 'pr' of both group SZ and HC 
% webutils.htmlrenderer('basic'); 
% Group the subjects into SZ and HC 
% For each pair, 
%     - take all the time series (across the subjects) within a group, 
%     - Extract local motifs from each time series
%     - Store all required information for local motifs 
%     - run EMD -> tSNE -> KDE -> Peakfinder to 
%     - find dominants/frequent motifs across the group given that pair  
% Return: extracted motifs across all SZ/HC for a given pair pr, other
% informations, EMD scores of all the shapes 
% Usually, we run the script on the server using 'arrayjob' properties of 'SLURM' workload manager  
datadir = '/data/mialab/users/mrahaman/StateletFramework';                    % For the server
%datadir = '/Users/mrahaman1/Documents/Statelet_V2/StateletFramework';        % For local 
addpath(genpath(fullfile(datadir,'utilities')));
addpath(genpath(fullfile(datadir,'main_files','step2_tSNE')));
addpath(genpath(fullfile(datadir,'main_files','step3_Probability Density(KDE)')));
addpath(genpath(fullfile(datadir,'main_files','step4_finding dominants')));

% Loading the time series data from the data directory (dFNC time course) 
load (fullfile(datadir,'data','dFNCs.mat'));

load (fullfile(datadir,'data','candidates_30to50.mat')); % Candidate for length 30 to 50
load (fullfile(datadir,'data','SZsubjectIdx.mat'));      % Subeject's label SZ
load (fullfile(datadir,'data','HCsubjectIdx.mat'));      % Subejct's label HC
TotalWindows = 136;

%% ----------------------------------------- SZ group ----------------------------------

data_SZ  = squeeze(rFdyn(SZIdx,:,pr));       % all subjects' data for given pair 'pr'
c = 1;                                       % shapes counter
shapes = initshape();                        % Initializing the 'shapes' structure
for sz = 1:size(data_SZ,1) 
tdata_SZ = squeeze(data_SZ(sz,:));           % subject's data
sizeWISEminSigma = [];
minSigmaidx = [];
allSigmoids = {};

% ----------------------------- check all candidates of different lengths ---------------------------
for j = 1:length(candidates)                
    cand_j={};
    cand_j = candidates{:,j};
    slices = Slicer(length(cand_j{1}),TotalWindows);
    meanSigma = zeros(1,length(cand_j));
 for k = 1:length(cand_j)
    sigmoid   = zeros(length(cand_j)-1,1);
    cand_jj = [];
    c1 = 1;
    for m = 1:length(cand_j)
        if(m~=k)
            cand_jj(c1) = m;
           %--------- Normalization-----------------
            shape1 = tdata_SZ(cand_j{k,1});
            shape1 = shape1-min(shape1);
            shape1 = shape1/max(shape1);
            shape2 = tdata_SZ(cand_j{m,1});
            shape2 = shape2-min(shape2);
            shape2 = shape2/max(shape2);
            %------------------------------
            sigmoid(c1) = dEMD(shape1,shape2);
            c1 = c1+1;
        end
    end
lMinimas     = findLocalMinimas(slices,cand_j,cand_jj,sigmoid);
meanSigma(k) = mean(lMinimas);
allSigmoids{j,k} = sigmoid;
end
[sizeWISEminSigma(j,1),minSigmaidx(j,1)] = min(meanSigma);
end
%------------------------------- end of checking candidates of different length -----------------------------
% Get the target size and corresponding candidate.
%[~,nearbySigma]           = getNearbySigma(sizeWISEminSigma,5)
[~,nearbySigma]           = min(sizeWISEminSigma);
ind_BW                    = minSigmaidx(nearbySigma);
allWavelets               = candidates{nearbySigma};                            % All wavelets of that size
Sigmoids1                 = allSigmoids{nearbySigma,ind_BW};
Sigmoids = zeros(1,size(allWavelets,1));

if(ind_BW==1)
    Sigmoids(ind_BW)=0;
    Sigmoids((ind_BW+1):end)=Sigmoids1;
elseif(ind_BW==size(allWavelets,1))
    Sigmoids(ind_BW)=0;
    Sigmoids(1:(ind_BW-1)) = Sigmoids1;
else
tempSig1 = Sigmoids1(1:ind_BW-1);
tempSig2 = Sigmoids1(ind_BW:end);
Sigmoids(1:length(tempSig1)) = tempSig1;
Sigmoids(ind_BW) = 0;
Sigmoids((ind_BW+1):end) = tempSig2;
end
%--------------------------------

[sortedSigmoids,SSidx] = sort(Sigmoids);
Index_of_hSigmoids     = SSidx(sortedSigmoids<=mean(sortedSigmoids));
maxidx = [];
maxidx(1) = Index_of_hSigmoids(1);
for i = 2:length(Index_of_hSigmoids)
    flag = 0;
    for j =1:length(maxidx)
        if(~isempty(intersect(allWavelets{Index_of_hSigmoids(i),1},allWavelets{maxidx(j),1})))
            flag = 1;
            break;
        end
    end
    if(flag~=1)
        maxidx(end+1) = Index_of_hSigmoids(i);
    end
end
cands_p = candidates{1,nearbySigma};  
maxidx = unique(maxidx);                                          
occur = maxidx;                                                       % occurences
                         
for j =1:length(occur)
    x = tdata_SZ(cands_p{occur(j)});                    
    shapes(c).real_length  = x;                                       
    if(length(x)<50)
        x = interp1( x, linspace( 1, length(x), 50 ) );                 
        x = x';
    end
  allshapes_pr{c,1} = x;                                              
  shapes(c).extrapolated     = x;                                    
  shapes(c).occurences       = occur;                                 % occurences 
  shapes(c).whichsubject     = SZIdx(sz);                             % subject
  shapes(c).sizeWISEminsigma = sizeWISEminSigma;                      
  shapes(c).locExactshapelet = ind_BW;                                
  shapes(c).isLeader         = 0;                                     
  if(j==1)
      shapes(c).isLeader = 1;                                         
  end
  c = c+1;
end
end
save(fullfile(datadir,'results','SZ_shapes',['pair_' num2str(pr,'%04.f') '.mat']),'shapes');
%% get the EMD score matrix across all the shapes of a pair
simscore = zeros(size(allshapes_pr,1));
for j = 1:size(allshapes_pr,1)
  x = allshapes_pr{j,1};  
  for k = 1:size(allshapes_pr,1)
    if j~=k
    y = allshapes_pr{k,1};
    x = x-min(x);
    x = x/sum(x);
    y = y-min(y);
    y = y/sum(y);
    simscore(j,k) = dEMD(x,y); 
    end 
  end
end
save(fullfile(datadir,'results','SZ_EMD_matrices',['pair_' num2str(pr,'%04.f') '.mat']),'simscore');
%% Running the tSNE
cd(fullfile(datadir,'main_files','step2_tSNE'))
exce_str = sprintf('python tSNE_for_pairwise_dominant_motifs.py %d %d',pr,0);
system(exce_str)
%% Running KDE pairwise 

pairwise_KDE(pr,0) 

%% Running pairwise peak finder.SZ = 0 

pairwise_peak_shapes(pr,0)

%%   ----------------------------------- HC group -----------------------------
%{
datadir = '/data/mialab/users/mrahaman/StateletFramework';                  % For the server
%datadir = '/Users/mrahaman1/Documents/Statelet_V2/StateletFramework';        % For local 
addpath(genpath(fullfile(datadir,'utilities')));
addpath(genpath(fullfile(datadir,'main_files','step2_tSNE')));
addpath(genpath(fullfile(datadir,'main_files','step3_Probability Density(KDE)')));
addpath(genpath(fullfile(datadir,'main_files','step4_finding dominants')));
load (fullfile(datadir,'data','dFNCs.mat'));
load (fullfile(datadir,'data','candidates_30to50.mat'));
load (fullfile(datadir,'data','SZsubjectIdx.mat'));
load (fullfile(datadir,'data','HCsubjectIdx.mat'));
TotalWindows = 136;
%%
%}
data_HC  = squeeze(rFdyn(CtIdx,:,pr));              % all subjects' data for given pair 'pr'
c = 1;                                            	% shapes counter
shapes = initshape();                               % Initializing the 'shapes' structure
for hc = 1:size(data_HC,1)                      
tdata_HC = squeeze(data_HC(hc,:));                 % subject's data
sizeWISEminSigma = [];
minSigmaidx = [];
allSigmoids = {};

% ----------------------------- check all candidaets of different lengths ---------------------------
for j = 1:length(candidates)                
    cand_j={};
    cand_j = candidates{:,j};
    slices = Slicer(length(cand_j{1}),TotalWindows);
    meanSigma = zeros(1,length(cand_j));
 for k = 1:length(cand_j)
    sigmoid_HC   = zeros(length(cand_j)-1,1);
    cand_jj = [];
    c1 = 1;
    for m = 1:length(cand_j)
        if(m~=k)
            cand_jj(c1) = m;
%           --------- Normalization-----------------
            shape1 = tdata_HC(cand_j{k,1});
            shape1 = shape1-min(shape1);
            shape1 = shape1/max(shape1);
            shape2 = tdata_HC(cand_j{m,1});
            shape2 = shape2-min(shape2);
            shape2 = shape2/max(shape2);
%            ------------------------------
            sigmoid_HC(c1) = dEMD(shape1,shape2);
            c1 = c1+1;
        end
    end
lMinimas     = findLocalMinimas(slices,cand_j,cand_jj,sigmoid_HC);
meanSigma(k) = mean(lMinimas);
allSigmoids{j,k} = sigmoid_HC;
end
[sizeWISEminSigma(j,1),minSigmaidx(j,1)] = min(meanSigma);
end
% ----
%% Get the target size and corresponding candidate.
[~,nearbySigma]           = min(sizeWISEminSigma);
ind_BW                    = minSigmaidx(nearbySigma);
allWavelets               = candidates{nearbySigma};                            % All wavelets of that size
Sigmoids1                 = allSigmoids{nearbySigma,ind_BW};
Sigmoids_HC = zeros(1,size(allWavelets,1));

if(ind_BW==1)
    Sigmoids_HC(ind_BW)=0;
    Sigmoids_HC((ind_BW+1):end)=Sigmoids1;
elseif(ind_BW==size(allWavelets,1))
    Sigmoids_HC(ind_BW)=0;
    Sigmoids_HC(1:(ind_BW-1)) = Sigmoids1;
else
tempSig1 = Sigmoids1(1:ind_BW-1);
tempSig2 = Sigmoids1(ind_BW:end);
Sigmoids_HC(1:length(tempSig1)) = tempSig1;
Sigmoids_HC(ind_BW) = 0;
Sigmoids_HC((ind_BW+1):end) = tempSig2;
end
%--------------------------------

[sortedSigmoids,SSidx] = sort(Sigmoids_HC);
Index_of_hSigmoids     = SSidx(sortedSigmoids<=mean(sortedSigmoids));
maxidx = [];
maxidx(1) = Index_of_hSigmoids(1);
for i = 2:length(Index_of_hSigmoids)
    flag = 0;
    for j =1:length(maxidx)
        if(~isempty(intersect(allWavelets{Index_of_hSigmoids(i),1},allWavelets{maxidx(j),1})))
            flag = 1;
            break;
        end
    end
    if(flag~=1)
        maxidx(end+1) = Index_of_hSigmoids(i);
    end
end
cands_p = candidates{1,nearbySigma};  
maxidx = unique(maxidx);                                           
occur = maxidx;                                                      
                         
for j =1:length(occur)
  x = tdata_HC(cands_p{occur(j)});                    
  shapes(c).real_length  = x;                                         
  if(length(x)<50)
      x = interp1( x, linspace( 1, length(x), 50 ) );                 
      x = x';
  end
  allshapes_pr_HC{c,1} = x; 
  shapes(c).extrapolated     = x;                                     
  shapes(c).occurences       = occur;                                 % occurences 
  shapes(c).whichsubject     = CtIdx(hc);                             % Which subject
  shapes(c).sizeWISEminsigma = sizeWISEminSigma;                      % min Sigma for that size 
  shapes(c).locExactshapelet = ind_BW;                                % which exact candidate of that length performed best
  shapes(c).isLeader         = 0;                                     % default not leader 
  if(j==1)
      shapes(c).isLeader = 1;                                         
  end
  c = c+1;
end
end
save(fullfile(datadir,'results','HC_shapes',['pair_' num2str(pr,'%04.f') '.mat']),'shapes');
%% get the EMD score matrix across all the shapes of a pair
simscore = zeros(size(allshapes_pr_HC,1));
for j = 1:size(allshapes_pr_HC,1)
  x = allshapes_pr_HC{j,1};  
  for k = 1:size(allshapes_pr_HC,1)
    if j~=k
    y = allshapes_pr_HC{k,1};
    x = x-min(x);
    x = x/sum(x);
    y = y-min(y);
    y = y/sum(y);
    simscore(j,k) = dEMD(x,y); 
    end 
  end
end
save(fullfile(datadir,'results','HC_EMD_matrices',['pair_' num2str(pr,'%04.f') '.mat']),'simscore');

%% Running tSNE
cd(fullfile(datadir,'main_files','step2_tSNE'))
exce_str1 = sprintf('python tSNE_for_pairwise_dominant_motifs.py %d %d',pr,1);
system(exce_str1)

%% Running KDE pairwise 

pairwise_KDE(pr,1) 

%% Running pairwise peak finder 

pairwise_peak_shapes(pr,1)


%% test plotting
%{
colors = 'r';
figure(3); cla;
plot(1:length(tdata),tdata)
hold on
for p = 1: length(maxidx)
    templine= NaN(1,136);
    templine(allWavelets{maxidx(p)})= tdata(allWavelets{maxidx(p)});
    plot(1:length(tdata),templine,'Color',colors(1));
    hold on
end
%}
%% 


