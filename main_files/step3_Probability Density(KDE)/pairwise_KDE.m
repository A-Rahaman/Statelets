function fr = pairwise_KDE(pr,gr)
% The script computes the probability density of all the motifs extratced from pair 'pr' of group 'gr'  
% Input: EMD matrix [NXN] where N is the number of shapes int he collection  
% Use a different bandwidth for each given pair
% It computes the mean stds across all the shapes of the pair and use this mean std to compute the h (bandwidths)
% Then using gaussian equation it label each shape with probability density
% Returns: PD for each shape (consider it as a point in the porbability space)
% gr = group, e.g., gr = 1 for HC and 0 for  SZ
% pr = pair  e.g.,  pr = 1 to 1081

% Usually, we run the script on the server using 'arrayjob' properties of 'SLURM' workload manager  

%maindir = '/Users/mrahaman1/Documents/Statelet_V2/StateletFramework/results';
maindir = '/data/mialab/users/mrahaman/StateletFramework/results';
if(gr == 1)
    datadir = fullfile(maindir,'HC_EMD_matrices');
    outdir = fullfile(maindir,'HC_pairwise_PD');
else
    datadir = fullfile(maindir,'SZ_EMD_matrices');
    outdir = fullfile(maindir,'SZ_pairwise_PD');
end
D       = 1;     % dimensions
%% shape wise porbability density 
load(fullfile(datadir,['pair_' num2str(pr,'%04.f') '.mat'])); % load the EMD matrix for pair 'pr' of group 'gr'
N = size(simscore,1);
std   = sqrt(mean(simscore(:).^2));    % Standard deviation across the distances  
h     = 1.06*std*N^-(1/5);             % bandwidth
for j = 1:N
    dists = squeeze(simscore(j,:));
    Gau = zeros(1,N-1);
    for k = 1:N-1
        Gau(k)= ((2*pi*h^2)^-(D/2))*exp(-(dists(1,k)^2)/(2*h^2));  % Gaussian kerenel on point(shapes) 'j' of pair pr 
    end
    p(j)  = mean(Gau);
end                                                         
save(fullfile(outdir,['pair_' num2str(pr,'%04.f') '.mat']),'p')      % Storing pd for each of the pair across all subjects within SZ/HC
%% --------- 