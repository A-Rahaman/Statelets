function fr = groupwise_KDE(gr)
% Input: EMD matrix [NXN] where N is the number of shapes in the collection  
% Returns the probability density for each motif 
% Use a different bandwidth for each given pair
% It computes the mean stds across all the motifs of the pair and use this mean std to compute the h (bandwidths)
% Then using gaussian equation it label each motif with probability density
% Returns: PD for each shape (considers it as a point in the porbability space)
% gr = group, 1 = HC, 0 = SZ
% Assign appropriate directory for groupwise results computed in earleir
% steps. We need EMD matrix for all the motifs in the group.   
maindir = '/Users/mrahaman1/Documents/Statelet_V2/StateletFramework/results/GroupwiseResults';
%maindir = '/data/mialab/users/mrahaman/StateletFramework/results/GroupwiseResults';
if(gr == 1)
    datadir = fullfile(maindir,'group_wise_EMD_HC.mat');
    g = 'HC';
else
    datadir = fullfile(maindir,'group_wise_EMD_SZ.mat');
    g = 'SZ';
end
D       = 1;     % dimensions
%% shape wise porbability density 
load(datadir); % EMD matrix for all the motifs in the group.   
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
save(fullfile(maindir,['allshapesPD_' g '.mat']),'p')      % Storing pd for each of the pair across all subjects within SZ/HC
%% --------- 