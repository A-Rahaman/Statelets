% rFdyn is M (subj) x Nwin (windows) x nCC (ICNpairs) dFNC matrix

[M, Nwin, nCC] = size(rFdyn);
%%

SP = cell(M,1);
k1_peaks = zeros(1,M);
for ii = 1:M
    fprintf('Working on subject %d of %d\n', ii, M)
    submat = squeeze(rFdyn(ii,:,:));    
    DEV = std(submat, [], 2);
    %DEV =std(submat-repmat(mean(submat),Nwin,1), [], 2);
    [xmax,imax,xmin,imin] = extrema(DEV); % find the extrema
    pIND = sort(imax);
    k1_peaks(ii) = length(pIND);
    SP{ii} = submat(pIND,:);
end
%%
% % % Determine approximate k2 based on the subject clusters 
ktest = 2:20;dmethod = 'corr';
ntry = 10;
bestx = zeros(1,ntry);
%nrep = 150;
for pp = 1:ntry
    sublist = randperm(M);
    SPsamp = cell2mat(SP(sublist));
    SPsamp = SPsamp(1:10:end,:);
    SSE = zeros(1,length(ktest));
    R = zeros(1,length(ktest));
    nrep=5;
    for k2 = ktest
        fprintf('---------------- k = %d ------------------\n', k2)
        [IDX, C, SUMD, D] = kmeans(SPsamp, k2, 'distance', dmethod, 'Replicates', nrep, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop');
        [SSE(k2-1), R(k2-1)] = cluster_goodness(D, IDX);
    end
    
    [bestx(pp), FFF] = fit_L_tocurve_area(ktest,R, 1);
    drawnow;
end
%% Obtain centroids for all subject exemplar window data 
SPflat = cell2mat(SP);
for k2 = 2:8 
    [IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop');
 
    save(fullfile(['FNC_group_clusters_k_' num2str(k2) ]), 'SPflat', 'IDXp','Cp');

    G = figure('Color', 'w', 'Name', ['FNC_group_clusters_k' num2str(k2) 'subj_exemplar'], 'Position', CPOS);
    for ii = 1:k2, 
        subplot(1,k2,ii);
        stmat = vec2mat(median(SPflat(IDXp == ii ,:)),1);
        imagesc(stmat, [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp == ii),  round(100*sum(IDXp == ii)/length(IDXp))));
    end
    
    clear IDXp  Cp SUMDp Dp
end
%%
