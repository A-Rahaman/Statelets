clear all; clc
%%addpath(genpath('/export/mialab/users/eallen/export_fig'))
%inname = '/export/mialab/hcp/dynamics/FNC_dynamicSWTukey18_zscored.mat';
%FIGoutputdir = '/export/mialab/hcp/dynamics/figures/No_L1';
addpath(genpath('/export/mialab/users/eswar/fbirn_p3/scripts/cluster'))
addpath(genpath('/export/mialab/users/eswar/fbirn_p3/scripts/network_estimation'))
addpath(genpath('/export/mialab/users/eswar/software/nic'))
%inname = '/export/mialab/hcp/dynamics/FNC_dynamicSWTukey18_L1ICOV_zscored';
%FIGoutputdir = '/export/mialab/hcp/dynamics/figures/L1';

FIGoutputdir = '/export/mialab/users/eswar/fbirn_p3/results_C100_vn/figures/dyn_L1';
inname = '/export/mialab/users/eswar/fbirn_p3/results_C100_vn/FNCdynamics/FNC_dynamicSWTukey22_L1ICOV_zscored_scr';

TIMEPOS = [ 680   752   497   346];
CPOS = [ 180         886        1717         212];
load(inname)
% loads FNCdynflat [Nwin, M, C*(C-1)/2]
getICA_SubjGroupID;
FNCdynflat = FNCdynflat(:,keepID,:);
[Nwin, M, nCC] = size(FNCdynflat);
%% keep only the RSNs
[L, RSN_I, ART_I, fMOD] = comp_labels_fb_C100_vn;
Fdyn = zeros(M, Nwin, length(RSN_I)*(length(RSN_I)-1)/2);
for ii = 1:M;
    temp = squeeze(FNCdynflat(:,ii,:));
    temp = vec2mat(temp);
    temp = temp(:,RSN_I,RSN_I);
    temp = mat2vec(temp);
    Fdyn(ii,:,:) = temp;
end
clear FNCdynflat

% Fdyn [M, Nwin, C*(C-1)/2]
%% remove the mean across time from the matrix
% dFdyn = Fdyn - repmat(mean(Fdyn,2), [1, Nwin, 1]);
% Fdynmean = 

%% Find peaks in the data for each subject
SP = cell(M,1);
k1_peaks = zeros(1,M);
for ii = 1:M
    fprintf('Working on subject %d of %d\n', ii, M)
    submat = squeeze(Fdyn(ii,:,:));    
    DEV = std(submat, [], 2);
    %DEV =std(submat-repmat(mean(submat),Nwin,1), [], 2);
    [xmax,imax,xmin,imin] = extrema(DEV); % find the extrema
    pIND = sort(imax);
    k1_peaks(ii) = length(pIND);
    SP{ii} = submat(pIND,:);
end
%%
% % % Determine approximate k2 based on the subject clusters (roughly 2000 samples)
ktest = 2:20;
for pp = 1:3
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
    
    [bestx, F] = fit_L_tocurve_area(ktest,R, 1);
    drawnow;
end
%% Cluster
dmethod = 'city';
nrep = 150;
% --------------------------------------------------------------------------
% All the data
SPflat = cell2mat(SP);
for k2 = 2:9
[IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop');
save(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) 'scrubbed']), 'SPflat', 'IDXp','Cp');
%load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) 'scrubbed']), 'SPflat', 'IDXp','Cp');

G = figure('Color', 'w', 'Name', ['FNC_group_clusters_k' num2str(k2) 'scrubbed'], 'Position', CPOS);
for ii = 1:k2, 
    subplot(1,k2,ii);
    imagesc(vec2mat(median(SPflat(IDXp == ii ,:)),1), [-.5,.5]); axis square;
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%d (%d%%)', sum(IDXp == ii),  round(100*sum(IDXp == ii)/length(IDXp))));
end
IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
end
% --------------------------------------------------------------------------

%%
% Half-split analysis
k2 = 7;
sublist = randperm(M);
SPflat1 = cell2mat(SP(sublist(1:floor(M/2))));
SPflat2 = cell2mat(SP(sublist(floor(M/2)+1:end)));
load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) 'scrubbed']), 'SPflat', 'IDXp','Cp');

[IDXp1, Cp1, SUMDp1, Dp1] = kmeans(SPflat1, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop', 'Start', Cp);

save(fullfile(fileparts(inname), ['FNC_half1_clusters_k2_' num2str(k2) 'scrubbed']), 'SPflat1', 'IDXp1','Cp1', 'sublist');

[IDXp2, Cp2, SUMDp2, Dp2] = kmeans(SPflat2, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop', 'Start', Cp);

save(fullfile(fileparts(inname), ['FNC_half2_clusters_k2_' num2str(k2) 'scrubbed']), 'SPflat2', 'IDXp2','Cp2', 'sublist');

G = figure('Color', 'w', 'Name', ['FNC_halfgroup_clusters_k' num2str(k2) 'scrubbed'], 'Position', CPOS);
for ii = 1:k2,
    subplot(2,k2,ii);
    imagesc(vec2mat(median(SPflat1(IDXp1 == ii ,:)),1), [-.5,.5]); axis square;
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%d (%d%%)', sum(IDXp1 == ii),  round(100*sum(IDXp1 == ii)/length(IDXp1))));

    subplot(2,k2,ii+k2);
    imagesc(vec2mat(median(SPflat2(IDXp2 == ii ,:)),1), [-.5,.5]); axis square;
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%d (%d%%)', sum(IDXp2 == ii),  round(100*sum(IDXp2 == ii)/length(IDXp2))));
end
IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);

%--------------------------------------------------------------------------
%%  HC SZ
%sublist = randperm(M);
HCid = find(ica_ursi_group_keepID == 0);
SZid = find(ica_ursi_group_keepID == 1);

SPflatHC = cell2mat(SP(ica_ursi_group_keepID==0));
SPflatSZ = cell2mat(SP(ica_ursi_group_keepID==1));
%%
for k2 = 2:9
    [IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, k2, 'distance', dmethod, 'Replicates', 10, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop');
    
    [IDXp1, Cp1, SUMDp1, Dp1] = kmeans(SPflatHC, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop', 'Start', Cp);

    save(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatHC', 'IDXp1','Cp1', 'HCid');
    
    [IDXp2, Cp2, SUMDp2, Dp2] = kmeans(SPflatSZ, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop', 'Start', Cp);

    save(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid');

    G = figure('Color', 'w', 'Name', ['FNC_DIAGgroup_clusters_k' num2str(k2) '_scrubbed'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        imagesc(vec2mat(median(SPflatHC(IDXp1 == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp1 == ii),  round(100*sum(IDXp1 == ii)/length(IDXp1))));
        
        subplot(2,k2,ii+k2);
        imagesc(vec2mat(median(SPflatSZ(IDXp2 == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp2 == ii),  round(100*sum(IDXp2 == ii)/length(IDXp2))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
end

%%
% %% Cluster each subject to find the Bestk
ktest = 2:15;
nrep=5;
Bestk = zeros(1,M);
dmethod = 'city';

for ii = 1:M
    fprintf('-Subject %d-, k = ', ii)
    subF = squeeze(Fdyn(ii,:,:));
    SSE = zeros(1,length(ktest));
    R = zeros(1,length(ktest)); 
    for k1 = ktest
        fprintf('%d, ', k1)
        [IDX, C, SUMD, D] = kmeans(subF, k1, 'distance', dmethod, 'Replicates', nrep, 'MaxIter', 150, 'empty', 'drop');
        [SSE(k1-1), R(k1-1)] = cluster_goodness(D, IDX);
    end
    [bestx] = fit_L_tocurve_area(ktest,R, 0);
    Bestk(ii) = bestx;
    fprintf('best k = %d\n', bestx)
end

% %save(fullfile(fileparts(inname), 'FNC_subject_bestk'), 'Bestk', 'ktest', 'nrep')
 save(fullfile(fileparts(inname), 'FNC_subject_bestk_scrubbed'), 'Bestk', 'ktest', 'nrep')
 %%
 numPeaks = zeros(M,1);
 for ii = 1:M,
     numPeaks(ii) = size(SP{ii},1);
 end
 cnumPeaks = cumsum(numPeaks);
 numPeaks_HC = numPeaks(find(ica_ursi_group_keepID == 0));
 numPeaks_SZ = numPeaks(find(ica_ursi_group_keepID == 1));
 cnumPeaks_HC = cumsum(numPeaks_HC);
 cnumPeaks_SZ = cumsum(numPeaks_SZ);
 %%
 for k2 = 2:8;
     KoutHC = load(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatHC', 'IDXp1','Cp1', 'HCid');
     
     nHC = length(find(ica_ursi_group_keepID == 0));
     nSZ = length(find(ica_ursi_group_keepID == 1));
     
     HCclustmeans = nan(k2, nHC, size(SPflatHC,2));
     HCclass = cell(nHC,1);
     for ii = 1:nHC
         if ii == 1
             kclass = KoutHC.IDXp1(1:cnumPeaks_HC(ii));
             subjkwin = SPflatHC(1:cnumPeaks_HC(ii),:);
         else
             kclass = KoutHC.IDXp1(cnumPeaks_HC(ii-1)+1:cnumPeaks_HC(ii));
             subjkwin = SPflatHC(cnumPeaks_HC(ii-1)+1:cnumPeaks_HC(ii),:);
         end
         ukclass = unique(kclass);
         HCclass{ii} = ukclass;
         for jj = 1:length(ukclass)
             ck = ukclass(jj);
             subjkmean = mean(subjkwin(kclass == ck,:));
             HCclustmeans(ck,ii,:) = subjkmean;
         end
     end
     
     %%
     KoutSZ = load(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid');
     
     SZclustmeans = nan(k2, nSZ, size(SPflatSZ,2));
     SZclass = cell(nSZ,1);
     for ii = 1:nSZ
         if ii == 1
             kclass = KoutSZ.IDXp2(1:cnumPeaks_SZ(ii));
             subjkwin = SPflatSZ(1:cnumPeaks_SZ(ii),:);
         else
             kclass = KoutSZ.IDXp2(cnumPeaks_SZ(ii-1)+1:cnumPeaks_SZ(ii));
             subjkwin = SPflatSZ(cnumPeaks_SZ(ii-1)+1:cnumPeaks_SZ(ii),:);
         end
         ukclass = unique(kclass);
         SZclass{ii} = ukclass;
         for jj = 1:length(ukclass)
             ck = ukclass(jj);
             subjkmean = mean(subjkwin(kclass == ck,:));
             SZclustmeans(ck,ii,:) = subjkmean;
         end
     end
     %%
     list_k_HC = cell(k2,1);
     for jj = 1:k2
         for ii = 1:nHC,
             
             if sum(ismember(HCclass{ii},jj))
                 list_k_HC{jj} = [list_k_HC{jj} ii];
             end
         end
     end
     list_k_SZ = cell(k2,1);
     for ii = 1:nSZ,
         for jj = 1:k2
             if sum(ismember(SZclass{ii},jj))
                 list_k_SZ{jj} = [list_k_SZ{jj} ii];
             end
         end
     end
     save(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatHC', 'IDXp1','Cp1', 'HCid','HCclustmeans','list_k_HC');
     save(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid','SZclustmeans','list_k_SZ');
 end
%%

cMOD = cumsum(fMOD);
T = '-log_1_0 (P)*sign(t)';
T1 = 'r';
CLIM = [-.5,.5];
CPOSn = [180         886        1717         612];
%%
for k2 = 2:9
    load(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatHC', 'IDXp1','Cp1', 'HCid','HCclustmeans','list_k_HC');
    load(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid','SZclustmeans','list_k_SZ');
    
    G = figure('Color', 'w', 'Name', ['FNC_clustersGroupDiff_k' num2str(k2) '_scrubbed'], 'Position', CPOSn);
    for jj = 1:k2
        [hh pp ci tst] = ttest2(squeeze(HCclustmeans(jj,list_k_HC{jj},:)),squeeze(SZclustmeans(jj,list_k_SZ{jj},:)));
        
     %%
        subplot(3,k2,jj);
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), vec2mat(squeeze(mean(HCclustmeans(jj,list_k_HC{jj},:),2))));
        %imagesc(vec2mat(squeeze(mean(HCclustmeans(jj,list_k_HC{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('HC: (N = %d)',length(list_k_HC{jj})));
        set(gca, 'clim', CLIM); axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(c(2*(length(RSN_I)+1)-cMOD), 'top')
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        axis square
        
        subplot(3,k2,jj+k2);
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), vec2mat(squeeze(mean(SZclustmeans(jj,list_k_SZ{jj},:),2))));
        %imagesc(vec2mat(squeeze(mean(SZclustmeans(jj,list_k_SZ{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('SZ: (N = %d)',length(list_k_SZ{jj})));
        set(gca, 'clim', CLIM); axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(c(2*(length(RSN_I)+1)-cMOD), 'top')
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        axis square
             
        subplot(3,k2,jj+2*k2);
        %imagesc(vec2mat(-log10(pp).*sign(tst.tstat))); 
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I),vec2mat(-log10(pp).*sign(tst.tstat)));
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('HC - SZ; maxP=%1.2g', max(-log10(pp))));
        CLIM1 = [-abs(max(-log10(pp))) abs(max(-log10(pp)))];
        set(gca, 'clim', CLIM1);axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(c(2*(length(RSN_I)+1)-cMOD), 'top')
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T)
        axis square
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
end
%% mancovan
load DEMO_ICA_fBIRNp3_SZ_wPANSS_v3
varnames = DEMO.demo_col_var;

[MODEL] = make_design_matrixSZwPANSSV3(varnames);
%%
for k2 = 2:9
    load(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid','SZclustmeans','list_k_SZ');
    %%
    for jj = 3:k2
        kSZid = list_k_SZ{jj};
        SZkeepid = setdiff(1:151,DEMO.dropSZind);
        [tf loc] = ismember(kSZid,SZkeepid);
        [tf2 loc2] = ismember(SZkeepid,kSZid);
        cSZid = loc(loc > 0);
        kSZid2 = setdiff(kSZid,DEMO.dropSZind);
        data = squeeze(SZclustmeans(jj,kSZid2,:));
        cMODEL = MODEL;
        cMODEL.X = cMODEL.X(cSZid,:);
        [DEMO1, MULT1, UNI1] = run_model_fbirn_SZ_PANSS_V3(cMODEL, 0.01, data , [], [], 1, 1,2);%,40);BIC estimate = 13;

        save([fullfile(fileparts(inname)) filesep 'FNCstats_dynWin_k2_' num2str(k2) '_j_' num2str(jj) '_allcov_scr_SZwPANSS_v3'], 'DEMO1', 'MULT1', 'UNI1');
    end
end
%% Cluster NIC 
stParHC.dEps = 1/size(SPflatHC,1);
stParHC.iNumIter = 10;
stParSZ.dEps = 1/size(SPflatSZ,1);
stParSZ.iNumIter = 10;

for k2 = 2:9
   
    stParHC.iNumClus = k2;
    stParSZ.iNumClus = k2;
    [IDXp1,stResHC] = nic(SPflatHC,stParHC);
    [IDXp2,stResSZ] = nic(SPflatSZ,stParSZ);
   
    save(fullfile(fileparts(inname), ['FNC_HC_clustersNIC_k2_' num2str(k2) '_scrubbed']), 'SPflatHC', 'IDXp1','stResHC', 'HCid');    
    
    save(fullfile(fileparts(inname), ['FNC_SZ_clustersNIC_k2_' num2str(k2) '_scrubbed']), 'SPflatSZ', 'IDXp2','stResHC', 'SZid');

    G = figure('Color', 'w', 'Name', ['FNC_DIAGgroup_clustersNIC_k' num2str(k2) '_scrubbed'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        imagesc(vec2mat(median(SPflatHC(IDXp1 == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp1 == ii),  round(100*sum(IDXp1 == ii)/length(IDXp1))));
        
        subplot(2,k2,ii+k2);
        imagesc(vec2mat(median(SPflatSZ(IDXp2 == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp2 == ii),  round(100*sum(IDXp2 == ii)/length(IDXp2))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
end
%% all data together NIC

SPflatALL = [SPflatHC;SPflatSZ];
stPar.dEps = 1/size(SPflatALL,1);
stPar.iNumIter = 10;
nGrp = [size(SPflatHC,1) size(SPflatSZ,1)];
for k2 = 2:9
   
    stPar.iNumClus = k2;
    
    [IDXp,stRes] = nic(SPflatALL,stPar);    
   
    save(fullfile(fileparts(inname), ['FNC_combined_clustersNICcombined_k' num2str(k2) '_scrubbed']), 'SPflatALL', 'IDXp','stRes', 'nGrp');       
   

    G = figure('Color', 'w', 'Name', ['FNC_DIAGgroup_clustersNICcombined_k' num2str(k2) '_scrubbed'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        imagesc(vec2mat(median(SPflatALL(IDXp(1:nGrp(1)) == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp(1:nGrp(1)) == ii),  round(100*sum(IDXp(1:nGrp(1)) == ii)/length(IDXp))));
        
        subplot(2,k2,ii+k2);
        imagesc(vec2mat(median(SPflatALL(IDXp(nGrp(1)+1:end) == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp(nGrp(1)+1:end) == ii),  round(100*sum(IDXp(nGrp(1)+1:end) == ii)/length(IDXp))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
end

%% single subject example
% ktest = 2:15;
% sub = 297;
% subF = squeeze(Fdyn(sub,:,:));
% SSE = zeros(1,length(ktest));
% R = zeros(1,length(ktest));
% nrep=5;
% for k1 = ktest
%     fprintf('%d, ', k1)
%     [IDX, C, SUMD, D] = kmeans(subF, k1, 'distance', dmethod, 'Replicates', nrep, 'MaxIter', 150, 'empty', 'drop');
%     [SSE(k1-1), R(k1-1)] = cluster_goodness(D, IDX);
% end
% [bestx] = fit_L_tocurve_area(ktest,R, 1);
% set(gcf, 'Name', ['sub' num2str(sub) '_bestk'])
% IM = export_fig(fullfile(FIGoutputdir, [get(gcf, 'Name') '.pdf']), gcf);
% 
% 
% topwin = 10;
% nrep = 50;
% [IDX, C, SUMD, D] = kmeans(subF, bestx, 'distance', dmethod, 'Replicates', nrep, 'MaxIter', 150, 'empty', 'drop');
% F = figure('Color', 'w', 'Name', ['FNC_sub' num2str(sub) 'distance_to_centroid_k' num2str(bestx)]); 
% TR=2;
% timeline = 0:(Nwin-1); timeline = timeline + wsize/2; timeline = timeline*TR;
% plot(D, '.-'); box off; set(gca, 'TickDir', 'out')
% xlabel('Time (s)')
% ylabel('Distance from centroid')
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% CPOS = [ 180         886        1717         212];
% [Y,I] = sort(D,1);
% I = I(1:topwin,:);
% G = figure('Color', 'w', 'Name', ['FNC_sub' num2str(sub) 'centroids_k' num2str(bestx)], 'Position', CPOS);
% for kk = 1:bestx
%     % keep the top windows that are classified in the cluster
%     Ikeep = intersect(I(:,kk), find(IDX == kk));
%     count = count+1;
%     CEN = mean(subF(Ikeep,:));
%     subplot(1,bestx,kk);
%     imagesc(vec2mat(CEN,1), [-.75,.75]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
% end
% IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);


%% Now cluster them at that size
% nrep = 20;
% topwin = 10;
% nCL = sum(Bestk);
% SC = zeros(nCL,size(Fdyn,3));
% count = 0;
% for ii = 1:M
%     fprintf('Subject %d\n', ii)
%     subF = squeeze(Fdyn(ii,:,:));
%     [IDX, C, SUMD, D] = kmeans(subF, Bestk(ii), 'distance', dmethod, 'Replicates', nrep, 'MaxIter', 150, 'empty', 'drop');
%     [Y,I] = sort(D,1);
%     I = I(1:topwin,:);
%     for kk = 1:Bestk(ii)
%         % keep the top windows that are classified in the cluster
%         Ikeep = intersect(I(:,kk), find(IDX == kk));
%         count = count+1;
%         SC(count,:) = mean(subF(Ikeep,:));
%     end
% end
% %save(fullfile(fileparts(inname), 'FNC_subject_clusters_bestk'), 'SC', 'Bestk', 'nrep')
% save(fullfile(fileparts(inname), 'FNC_subject_clusters_bestk_scrubbed'), 'SC', 'Bestk', 'nrep')

%% Determine approximate k2 based on the subject clusters (roughly 2000 samples)
% ktest = 2:10;
% SSE = zeros(1,length(ktest));
% R = zeros(1,length(ktest));
%  nrep=1;
% for k2 = ktest
%     fprintf('---------------- k = %d ------------------\n', k2)
%     [IDX, C, SUMD, D] = kmeans(SC, k2, 'distance', dmethod, 'Replicates', nrep, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop');
%     [SSE(k2-1), R(k2-1)] = cluster_goodness(D, IDX);
% end
% [bestx, OBJ] = fit_L_tocurve(ktest,R, 1);



%% Cluster, based on the subject clusters
% k2 = 7;
% CPOS = [ 180         886        1717         212];
% 
% SCtemp1 =  SCf(size(SCf,1)/2+1:end,:);%SCf(1:size(SCf,1)/2,:);
% [IDXp1, Cp1, SUMDp1, Dp1] = kmeans(SCtemp1, k2, 'distance', dmethod, 'Replicates', 20, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop');
% 
% %save(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2)]), 'SC', 'IDXp','Cp');
% 
% %%
% G = figure('Color', 'w', 'Name', ['FNC_clusters_subclusters_k' num2str(k2)], 'Position', CPOS);
% for ii = 1:k2, 
%     subplot(1,k2,ii);
%     imagesc(vec2mat(median(SCtemp1(IDXp1 == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXp1 == ii),  round(100*sum(IDXp1 == ii)/length(IDXp1))));
% end
% 
% SCtemp2 =  SCf(1:size(SCf,1)/2,:);
% [IDXp2, Cp2, SUMDp2, Dp2] = kmeans(SCtemp2, k2, 'distance', dmethod, 'Replicates', 20, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop');
% 
% %save(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2)]), 'SC', 'IDXp','Cp');
% 
% %%
% G = figure('Color', 'w', 'Name', ['FNC_clusters_subclusters_k' num2str(k2)], 'Position', CPOS);
% for ii = 1:k2, 
%     subplot(1,k2,ii);
%     imagesc(vec2mat(median(SCtemp2(IDXp2 == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXp2 == ii),  round(100*sum(IDXp2 == ii)/length(IDXp2))));
% end
% 
% 
% IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);

%% Use as a starting point to cluster all the data
Fdynflat = reshape(Fdyn, M*Nwin , size(Fdyn,3));
for k2 = [4 5 7 9] 
    
    load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflat', 'IDXp','Cp');
    
    
    [IDXall, Call, SUMDall, Dall] = kmeans(Fdynflat, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
    
    aIND = reshape(IDXall, M, Nwin);
    aIND = aIND';
    
    save(fullfile(fileparts(inname), ['FNC_clusters_all_k' num2str(k2) '_scrubbed']), 'IDXall', 'Call')
    save(fullfile(fileparts(inname), ['FNC_all_states_k' num2str(k2) '_scrubbed']), 'aIND')
    
    G = figure('Color', 'w', 'Name', ['FNC_clusters_all_k' num2str(k2) '_scrubbed'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(1,k2,ii);
        imagesc(vec2mat(median(Fdynflat(IDXall == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXall == ii),  round(100*sum(IDXall == ii)/length(IDXall))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
    
    TR=2;
    timeline = 0:(Nwin-1); timeline = timeline + wsize/2; timeline = timeline*TR;
    
    H = figure('Color', 'w', 'Name', ['FNC_clusters_occurrence_k' num2str(k2) '_scrubbed'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(1,k2,ii);
        for bb = 1:100
            pickME = ceil(M*rand(M,1));
            octemp = 100*mean(aIND(:,pickME) == ii,2);
            plot(timeline, octemp, 'Color', [.8 .8 .8])
            hold on
        end
        hold on
        oc = 100*mean(aIND == ii,2);
        plot(timeline, oc, 'b')
        xlabel('time (s)'); ylabel('frequency'); box off; set(gca, 'TickDir', 'out')
        %set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXall == ii),  round(100*sum(IDXall == ii)/length(IDXall))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(H, 'Name') '.pdf']), H);
end
%%
FdynflatHC = reshape(Fdyn(HCid,:,:), nHC*Nwin , size(Fdyn,3));
FdynflatSZ = reshape(Fdyn(SZid,:,:), nSZ*Nwin , size(Fdyn,3));
for k2 = [4 5 6 7 8 9] 
    %[IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, k2, 'distance', dmethod, 'Replicates', 10, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop');
    %load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) 'scrubbed']), 'SPflat', 'IDXp','Cp');
    
    load(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatHC', 'IDXp1','Cp1', 'HCid');
    load(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_scrubbed']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid');

    [IDXallHC, CallHC, SUMDallHC, DallHC] = kmeans(FdynflatHC, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 60, 'empty', 'drop', 'Start', Cp1);
    [IDXallSZ, CallSZ, SUMDallSZ, DallSZ] = kmeans(FdynflatSZ, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 60, 'empty', 'drop', 'Start', Cp2);

    aINDHC = reshape(IDXallHC, nHC, Nwin);
    aINDHC = aINDHC';
    aINDSZ = reshape(IDXallSZ, nSZ, Nwin);
    aINDSZ = aINDSZ';
    
    save(fullfile(fileparts(inname), ['FNC_clusters_allHC_k' num2str(k2) '_scrubbed_HCcp1']), 'IDXallHC', 'CallHC')
    save(fullfile(fileparts(inname), ['FNC_allHC_states_k' num2str(k2) '_scrubbed_HCcp1']), 'aINDHC')
    save(fullfile(fileparts(inname), ['FNC_clusters_allSZ_k' num2str(k2) '_scrubbed_SZcp1']), 'IDXallSZ', 'CallSZ')
    save(fullfile(fileparts(inname), ['FNC_allSZ_states_k' num2str(k2) '_scrubbed_SZcp1']), 'aINDSZ')

    G = figure('Color', 'w', 'Name', ['FNC_clusters_allHCSZ_k' num2str(k2) '_scrubbed_GRPcp'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        imagesc(vec2mat(median(FdynflatHC(IDXallHC == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXallHC == ii),  round(100*sum(IDXallHC == ii)/length(IDXallHC))));
        
        subplot(2,k2,ii+k2);
        imagesc(vec2mat(median(FdynflatSZ(IDXallSZ == ii ,:)),1), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXallSZ == ii),  round(100*sum(IDXallSZ == ii)/length(IDXallSZ))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
    
    TR=2;
    timeline = 0:(Nwin-1); timeline = timeline + wsize/2; timeline = timeline*TR;
    
    H = figure('Color', 'w', 'Name', ['FNC_clusters_occurrenceHCSZ_k' num2str(k2) '_scrubbed_GRPcp'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        for bb = 1:100
            pickME1 = ceil(nHC*rand(nHC,1));
            octemp = 100*mean(aINDHC(:,pickME1) == ii,2);
            plot(timeline, octemp, 'Color', [.8 .8 .8])
            hold on
        end
        hold on
        oc = 100*mean(aINDHC == ii,2);
        plot(timeline, oc, 'b')
        xlabel('time (s)'); ylabel('frequency'); box off; set(gca, 'TickDir', 'out')
        %set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXallHC == ii),  round(100*sum(IDXallHC == ii)/length(IDXallHC))));
        
        subplot(2,k2,ii+k2);
        for bb = 1:100
            pickME2 = ceil(nSZ*rand(nSZ,1));
            octemp = 100*mean(aINDSZ(:,pickME2) == ii,2);
            plot(timeline, octemp, 'Color', [.8 .8 .8])
            hold on
        end
        hold on
        oc = 100*mean(aINDSZ == ii,2);
        plot(timeline, oc, 'b')
        xlabel('time (s)'); ylabel('frequency'); box off; set(gca, 'TickDir', 'out')
        %set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXallSZ == ii),  round(100*sum(IDXallSZ == ii)/length(IDXallSZ))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(H, 'Name') '.pdf']), H);
end

%% Cluster all the data with random starting point (for comparison) -- makes little difference
% m = 100;
% Fdynflattemp = reshape(Fdyn(1:m,:,:), m*Nwin , size(Fdyn,3));
% 
% k2=6;
% [IDXALL, CALL, SUMDALL, DALL] = kmeans(Fdynflattemp, k2, 'distance', dmethod, 'Replicates', 5, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop');
% 
% G = figure('Color', 'w', 'Name', ['FNC_clusters_k' num2str(k2)], 'Position', CPOS);
% 
% for ii = 1:k2,
%     subplot(1,k2,ii);
%     imagesc(vec2mat(median(Fdynflattemp(IDXALL == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXALL == ii),  round(100*sum(IDXALL == ii)/length(IDXALL))));
% end
% 
% aIND = reshape(IDXALL, M, Nwin);
% aIND = aIND';

% F = figure('Color', 'w', 'Name', 'NumberofClusters_vs_R', 'Position', [680   824   345   274]);
% plot(ktest, R, 'k.-'); box off; set(gca, 'TickDir', 'out')
% xlabel('Number of clusters (k)')
% ylabel('Intra-cluster distance/Extra-cluster distance')
% title(['Distance metric: ' dmethod])
% set(gca, 'Xlim', [ktest(1)-1, ktest(end)+1])
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% 
% %% Cluster based on peaks
% dmethod = 'city';
% k2 = 7;
% CPOS = [ 180         886        1717         212];
% 
% % plot mean
% temp = mean(SP);
% [F,A,C,I] = plot_FNC(temp, [-.5, .5], L, [RSN_I], [], 'average FNC');
% set(F, 'Name', 'FNC_mean_peaks', 'Position', [680          79        1041        1019]);
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% 
% [IDXp, Cp, SUMDp, Dp] = kmeans(dSP, k2, 'distance', dmethod, 'Replicates', 500, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop');
% 
% G = figure('Color', 'w', 'Name', ['FNC_clusters_peak_k' num2str(k2)], 'Position', CPOS);
% 
% for ii = 1:k2,
%     subplot(1,k2,ii);
%     imagesc(vec2mat(median(SP(IDXp == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXp == ii),  round(100*sum(IDXp == ii)/length(IDXp))));
% end
% 
% save(fullfile(fileparts(inname), ['FNC_clusters_500reps_peak_k' num2str(k2)]), 'IDXp', 'Cp', 'dSP', 'SP')
% 
% IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
% 
% %% Cluster based on random selection (equal to or greater number of obs as peaks)
% dSR = zeros(M*nW,size(dFdyn,3)); 
% SR = zeros(M*nW,size(dFdyn,3));
% 
% rng(111); % set the randomness
% for ii = 1:M
%     fprintf('Working on subject %d of %d\n', ii, M)
%     dsubmat = squeeze(dFdyn(ii,:,:));
%     submat = squeeze(Fdyn(ii,:,:));
%     pIND = ceil(Nwin*rand(1,nW));
%     dSR((nW*(ii-1) + 1):nW*ii, : ) = dsubmat(pIND,:);
%     SR((nW*(ii-1) + 1):nW*ii, : ) = submat(pIND,:);
% end
% 
% % plot mean
% temp = mean(SR);
% [F,A,C,I] = plot_FNC(temp, [-.5, .5], L, [RSN_I], [], 'average FNC');
% set(F, 'Name', 'FNC_mean_rand', 'Position', [680          79        1041        1019]);
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% 
% %cluster
% [IDXr, Cr, SUMDr, Dr] = kmeans(dSR, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
% 
% G = figure('Color', 'w', 'Name', ['FNC_clusters_random_k' num2str(k2)], 'Position', CPOS);
% for ii = 1:k2,
%     subplot(1,k2,ii);
%     imagesc(vec2mat(mean(SR(IDXr == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXr == ii),  round(100*sum(IDXr == ii)/length(IDXr))));
% end
% IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
% 
% %% Cluster based on subject clusters
% dSC = zeros(M*nW,size(dFdyn,3)); 
% SC = zeros(M*nW,size(dFdyn,3));
% for ii = 1:M
%     fprintf('Working on subject %d of %d\n', ii, M)
%     dsubmat = squeeze(dFdyn(ii,:,:));
%     submat = squeeze(Fdyn(ii,:,:));
%     [IDXss, Css, SUMDss, Dss] = kmeans(dsubmat, nW, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'empty', 'drop');
%     [~, inds]=min(Dss);
%     dSC((nW*(ii-1) + 1):nW*ii, : ) = dsubmat(inds,:);
%     SC((nW*(ii-1) + 1):nW*ii, : ) = submat(inds,:);
% end
% 
% temp = mean(SC);
% [F,A,C,I] = plot_FNC(temp, [-.5, .5], L, [RSN_I], [], 'average FNC');
% set(F, 'Name', 'FNC_mean_subclust', 'Position', [680          79        1041        1019]);
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% %cluster
% [IDXc, Cc, SUMDc, Dc] = kmeans(dSC, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
% 
% G = figure('Color', 'w', 'Name', ['FNC_clusters_subclust_k' num2str(k2)], 'Position', CPOS);
% for ii = 1:k2,
%     subplot(1,k2,ii);
%     imagesc(vec2mat(mean(SC(IDXc == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXc == ii),  round(100*sum(IDXc == ii)/length(IDXc))));    
% end
% IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
% 
% %% Split subjects in half -- use the top 10 peaks
% 
% dSPh = zeros(M*nW*2,size(dFdyn,3)); % overestimate space needed
% SPh = zeros(M*nW*2,size(dFdyn,3));
% k1_peaksh = zeros(1,M);
% for ii = 1:M
%     fprintf('Working on subject %d of %d\n', ii, M)
%     dsubmat = squeeze(dFdyn(ii,:,:));
%     submat = squeeze(Fdyn(ii,:,:));
%     DEV = sum(abs(dsubmat),2); % sum across bins
%     [xmax,imax,xmin,imin] = extrema(DEV); % find the extrema
%     [xmax, sorted] = sort(xmax, 1, 'descend');
%     imax = imax(sorted);
%     pIND = imax(1:min(nW*2, length(imax)));
%     k1_peaksh(ii) = length(pIND);
%     dSPh((nW*2*(ii-1) + 1):nW*2*(ii-1) + k1_peaksh(ii), : ) = dsubmat(pIND,:);
%     SPh((nW*2*(ii-1) + 1):nW*2*(ii-1) + k1_peaksh(ii), : ) = submat(pIND,:);
% end
% dSPh  = dSPh(sum(dSPh,2) ~= 0,:);
% SPh  =  SPh(sum(SPh,2) ~= 0,:);
% 
% for hh = 1:2
%     if hh == 1
%         sub = floor(M/2);
%         sIND = 1:sum(k1_peaksh(1:sub));
%     else
%         sub = floor(M/2);
%         sIND = (sum(k1_peaksh(1:sub))+1):size(SPh,1);
%     end
%     temp = mean(SPh(sIND,:));
%     [F,A,C,I] = plot_FNC(temp, [-.5, .5], L, [RSN_I], [], 'average FNC');
%     set(F, 'Name', ['FNC_mean_half_h' num2str(hh)], 'Position', [680          79        1041        1019]);
%     IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% end
% 
% 
% for hh = 1:2
%     if hh == 1
%         sub = floor(M/2);
%         sIND = 1:sum(k1_peaksh(1:sub));
%     else
%         sub = floor(M/2);
%         sIND = (sum(k1_peaksh(1:sub))+1):size(SPh,1);
%     end
%     temp = SPh(sIND,:);
%     dtemp = dSPh(sIND,:);
%     [IDXh, Ch, SUMDh, Dh] = kmeans(temp, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
%     G = figure('Color', 'w', 'Name', ['FNC_clusters_half_' num2str(hh) '_k' num2str(k2)], 'Position', CPOS);
%     for ii = 1:k2,
%     subplot(1,k2,ii);
%     imagesc(vec2mat(mean(temp(IDXh == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXh == ii),  round(100*sum(IDXh == ii)/length(IDXh))));    
%     end
%     IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
% end
% 
% %% Now look at surrogate datasets
% % 1. Phases are shifted consistently across components (correlation structure fixed)
% [Sdyn, dSdyn] = create_surrogate_FNCdyn(Fdyn, 'consistent');
% % surrogate, consistent
% % Find the peaks, exactly like in the real data
% dsSP = zeros(M*nW,size(dSdyn,3)); % overestimate space needed
% sSP = zeros(M*nW,size(dSdyn,3));
% sk1_peaks = zeros(1,M);
% for ii = 1:M
%     fprintf('Working on subject %d of %d\n', ii, M)
%     dsubmat = squeeze(dSdyn(ii,:,:));
%     submat = squeeze(Sdyn(ii,:,:));
%     DEV = sum(abs(dsubmat),2); % sum across bins
%     [xmax,imax,xmin,imin] = extrema(DEV); % find the extrema
%     [xmax, sorted] = sort(xmax, 1, 'descend');
%     imax = imax(sorted);
%     pIND = imax(1:min(nW, length(imax)));
%     sk1_peaks(ii) = length(pIND);
%     dsSP((nW*(ii-1) + 1):nW*(ii-1) + sk1_peaks(ii), : ) = dsubmat(pIND,:);
%     sSP((nW*(ii-1) + 1):nW*(ii-1) + sk1_peaks(ii), : ) = submat(pIND,:);
% end
% dsSP  = dsSP(sum(dsSP,2) ~= 0,:);
% sSP  =  sSP(sum(sSP,2) ~= 0,:);
% 
% 
% temp = mean(sSP);
% [F,A,C,I] = plot_FNC(temp, [-.5, .5], L, [RSN_I], [], 'average FNC');
% set(F, 'Name', 'FNC_mean_surr_con', 'Position', [680          79        1041        1019]);
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% 
% % cluster surrogate, consistent
% [IDXsc, Csc, SUMDsc, Dsc] = kmeans(dsSP, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
% 
% G = figure('Color', 'w', 'Name', ['FNC_clusters_surr_con_k' num2str(k2)], 'Position', CPOS);
% for ii = 1:k2,
%     subplot(1,k2,ii);
%     imagesc(vec2mat(mean(sSP(IDXsc == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXsc == ii),  round(100*sum(IDXsc == ii)/length(IDXsc))));    
% end
% IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
% 
% 
% % 2. Phases are shifted inconsistently across components (correlation structure changes)
% [Sdyni, dSdyni] = create_surrogate_FNCdyn(Fdyn, 'inconsistent');
% % Find the peaks, exactly like in the real data
% dsSPi = zeros(M*nW,size(dSdyni,3)); % overestimate space needed
% sSPi = zeros(M*nW,size(dSdyni,3));
% sk1i_peaks = zeros(1,M);
% for ii = 1:M
%     fprintf('Working on subject %d of %d\n', ii, M)
%     dsubmat = squeeze(dSdyni(ii,:,:));
%     submat = squeeze(Sdyni(ii,:,:));
%     DEV = sum(abs(dsubmat),2); % sum across bins
%     [xmax,imax,xmin,imin] = extrema(DEV); % find the extrema
%     [xmax, sorted] = sort(xmax, 1, 'descend');
%     imax = imax(sorted);
%     pIND = imax(1:min(nW, length(imax)));
%     sk1i_peaks(ii) = length(pIND);
%     dsSPi((nW*(ii-1) + 1):nW*(ii-1) + sk1i_peaks(ii), : ) = dsubmat(pIND,:);
%     sSPi((nW*(ii-1) + 1):nW*(ii-1) + sk1i_peaks(ii), : ) = submat(pIND,:);
% end
% dsSPi  = dsSPi(sum(dsSPi,2) ~= 0,:);
% sSPi  =  sSPi(sum(sSPi,2) ~= 0,:);
% 
% 
% % picture of the mean
% temp = mean(sSPi);
% [F,A,C,I] = plot_FNC(temp, [-.5, .5], L, [RSN_I], [], 'average FNC');
% set(F, 'Name', 'FNC_mean_surr_incon', 'Position', [680          79        1041        1019]);
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% %cluster
% [IDXsi, Csi, SUMDsi, Dsi] = kmeans(dsSPi, k2, 'distance', dmethod, 'Replicates', 10, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop');% 'Start', Cp);
% 
% 
% G = figure('Color', 'w', 'Name', ['FNC_clusters_surr_incon_k' num2str(k2)], 'Position', CPOS);
% for ii = 1:k2,
%     subplot(1,k2,ii);
%     imagesc(vec2mat(median(sSPi(IDXsi == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(sprintf('%d (%d%%)', sum(IDXsi == ii),  round(100*sum(IDXsi == ii)/length(IDXsi))));    
% end
% IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);

% %% Evaluate stability of cluster-based classification
% [dLabel, dLabelNull] = cluster_instability(SC, k2, dmethod, 100);
% % should redo this using more data samples
% 
% F = figure('Color', 'w', 'Name', 'Label_Distance', 'Position', [100         948        1106         150]);
% [N,bc] = hist(dLabelNull{1}, 75);
% B = bar(bc, N, 1); box off; 
% set(B, 'EdgeColor', 'none', 'FaceColor', 'k')
% set(gca, 'TickDir', 'out')
% xlabel('Distance between classifiers')
% ylabel('Frequency')
% hold on
% plot(dLabel, 2, 'r*')
% TO = text(dLabel, 50, 'observed'); set(TO, 'Color', 'r', 'HorizontalAlignment', 'center')
% TN = text(mean(dLabelNull{1}), max(N)+50, 'null'); set(TN, 'Color', 'k', 'HorizontalAlignment', 'center');
% axis tight;
% drawaxismargins(gca, .1,0)
% title(sprintf('Z(obs) = %0.1f', (dLabel-mean(dLabelNull{1}))/std(dLabelNull{1})))
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);

%% Hierarchical clustering based on the peaks
% Y = pdist(dSP, 'cosine');
% % %       'single'    --- nearest distance (default)
% % %        'complete'  --- furthest distance
% % %        'average'   --- unweighted average distance (UPGMA) (also known as
% % %                        group average)
% % %        'weighted'  --- weighted average distance (WPGMA)
% % %        'centroid'  --- unweighted center of mass distance (UPGMC)
% % %        'median'    --- weighted center of mass distance (WPGMC)
% % %        'ward'      --- inner squared distance (min variance algorithm)
% Z = linkage(Y, 'ward'); % or 'average', 'ward' 'complete'
% maxc = 6;
% F = figure('Color', 'w', 'Name', 'FNC_HIERARCHY_DENDROGRAM');
% subplot(1,2,1)
% [H,T,perm]=dendrogram(Z, 50, 'COLORTHRESHOLD',28);
% subplot(1,2,2)
% [H,T,perm]=dendrogram(Z, maxc, 'COLORTHRESHOLD',31);
% for cc = 1:maxc
%     n(cc) = sum(T == perm(cc));
% end
% set(gca, 'XTicklabel', num2str(n'))
% 
% G = figure('Color', 'w', 'Name', 'FNC_HIERARCHY_CLUSTERS');
% for mclust = 2:maxc
%     T = cluster(Z, 'maxclust', mclust);
%     kc = length(unique(T));
%     for ii = 1:kc
%     subplot(maxc-1, maxc, maxc*(mclust-2) + ii);
%     imagesc(vec2mat(mean(SP(T == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(num2str(sum(T == ii)));  
%     end
% end
% 
% G = figure('Color', 'w', 'Name', 'FNC_HIERARCHY_CLUSTERS');
% mclust = maxc
%     T = cluster(Z, 'maxclust', mclust);
%     kc = length(unique(T));
%     for ii = 1:kc
%     subplot(1, mclust, ii);
%     imagesc(vec2mat(mean(SP(T == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(num2str(sum(T == ii)));  
%     end

% %% USING EUCLIDEAN DISTANCE
% Y = pdist(dSP, 'euclidean');
% Z = linkage(Y, 'ward'); % or 'complete'
% maxc = 8;
% F = figure('Color', 'w', 'Name', 'FNC_HIERARCHY_DENDROGRAM');
% subplot(1,2,1)
% [H,T,perm]=dendrogram(Z, 50, 'COLORTHRESHOLD',28);
% subplot(1,2,2)
% [H,T,perm]=dendrogram(Z, maxc, 'COLORTHRESHOLD',31);
% for cc = 1:maxc
%     n(cc) = sum(T == perm(cc));
% end
% set(gca, 'XTicklabel', num2str(n'))
% 
% G = figure('Color', 'w', 'Name', 'FNC_HIERARCHY_CLUSTERS');
% for mclust = 2:maxc
%     T = cluster(Z, 'maxclust', mclust);
%     kc = length(unique(T));
%     for ii = 1:kc
%     subplot(maxc-1, maxc, maxc*(mclust-2) + ii);
%     imagesc(vec2mat(mean(SP(T == ii ,:)),1), [-.5,.5]); axis square;
%     set(gca, 'XTick', [], 'YTick', [])
%     title(num2str(sum(T == ii)));  
%     end
% end

%% Network properties of the centrotypes
%                 qtype,  modularity type (see Rubinov and Sporns, 2011)
%                             'sta',  Q_* (default if qtype is not specified)
%                             'pos',  Q_+
%                             'smp',  Q_simple
%                             'gja',  Q_GJA
%                             'neg',  Q_-
% for ii = 1:k2
%     swin = vec2mat(median(Fdynflat(IDXall == ii ,:)),1);
%     swin(isnan(swin)) = 0;
%     [Ci, Q(ii)] = modularity_louvain_und_sign(swin,'sta');
%     nMod(ii) = length(unique(Ci));
% end
% 
% Qs = cell(1,k2);
% nMods = cell(1,k2);
% for ii = 1:k2
%     fprintf('Cluster %d of %d\n', ii, k2)
%     cIND = find(IDXall == ii);
%     for jj = 1:length(cIND)
%         swin = vec2mat(Fdynflat(cIND(jj) ,:));
%         swin(isnan(swin)) = 0;
%         [Ci, Qs{ii}(jj)] = modularity_louvain_und_sign(swin,'sta');
%         nMods{ii}(jj) = length(unique(Ci));
%     end
% end
% %% Example of subject transition statistics
% sub = 297;%22%130;
% F = figure('color', 'w', 'Name', 'SUB297_STATEVECTOR', 'Position', TIMEPOS);
% plot(timeline, aIND(:,sub), 'k'); box off; set(gca, 'TickDir', 'out')
% set(gca, 'XLim', [0, ((Nwin-1) + wsize)*TR])
% set(gca, 'YLim', [0 k2+1], 'YTick', 1:k2)
% xlabel('Time (s)')
% ylabel('State (cluster index)')

% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% [FR, TM, MDT, NT] = statevector_stats(aIND(:,sub), k2);
% G = plot_statevector_stats(k2, FR, TM, MDT, NT);
% set(G, 'Name', 'SUB130_STATESTATS')
% IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);

%% Look at transition statistics for all subjects
aFR = zeros(M,k2);
aTM = zeros(M,k2,k2);
aMDT = zeros(M,k2);
aNT = zeros(M,1);
for ii = 1:M
    [FRii, TMii, MDTii, NTii] = statevector_stats(aIND(:,ii), k2);
    aFR(ii,:) = FRii;
    aTM(ii,:,:) = TMii;
    aMDT(ii,:) = MDTii;
    aNT(ii) = NTii;
end

G = figure('Color', 'w', 'Name', 'Frequency');
edges = linspace(0,1,20);
for ii = 1:k2
    subplot(k2,1,ii);
    N=histc(aFR(:, ii), edges);
    bar(edges,N, 'histc');
    box off;
    set(gca, 'xlim', [0,1], 'ylim', [0 200]);
end
distributionPlot(aFR,2,2,[],1,1,0);
hold on
B1 = boxplot(aFR, 'widths' , .1,  'symbol'  , '', 'whisker', 2);
hold on
A=barplot([1 2 3],mean(G),[mean(G)-std(G)/sqrt(n); mean(G)+std(G)/sqrt(n)], 'o', [.8 .8 .8], 0.6);
drawaxismargins(gca, .05, .045);
ax = axis;
ax(3:4) = [-3.4   15];
axis(ax)
box off
set(gca, 'XTick', 1:3, 'XTickLabel', glabel)
set(gca, 'TickDir', 'out')




G = plot_statevector_stats(k2, aFR, aTM, aMDT, aNT);
set(G, 'Name', 'Average_STATESTATS')
IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);

 A = squeeze(mean(aTM));
 statvec = compute_statprob(A);

%% relationship between transition statistics and demographics
a = loadlabel('age');
g = loadlabel('gender');
IQ = get_IQ;
loga = log(a);
loga = zscore(loga); 


IQall = IQ.zall;
keepIND = ~isnan(IQall);


% with IQ

[ T, p, stats ] = mStepwise(aFR(keepIND,1:end-1), g(keepIND), [loga(keepIND), IQall(keepIND)], 0.05, {'verbose', 'group-covariate','covariate-covariate'});
[ T, p, stats ] = mStepwise(aMDT(keepIND,1:end-1), g(keepIND), [loga(keepIND), IQall(keepIND)], 0.05, {'verbose', 'group-covariate','covariate-covariate'});
%[ T, p, stats ] = mStepwise(aNT(keepIND), g(keepIND), [loga(keepIND), IQall(keepIND)], 0.05, {'verbose', 'group-covariate','covariate-covariate'});

%-- removes gender and IQ terms
% univariate

[ tg, pg, statsg ] = mT(aFR(keepIND,:), g(keepIND), [loga(keepIND), IQall(keepIND)], 1, { 'verbose' });
[ ta, pa, statsa ] = mT(aFR(keepIND,:), g(keepIND), [loga(keepIND), IQall(keepIND)], 2, { 'verbose' });
[ ti, pi, statsi ] = mT(aFR(keepIND,:), g(keepIND), [loga(keepIND), IQall(keepIND)], 3, { 'verbose' });
figure; plot([tg; ta; ti]')
[mv, cIND] = max(ta);
%figure; plot(loga(keepIND), aFR(keepIND,cIND), 'b.')
Y = aFR(keepIND,cIND);

figure; plot(statsg.X(:,3), Y, 'r.')
Yadjusted = Y-statsg.X(:,[2 4])*statsg.B([2 4],cIND);
figure; plot(statsg.X(:,3), Yadjusted, 'r.')


[ tg, pg, statsg ] = mT(aFR(:,:), g(:), [loga(:)], 1, { 'verbose' });
[ ta, pa, statsa ] = mT(aFR(:,:), g(:), [loga(:)], 2, { 'verbose' });



page = [log(10) log(70)];
Bline = statsMAIN.B(1) + page*statsMAIN.B(age_uni_index+1);
Gline = statsMAIN.B(1) + statsMAIN.B(gender_uni_index+1) + page*statsMAIN.B(age_uni_index+1);
hold on
plot(page, Bline, 'b')
plot(page, Gline, 'r')



[ tg, pg, statsg ] = mT(aMDT(keepIND,:), g(keepIND), [loga(keepIND), IQall(keepIND)], 1, { 'verbose' });
[ ta, pa, statsa ] = mT(aMDT(keepIND,:), g(keepIND), [loga(keepIND), IQall(keepIND)], 2, { 'verbose' });
[ ti, pi, statsi ] = mT(aMDT(keepIND,:), g(keepIND), [loga(keepIND), IQall(keepIND)], 3, { 'verbose' });
figure; plot([tg; ta; ti]')

[ tg, pg, statsg ] = mT(aNT(keepIND), g(keepIND), [loga(keepIND), IQall(keepIND)], 1, { 'verbose' });
[ ta, pa, statsa ] = mT(aNT(keepIND), g(keepIND), [loga(keepIND), IQall(keepIND)], 2, { 'verbose' });
[ ti, pi, statsi ] = mT(aNT(keepIND), g(keepIND), [loga(keepIND), IQall(keepIND)], 3, { 'verbose' });





%% age difference - stat prob vector

minage = 17;
maxage = 24;

Af = squeeze(mean(aTM(a < minage,:,:)));
stvecYOUNG = compute_statprob(Af);

Am = squeeze(mean(aTM(a > maxage,:,:)));
stvecOLD = compute_statprob(Am);
gdist = sqrt(sum((stvecYOUNG-stvecOLD).^2));
adist = abs(stvecYOUNG-stvecOLD);
%gtdist = sqrt(sum(sum((Af-Am).^2))); 

aTMr = aTM((a>maxage | a<minage),:,:);
ar = a(a>maxage | a<minage);

nboot = 100;
stvecOLDboot = zeros(nboot,k2);
stvecYOUNGboot = zeros(nboot,k2);
for ii = 1:nboot
    pickme = ceil(length(ar)*rand(1,length(ar)));
    ab = ar(pickme);
    aTMb = aTMr(pickme,:,:);
    Afboot = squeeze(mean(aTMb(ab < minage,:,:)));
    stvecfboot = compute_statprob(Afboot);
    stvecYOUNGboot(ii,:) = stvecfboot;
    Amboot = squeeze(mean(aTMb(ab > maxage,:,:)));
    stvecmboot = compute_statprob(Amboot);
    stvecOLDboot(ii,:) = stvecmboot;
end
otext = sprintf('older subjects (%d to %d, n = %d)', maxage, max(a), sum(a>maxage));
ytext = sprintf('younger subjects (%d to %d, n = %d)', min(a), minage, sum(a<minage));

F = figure('Color', 'w', 'Name', 'StatProbVec_AGE', 'Position', [680   762   617   336]);
Og = plot(1:k2, stvecOLDboot); set(Og, 'Color', [.7 .7 .7])
hold on
Ogl = plot(1:k2, stvecOLD, 'k', 'LineWidth', 2);
Yg = plot(1:k2, stvecYOUNGboot); set(Yg, 'Color', [1 .7 .7])
hold on
Ygl = plot(1:k2, stvecYOUNG, 'r', 'LineWidth', 2);
set(gca, 'TickDir', 'out'); box off
legend([Ogl, Ygl], {otext,ytext}, 'Location', 'Northwest')
xlabel('State (cluster index)')
ylabel('Stationary Probability')
IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);

aTMr = aTM((a>maxage | a<minage),:,:);
ar = a(a>maxage | a<minage);
nperm = 1000;
gdistnull = zeros(1,nboot);
adistnull = zeros(k2,nboot);
for ii = 1:nperm
    aboot = ar(randperm(length(ar)));
    Afboot = squeeze(mean(aTMr(aboot < minage,:,:)));
    stvecfboot = compute_statprob(Afboot);
    Amboot = squeeze(mean(aTMr(aboot > maxage,:,:)));
    stvecmboot = compute_statprob(Amboot);
    gdistnull(ii) = sqrt(sum((stvecfboot-stvecmboot).^2)); 
    adistnull(:,ii) = abs(stvecfboot-stvecmboot);
    %gtdistnull(ii) = sqrt(sum(sum((Afboot-Amboot).^2))); 
end

F = figure('Color', 'w', 'Name', 'StatProbVec_AGE_Distance', 'Position', [100         375        490         300]);
[N,bc] = hist(gdistnull, 75);
B = bar(bc, N, 1); box off; 
set(B, 'EdgeColor', 'none', 'FaceColor', 'k')
set(gca, 'TickDir', 'out')
xlabel('Euclidean distance between state probability vectors')
ylabel('Frequency')
hold on
plot(gdist, 10, 'r*')
TO = text(gdist, 100, 'observed'); set(TO, 'Color', 'r', 'HorizontalAlignment', 'center')
TN = text(mean(gdistnull), max(N)+30, 'null'); set(TN, 'Color', 'k', 'HorizontalAlignment', 'center');
%axis tight;
%drawaxismargins(gca, .1,0)
title(sprintf('Z(obs) = %0.1f, P = %0.3f', (gdist-mean(gdistnull))/std(gdistnull), sum(gdistnull>gdist)/nperm))
IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);



%% distribution of ages
F = figure('Color', 'w', 'Name', 'SubjectAge', 'Position', [680   762   350   250]);
hist(a, min(a):max(a)); box off; set(gca, 'TickDir', 'out'); axis tight
drawaxismargins(gca, .05,0)
xlabel('Age')
ylabel('Frequency')
ax = axis;
hold on
plot([minage, minage], ax(3:4), 'r--')
plot([maxage, maxage], ax(3:4), 'r--')
IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);

%% gender difference/age differences for statvector not found -- come back to this

g_o = loadlabel('gender');
a = loadlabel('age');

aTMr = aTM;
g = g_o;
Af = squeeze(mean(aTMr(g == 1,:,:)));
stvecf = compute_statprob(Af);

Am = squeeze(mean(aTMr(g == -1,:,:)));
stvecm = compute_statprob(Am);
%gdist =  sum(stvecm.*log(stvecm./stvecf));
gdist = sqrt(sum((stvecf-stvecm).^2));

nboot = 1000;
gdistnull = zeros(1,nboot);
for ii = 1:nboot
    gboot = g(randperm(length(g)));
    Afboot = squeeze(mean(aTMr(gboot == 1,:,:)));
    stvecfboot = compute_statprob(Afboot);
    Amboot = squeeze(mean(aTMr(gboot == -1,:,:)));
    stvecmboot = compute_statprob(Amboot);
    gdistnull(ii) = sqrt(sum((stvecfboot-stvecmboot).^2));
    %gdistnull(ii) =  sum(stvecmboot.*log(stvecmboot./stvecfboot));
end