clear all; clc
%%addpath(genpath('/export/mialab/users/eallen/export_fig'))
%inname = '/export/mialab/hcp/dynamics/FNC_dynamicSWTukey18_zscored.mat';
%FIGoutputdir = '/export/mialab/hcp/dynamics/figures/No_L1';
addpath(genpath('/export/mialab/users/eswar/fbirn_p3/scripts/cluster'))
addpath(genpath('/export/mialab/users/eswar/fbirn_p3/scripts/network_estimation'))
addpath(genpath('/export/mialab/users/eswar/software/nic'))
%inname = '/export/mialab/hcp/dynamics/FNC_dynamicSWTukey18_L1ICOV_zscored';
%FIGoutputdir = '/export/mialab/hcp/dynamics/figures/L1';

FIGoutputdir = '/export/mialab/users/eswar/fbirn_p3/results_C100_vn/figures/modularity/dyn_L1';
inname = '/export/mialab/users/eswar/fbirn_p3/results_C100_vn/FNCdynamics/FNC_dynamicSWTukey22_L1ICOV_zscored_scr_cnew';
cd ..
load  orderedRSN_by_modularity
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
ktest = 2:20;dmethod = 'correlation';
%nrep = 150;
for pp = 1:30
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
    
    [bestx(pp), FFF(pp)] = fit_L_tocurve_area(ktest,R, 1);
    drawnow;
end
%% Cluster

% --------------------------------------------------------------------------
% All the data
SPflat = cell2mat(SP);
for k2 =   2:9
[IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop');

save(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat', 'IDXp','Cp');

%load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat', 'IDXp','Cp');

G = figure('Color', 'w', 'Name', ['FNC_group_clusters_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
for ii = 1:k2, 
    subplot(1,k2,ii);
    statemat = vec2mat(median(SPflat(IDXp == ii ,:)));
    imagesc(statemat(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%d (%d%%)', sum(IDXp == ii),  round(100*sum(IDXp == ii)/length(IDXp))));
end
IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
clear IDXp  Cp SUMDp Dp
end
% --------------------------------------------------------------------------

%%
% Half-split analysis
k2 = 5;
sublist = randperm(M);
SPflat1 = cell2mat(SP(sublist(1:floor(M/2))));
SPflat2 = cell2mat(SP(sublist(floor(M/2)+1:end)));
load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat', 'IDXp','Cp');

[IDXp1, Cp1, SUMDp1, Dp1] = kmeans(SPflat1, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop', 'Start', Cp);

save(fullfile(fileparts(inname), ['FNC_half1_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat1', 'IDXp1','Cp1', 'sublist');
%load(fullfile(fileparts(inname), ['FNC_half1_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat1', 'IDXp1','Cp1', 'sublist');

[IDXp2, Cp2, SUMDp2, Dp2] = kmeans(SPflat2, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop', 'Start', Cp);

save(fullfile(fileparts(inname), ['FNC_half2_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat2', 'IDXp2','Cp2', 'sublist');
%load(fullfile(fileparts(inname), ['FNC_half2_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat2', 'IDXp2','Cp2', 'sublist');

G = figure('Color', 'w', 'Name', ['FNC_halfgroup_clusters_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
for ii = 1:k2,
    subplot(2,k2,ii);
    statemat1 = vec2mat(median(SPflat1(IDXp1 == ii ,:)),1);
    imagesc(statemat1(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%d (%d%%)', sum(IDXp1 == ii),  round(100*sum(IDXp1 == ii)/length(IDXp1))));

    subplot(2,k2,ii+k2);
    statemat2 = vec2mat(median(SPflat2(IDXp2 == ii ,:)),1);
    imagesc(statemat2(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%d (%d%%)', sum(IDXp2 == ii),  round(100*sum(IDXp2 == ii)/length(IDXp2))));
end
IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);

%--------------------------------------------------------------------------
%%  HC SZ
%sublist = randperm(M);
HCid = find(ica_ursi_group_keepID == 0);
SZid = find(ica_ursi_group_keepID == 1);
nHC = length(HCid);
nSZ = length(SZid);

SPflatHC = cell2mat(SP(ica_ursi_group_keepID==0));
SPflatSZ = cell2mat(SP(ica_ursi_group_keepID==1));
%%
for k2 = 2:9
    
%     load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat', 'IDXp','Cp')
%     [IDXp1, Cp1, SUMDp1, Dp1] = kmeans(SPflatHC, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop', 'Start', Cp);
% 
%     save(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatHC', 'IDXp1','Cp1', 'HCid');
%     
%     [IDXp2, Cp2, SUMDp2, Dp2] = kmeans(SPflatSZ, k2, 'distance', dmethod, 'Replicates', 1, 'MaxIter', 150, 'Display', 'iter', 'empty', 'drop', 'Start', Cp);
% 
%     save(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid');
    load(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatHC', 'IDXp1','Cp1', 'HCid');
    load(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid');
    G = figure('Color', 'w', 'Name', ['FNC_DIAGgroup_clusters_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        statemat_HC = vec2mat(median(SPflatHC(IDXp1 == ii ,:)),1);
        imagesc(statemat_HC(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp1 == ii),  round(100*sum(IDXp1 == ii)/length(IDXp1))));
        
        subplot(2,k2,ii+k2);
        statemat_SZ = vec2mat(median(SPflatSZ(IDXp2 == ii ,:)),1);
        imagesc(statemat_SZ(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXp2 == ii),  round(100*sum(IDXp2 == ii)/length(IDXp2))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
    clear IDXp1 Cp1  SUMDp1  Dp1 IDXp2 Cp2  SUMDp2  Dp2
end

%%
% %% Cluster each subject to find the Bestk
ktest = 2:10;
nrep=5;
Bestk = zeros(1,M);
dmethod = 'correlation';

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
 save(fullfile(fileparts(inname), 'FNC_subject_bestk_corr_scrubbed_cnew'), 'Bestk', 'ktest', 'nrep')
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
 for k2 = 2:9;
     KoutHC = load(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatHC', 'IDXp1','Cp1', 'HCid');
     
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
     KoutSZ = load(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid');
     
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
     IDXp1 = KoutHC.IDXp1;Cp1 = KoutHC.Cp1;IDXp2 = KoutSZ.IDXp2;Cp2 = KoutSZ.Cp2;
     save(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatHC', 'IDXp1','Cp1', 'HCid','HCclustmeans','list_k_HC');
     save(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid','SZclustmeans','list_k_SZ');
     clear IDXp1 Cp1  SUMDp1  Dp1 IDXp2 Cp2  SUMDp2  Dp2
 end
%%

cMOD = cumsum(mfmod_final);
T = '-log_1_0 (P)*sign(t)';
T1 = 'r';
CLIM = [-.5,.5];
CPOSn = [180         886        1717         612];
%%
for k2 = 2:9
    load(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatHC', 'IDXp1','Cp1', 'HCid','HCclustmeans','list_k_HC');
    load(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid','SZclustmeans','list_k_SZ');
    
    G = figure('Color', 'w', 'Name', ['FNC_clustersGroupDiff_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOSn);
    for jj = 1:k2
        [hh pp ci tst] = ttest2(squeeze(HCclustmeans(jj,list_k_HC{jj},:)),squeeze(SZclustmeans(jj,list_k_SZ{jj},:)));
        
     %%
        subplot(3,k2,jj);
        cmean1 = vec2mat(squeeze(mean(HCclustmeans(jj,list_k_HC{jj},:),2)));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), cmean1(mnOrd_final_ord,mnOrd_final_ord));
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
        cmean2 = vec2mat(squeeze(mean(SZclustmeans(jj,list_k_SZ{jj},:),2)));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), cmean2(mnOrd_final_ord,mnOrd_final_ord));
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
        dmean = vec2mat(-log10(pp).*sign(tst.tstat));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), dmean(mnOrd_final_ord,mnOrd_final_ord));
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
    load(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid','SZclustmeans','list_k_SZ');
    %%
    for jj = 1:k2
        kSZid = list_k_SZ{jj};
        SZkeepid = setdiff(1:151,DEMO.dropSZind);
        [tf loc] = ismember(kSZid,SZkeepid);
        [tf2 loc2] = ismember(SZkeepid,kSZid);
        cSZid = loc(loc > 0);
        kSZid2 = setdiff(kSZid,DEMO.dropSZind);
        data = squeeze(SZclustmeans(jj,kSZid2,:));
        cMODEL = MODEL;
        cMODEL.X = cMODEL.X(cSZid,:);
        mn_CX = mean(cMODEL.X);
        cMODEL.X = cMODEL.X - repmat([0 mn_CX(2) 0 0 0 0 0 0 0 mn_CX(10) mn_CX(11) mn_CX(12)],length(cSZid),1);
        
        [DEMO1, MULT1, UNI1] = run_model_fbirn_SZ_PANSS_V3(cMODEL, 0.01, data , [], [], 1, 1,2);%,40);BIC estimate = 13;

        save([fullfile(fileparts(inname)) filesep 'FNCstats_dynWin_corr_k2_' num2str(k2) '_j_' num2str(jj) '_allcov_scr_SZwPANSS_v3_cnew'], 'DEMO1', 'MULT1', 'UNI1');
    end
end
%%
for k2 = 2:9
    for jj = 1:k2
        
        if exist([fullfile(fileparts(inname)) filesep 'FNCstats_dynWin_corr_k2_' num2str(k2) '_j_' num2str(jj) '_allcov_scr_SZwPANSS_v3_cnew.mat'],'file')
            disp('-------------------------------------------------\n')
            disp(['  k = ' num2str(k2) '; level = ' num2str(jj)])
            
            load([fullfile(fileparts(inname)) filesep 'FNCstats_dynWin_corr_k2_' num2str(k2) '_j_' num2str(jj) '_allcov_scr_SZwPANSS_v3_cnew.mat'])
            disp(MULT1.final_terms)
            disp('-------------------------------------------------\n')
        end
    end
end
            

%%
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
   
    save(fullfile(fileparts(inname), ['FNC_HC_clustersNIC_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatHC', 'IDXp1','stResHC', 'HCid');    
    
    save(fullfile(fileparts(inname), ['FNC_SZ_clustersNIC_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatSZ', 'IDXp2','stResHC', 'SZid');

    G = figure('Color', 'w', 'Name', ['FNC_DIAGgroup_clustersNIC_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
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
   
    save(fullfile(fileparts(inname), ['FNC_combined_clustersNICcombined_k' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatALL', 'IDXp','stRes', 'nGrp');       
   

    G = figure('Color', 'w', 'Name', ['FNC_DIAGgroup_clustersNICcombined_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
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
for k2 = 2:9 
    
%     load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat', 'IDXp','Cp');
%     
%     
%     [IDXall, Call, SUMDall, Dall] = kmeans(Fdynflat, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
%     
%     aIND = reshape(IDXall, M, Nwin);
%     aIND = aIND';
%     
%     save(fullfile(fileparts(inname), ['FNC_clusters_all_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXall', 'Call')
%     save(fullfile(fileparts(inname), ['FNC_all_states_k' num2str(k2) '_corr_scrubbed_cnew']), 'aIND')
    load(fullfile(fileparts(inname), ['FNC_clusters_all_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXall', 'Call')
     load(fullfile(fileparts(inname), ['FNC_all_states_k' num2str(k2) '_corr_scrubbed_cnew']), 'aIND')
    
    G = figure('Color', 'w', 'Name', ['FNC_clusters_all_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(1,k2,ii);
        astatemat = vec2mat(median(Fdynflat(IDXall == ii ,:)),1);
        imagesc(astatemat(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXall == ii),  round(100*sum(IDXall == ii)/length(IDXall))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
    
    TR=2;
    timeline = 0:(Nwin-1); timeline = timeline + wsize/2; timeline = timeline*TR;
    
    H = figure('Color', 'w', 'Name', ['FNC_clusters_occurrence_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
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
%% Use as a starting point to cluster all the data
FdynflatHC = reshape(Fdyn(HCid,:,:), nHC*Nwin , size(Fdyn,3));
FdynflatSZ = reshape(Fdyn(SZid,:,:), nSZ*Nwin , size(Fdyn,3));
for k2 = 2:9
    
%     load(fullfile(fileparts(inname), ['FNC_HC_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatHC', 'IDXp1','Cp1', 'HCid','HCclustmeans','list_k_HC');
%     load(fullfile(fileparts(inname), ['FNC_SZ_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflatSZ', 'IDXp2','Cp2', 'SZid','SZclustmeans','list_k_SZ');
%     
%     
%     [IDXallHC, CallHC, SUMDallHC, DallHC] = kmeans(FdynflatHC, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 60, 'empty', 'drop', 'Start', Cp1);
%     [IDXallSZ, CallSZ, SUMDallSZ, DallSZ] = kmeans(FdynflatSZ, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 60, 'empty', 'drop', 'Start', Cp2);
% 
%     aINDHC = reshape(IDXallHC, nHC, Nwin);
%     aINDHC = aINDHC';
%     aINDSZ = reshape(IDXallSZ, nSZ, Nwin);
%     aINDSZ = aINDSZ';
%     
%     save(fullfile(fileparts(inname), ['FNC_clusters_allHC_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallHC', 'CallHC')
%     save(fullfile(fileparts(inname), ['FNC_allHC_states_k' num2str(k2) '_corr_scrubbed_cnew']), 'aINDHC')
%     save(fullfile(fileparts(inname), ['FNC_clusters_allSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallSZ', 'CallSZ')   
%     save(fullfile(fileparts(inname), ['FNC_allSZ_states_k' num2str(k2) '_corr_scrubbed_cnew']), 'aINDSZ')
     load(fullfile(fileparts(inname), ['FNC_clusters_allHC_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallHC', 'CallHC')
     load(fullfile(fileparts(inname), ['FNC_allHC_states_k' num2str(k2) '_corr_scrubbed_cnew']), 'aINDHC')
     load(fullfile(fileparts(inname), ['FNC_clusters_allSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallSZ', 'CallSZ')   
     load(fullfile(fileparts(inname), ['FNC_allSZ_states_k' num2str(k2) '_corr_scrubbed_cnew']), 'aINDSZ')

    G = figure('Color', 'w', 'Name', ['FNC_clusters_allHCSZ_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        astatemat_HC = vec2mat(median(FdynflatHC(IDXallHC == ii ,:)),1);
        imagesc(astatemat_HC(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXallHC == ii),  round(100*sum(IDXallHC == ii)/length(IDXallHC))));
        
        subplot(2,k2,ii+k2);
        astatemat_SZ = vec2mat(median(FdynflatSZ(IDXallSZ == ii ,:)),1);
        imagesc(astatemat_SZ(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(IDXallSZ == ii),  round(100*sum(IDXallSZ == ii)/length(IDXallSZ))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
      
    TR=2;
    timeline = 0:(Nwin-1); timeline = timeline + wsize/2; timeline = timeline*TR;
    
    H = figure('Color', 'w', 'Name', ['FNC_clusters_occurrenceHCSZ_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
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
    clear IDXallHC  CallHC  SUMDallHC  DallHC IDXallSZ  CallSZ  SUMDallSZ  DallSZ
end

%%
nS_Kwin_allHC = cell(6,1);nS_Kwin_allSZ = cell(6,1);
mean_Kwin_allHC = cell(6,8);mean_Kwin_allSZ = cell(6,8);
mdn_Kwin_allHC = cell(6,8);mdn_Kwin_allSZ = cell(6,8);
for k2 = 3:8
    
    load(fullfile(fileparts(inname), ['FNC_allHC_states_k' num2str(k2) '_corr_scrubbed_cnew']), 'aINDHC')
    load(fullfile(fileparts(inname), ['FNC_allSZ_states_k' num2str(k2) '_corr_scrubbed_cnew']), 'aINDSZ')
    aINDHC = aINDHC';aINDSZ = aINDSZ';
    
    for jj = 1:k2
        %%
        tst1 = find(aINDHC == jj);
        tst2 = find(aINDSZ == jj);
        
        [i11 i12] = ind2sub(size(aINDHC),tst1);
        [i21 i22] = ind2sub(size(aINDSZ),tst2);
        nS_Kwin_allHC{k2-2} = [nS_Kwin_allHC{k2-2} length(unique(i11))];
        nS_Kwin_allSZ{k2-2} = [nS_Kwin_allSZ{k2-2} length(unique(i21))];        
        ksubjHC = unique(i11);ksubjSZ = unique(i21);
     %   mean_Kwin_allHC{k2,jj} = zeros(length(ksubjHC),size(FdynflatHC,2));
        for ii = 1:length(ksubjHC)
            skindHC = sub2ind(size(aINDHC),i11(find(i11 == ksubjHC(ii))),i12(find(i11 == ksubjHC(ii))));
            
            if length(skindHC) > 1
                mean_Kwin_allHC{k2-2,jj} = [mean_Kwin_allHC{k2-2,jj};mean(FdynflatHC(skindHC,:))];
                mdn_Kwin_allHC{k2-2,jj} = [mdn_Kwin_allHC{k2-2,jj};median(FdynflatHC(skindHC,:))];
            else
                mean_Kwin_allHC{k2-2,jj} = [mean_Kwin_allHC{k2-2,jj};FdynflatHC(skindHC,:)];
                mdn_Kwin_allHC{k2-2,jj} = [mdn_Kwin_allHC{k2-2,jj};(FdynflatHC(skindHC,:))];
            end
            clear skindHC
        end
        for ii = 1:length(ksubjSZ)
            skindSZ = sub2ind(size(aINDSZ),i21(find(i21 == ksubjSZ(ii))),i22(find(i21 == ksubjSZ(ii))));
            if length(skindSZ) > 1
                mean_Kwin_allSZ{k2-2,jj} = [mean_Kwin_allSZ{k2-2,jj};mean(FdynflatSZ(skindSZ,:))];mdn_Kwin_allSZ{k2-2,jj} = [mdn_Kwin_allSZ{k2-2,jj};median(FdynflatSZ(skindSZ,:))];
            else
                mean_Kwin_allSZ{k2-2,jj} = [mean_Kwin_allSZ{k2-2,jj};(FdynflatSZ(skindSZ,:))];mdn_Kwin_allSZ{k2-2,jj} = [mdn_Kwin_allSZ{k2-2,jj};(FdynflatSZ(skindSZ,:))];
            end
            clear skindSZ 
        end
        disp(['k2 : ' num2str(k2) '; jj : ' num2str(jj)  '  done..'])
        clear ksubjHC ksubjSZ 
    end
end
save(fullfile(fileparts(inname), ['mn_mdn_FNC_clusters_allByDiag_kall_corr_scrubbed_cnew']), 'mean_Kwin_allHC', 'mean_Kwin_allSZ' , 'mdn_Kwin_allHC', 'mdn_Kwin_allSZ','nS_Kwin_allHC','nS_Kwin_allSZ');
%% stats
load(fullfile(fileparts(inname), ['mn_mdn_FNC_clusters_allByDiag_kall_corr_scrubbed_cnew']), 'mean_Kwin_allHC', 'mean_Kwin_allSZ' , 'mdn_Kwin_allHC', 'mdn_Kwin_allSZ','nS_Kwin_allHC','nS_Kwin_allSZ');
for k2 = 5%3:8
    mx = zeros(1,k2);
    for jj = 1:k2
       [hh pp ci tst] = ttest2(mean_Kwin_allHC{k2-2,jj} ,mean_Kwin_allSZ{k2-2,jj});
       mx(jj) = abs(max(-log10(pp)));
       [hh1 pp1 ci1 tst1] = ttest2(mdn_Kwin_allHC{k2-2,jj} ,mdn_Kwin_allSZ{k2-2,jj});
       mdx(jj) = abs(max(-log10(pp1)));
       
       clear hh pp ci tst hh1 pp1 ci1 tst1
    end
    mmx = max(mx);
    mmdx = max(mdx);
    
    G1 = figure('Color', 'w', 'Name', ['FNC_clustersMNGroupDiff_allSubj_k' num2str(k2) '_corr_scrubbed_cnewN'], 'Position', CPOSn);
    for jj = 1:k2
        [hh pp ci tst] = ttest2(mean_Kwin_allHC{k2-2,jj} ,mean_Kwin_allSZ{k2-2,jj});
        
     %%
        subplot(3,k2,jj);
        cmat1 = vec2mat(squeeze(mean(mean_Kwin_allHC{k2-2,jj})));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I),cmat1(mnOrd_final_ord,mnOrd_final_ord) );
        %imagesc(vec2mat(squeeze(mean(HCclustmeans(jj,list_k_HC{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('HC: (N = %d)',nS_Kwin_allHC{k2-2}(jj)));
        set(gca, 'clim', CLIM); axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        axis square
        
        subplot(3,k2,jj+k2);
        cmat2 = vec2mat(squeeze(mean(mean_Kwin_allSZ{k2-2,jj})));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), cmat2(mnOrd_final_ord,mnOrd_final_ord) );
        %imagesc(vec2mat(squeeze(mean(SZclustmeans(jj,list_k_SZ{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('SZ: (N = %d)',nS_Kwin_allSZ{k2-2}(jj)));
        set(gca, 'clim', CLIM); axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')

        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        axis square
             
        subplot(3,k2,jj+2*k2);
        dmat = vec2mat(-log10(pp).*sign(tst.tstat));
        %imagesc(vec2mat(-log10(pp).*sign(tst.tstat))); 
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), dmat(mnOrd_final_ord,mnOrd_final_ord));
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('HC - SZ; maxP=%1.2g', max(-log10(pp))));
        CLIM1 = [-mmx mmx];
        set(gca, 'clim', CLIM1);axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
         set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(c(2*(length(RSN_I)+1)-cMOD), 'top')
        A = gca;        
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T)
        axis square
        p_masked = fnctb_fdr(pp,0.01);
        disp(p_masked)
        if p_masked == 0
            fdrlim = 0;
        else
            fdrlim =  -log10(p_masked);%-log10(truep{cc}(p_masked));
        end
        yt = get(C, 'YTick');
        if fdrlim > 0
            set(C, 'YTick', sort([yt, -fdrlim fdrlim]))
        end
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G1, 'Name') '.pdf']), G1);
  %  close(G1)
    
    G2 = figure('Color', 'w', 'Name', ['FNC_clustersMDNGroupDiff_allSubj_k' num2str(k2) '_corr_scrubbed_cnewN'], 'Position', CPOSn);
    for jj = 1:k2
        [hh pp ci tst] = ttest2(mdn_Kwin_allHC{k2-2,jj} ,mdn_Kwin_allSZ{k2-2,jj});
        
     %%
        subplot(3,k2,jj);
        cdmat1 = vec2mat(squeeze(mean(mdn_Kwin_allHC{k2-2,jj})));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), cdmat1(mnOrd_final_ord,mnOrd_final_ord));
        %imagesc(vec2mat(squeeze(mean(HCclustmeans(jj,list_k_HC{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('HC: (N = %d)',nS_Kwin_allHC{k2-2}(jj)));
        axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        cc = c(find(strcmp(get(c, 'Type'),'line')));
        set(cc(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(cc(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(cc(2*(length(RSN_I)+1)-cMOD), 'top')
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        set(gca, 'clim', CLIM); 
        
        
        subplot(3,k2,jj+k2);
        cdmat2 = vec2mat(squeeze(mean(mdn_Kwin_allSZ{k2-2,jj})));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), cdmat2(mnOrd_final_ord,mnOrd_final_ord));
        %imagesc(vec2mat(squeeze(mean(SZclustmeans(jj,list_k_SZ{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('SZ: (N = %d)',nS_Kwin_allSZ{k2-2}(jj)));
        set(gca, 'clim', CLIM); axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        cc = c(find(strcmp(get(c, 'Type'),'line')));
        set(cc(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(cc(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(cc(2*(length(RSN_I)+1)-cMOD), 'top')
        
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        
             
        subplot(3,k2,jj+2*k2);
        ddmat = vec2mat(-log10(pp).*sign(tst.tstat));
        
        %imagesc(vec2mat(-log10(pp).*sign(tst.tstat))); 
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I),ddmat(mnOrd_final_ord,mnOrd_final_ord));
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('HC - SZ; maxP=%1.2g', max(-log10(pp))));
        CLIM1 = [-mmdx mmdx];
        set(gca, 'clim', CLIM1);axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        cc = c(find(strcmp(get(c, 'Type'),'line')));
        set(cc(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(cc(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(cc(2*(length(RSN_I)+1)-cMOD), 'top')
        p_masked = fnctb_fdr(pp,0.01);
        disp(p_masked)
        if p_masked == 0
            fdrlim = 0;
        else
            fdrlim =  -log10(p_masked);%-log10(truep{cc}(p_masked));
        end
        
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T)
        axis square
        yt = get(C, 'YTick');
        if fdrlim > 0
            set(C, 'YTick', sort([yt, -fdrlim fdrlim]))
        end
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G2, 'Name') '.pdf']), G2);
    
   % close(G2)
end

%% Cluster all data using single kmean algorithm for HC and SZ

Fdynflat = reshape(Fdyn, M*Nwin , size(Fdyn,3));
%%
for k2 = 3:9
    %load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) '_corr_scrubbed_cnew']), 'SPflat', 'IDXp','Cp');
        
    %[IDXallone, Callone, SUMDallone, Dallone] = kmeans(Fdynflat, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
    

    %aINDallone = reshape(IDXallone, M,Nwin);    
    
    %save(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallone','aINDallone', 'Callone','SUMDallone','Dallone')
    load(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallone','aINDallone', 'Callone','SUMDallone','Dallone')
    aINDallone_SZ = aINDallone(SZid,:);
    aINDallone_HC = aINDallone(HCid,:);
    

    G = figure('Color', 'w', 'Name', ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        a1_statemat_HC = vec2mat(median(FdynflatHC(aINDallone_HC == ii ,:)),1);
        imagesc(a1_statemat_HC(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(sum(aINDallone_HC == ii)),  round(100*sum(sum(aINDallone_HC == ii))/length(aINDallone_HC(:)))));
        
        subplot(2,k2,ii+k2);
        a1_statemat_SZ = vec2mat(median(FdynflatSZ(aINDallone_SZ == ii ,:)),1);
        imagesc(a1_statemat_SZ(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(sum(aINDallone_SZ == ii)),  round(100*sum(sum(aINDallone_SZ == ii))/length(aINDallone_SZ(:)))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
      
    TR=2;
    timeline = 0:(Nwin-1); timeline = timeline + wsize/2; timeline = timeline*TR;
    aINDallone_HC = aINDallone_HC';
    aINDallone_SZ = aINDallone_SZ';
    H = figure('Color', 'w', 'Name', ['FNC_clusters_allone_occurrenceHCSZ_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        for bb = 1:100
            pickME1 = ceil(nHC*rand(nHC,1));
            octemp = 100*mean(aINDallone_HC(:,pickME1) == ii,2);
            plot(timeline, octemp, 'Color', [.8 .8 .8])
            hold on
        end
        hold on
        oc = 100*mean(aINDallone_HC == ii,2);
        plot(timeline, oc, 'b')
        xlabel('time (s)'); ylabel('frequency'); box off; set(gca, 'TickDir', 'out')
        %set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(aINDallone_HC(:) == ii),  round(100*sum(aINDallone_HC(:) == ii)/length(aINDallone_HC(:)))));
        
        subplot(2,k2,ii+k2);
        for bb = 1:100
            pickME2 = ceil(nSZ*rand(nSZ,1));
            octemp = 100*mean(aINDallone_SZ(:,pickME2) == ii,2);
            plot(timeline, octemp, 'Color', [.8 .8 .8])
            hold on
        end
        hold on
        oc = 100*mean(aINDallone_SZ == ii,2);
        plot(timeline, oc, 'b')
        xlabel('time (s)'); ylabel('frequency'); box off; set(gca, 'TickDir', 'out')
        %set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(aINDallone_SZ(:) == ii),  round(100*sum(aINDallone_SZ(:) == ii)/length(aINDallone_SZ(:)))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(H, 'Name') '.pdf']), H);
    clear IDXallone Callone SUMDallone Dallone IDXallone_HC  IDXallone_SZ  
end
%% stats allone
clear nS_Kwin_allHC  nS_Kwin_allSZ mean_Kwin_allHC mean_Kwin_allSZ mdn_Kwin_allHC mdn_Kwin_allSZ
nS_Kwin_allHC = cell(6,1);nS_Kwin_allSZ = cell(6,1);
mean_Kwin_allHC = cell(6,8);mean_Kwin_allSZ = cell(6,8);
mdn_Kwin_allHC = cell(6,8);mdn_Kwin_allSZ = cell(6,8);
for k2 = 3:8   
    
    load(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallone','aINDallone', 'Callone','SUMDallone','Dallone')
    aINDallone_SZ = aINDallone(SZid,:);
    aINDallone_HC = aINDallone(HCid,:);
    aINDHC = aINDallone_HC;aINDSZ = aINDallone_SZ;
    
    for jj = 1:k2
        %%
        tst1 = find(aINDHC == jj);
        tst2 = find(aINDSZ == jj);
        
        [i11 i12] = ind2sub(size(aINDHC),tst1);
        [i21 i22] = ind2sub(size(aINDSZ),tst2);
        nS_Kwin_allHC{k2-2} = [nS_Kwin_allHC{k2-2} length(unique(i11))];
        nS_Kwin_allSZ{k2-2} = [nS_Kwin_allSZ{k2-2} length(unique(i21))];        
        ksubjHC = unique(i11);ksubjSZ = unique(i21);
     %   mean_Kwin_allHC{k2,jj} = zeros(length(ksubjHC),size(FdynflatHC,2));
        for ii = 1:length(ksubjHC)
            skindHC = sub2ind(size(aINDHC),i11(find(i11 == ksubjHC(ii))),i12(find(i11 == ksubjHC(ii))));
            
            if length(skindHC) > 1
                mean_Kwin_allHC{k2-2,jj} = [mean_Kwin_allHC{k2-2,jj};mean(FdynflatHC(skindHC,:))];
                mdn_Kwin_allHC{k2-2,jj} = [mdn_Kwin_allHC{k2-2,jj};median(FdynflatHC(skindHC,:))];
            else
                mean_Kwin_allHC{k2-2,jj} = [mean_Kwin_allHC{k2-2,jj};FdynflatHC(skindHC,:)];
                mdn_Kwin_allHC{k2-2,jj} = [mdn_Kwin_allHC{k2-2,jj};(FdynflatHC(skindHC,:))];
            end
            clear skindHC
        end
        for ii = 1:length(ksubjSZ)
            skindSZ = sub2ind(size(aINDSZ),i21(find(i21 == ksubjSZ(ii))),i22(find(i21 == ksubjSZ(ii))));
            if length(skindSZ) > 1
                mean_Kwin_allSZ{k2-2,jj} = [mean_Kwin_allSZ{k2-2,jj};mean(FdynflatSZ(skindSZ,:))];mdn_Kwin_allSZ{k2-2,jj} = [mdn_Kwin_allSZ{k2-2,jj};median(FdynflatSZ(skindSZ,:))];
            else
                mean_Kwin_allSZ{k2-2,jj} = [mean_Kwin_allSZ{k2-2,jj};(FdynflatSZ(skindSZ,:))];mdn_Kwin_allSZ{k2-2,jj} = [mdn_Kwin_allSZ{k2-2,jj};(FdynflatSZ(skindSZ,:))];
            end
            clear skindSZ 
        end
        disp(['k2 : ' num2str(k2) '; jj : ' num2str(jj)  '  done..'])
        clear ksubjHC ksubjSZ 
    end
end
save(fullfile(fileparts(inname), ['mn_mdn_FNC_clusters_allone_ByDiag_kall_corr_scrubbed_cnew']), 'mean_Kwin_allHC', 'mean_Kwin_allSZ' , 'mdn_Kwin_allHC', 'mdn_Kwin_allSZ','nS_Kwin_allHC','nS_Kwin_allSZ');
%%
load(fullfile(fileparts(inname), ['mn_mdn_FNC_clusters_allone_ByDiag_kall_corr_scrubbed_cnew']), 'mean_Kwin_allHC', 'mean_Kwin_allSZ' , 'mdn_Kwin_allHC', 'mdn_Kwin_allSZ','nS_Kwin_allHC','nS_Kwin_allSZ');
for k2 = 5%3:8
    mx = zeros(1,k2);
    for jj = 1:k2
       [hh pp ci tst] = ttest2(mean_Kwin_allHC{k2-2,jj} ,mean_Kwin_allSZ{k2-2,jj});
       mx(jj) = abs(max(-log10(pp)));
       [hh1 pp1 ci1 tst1] = ttest2(mdn_Kwin_allHC{k2-2,jj} ,mdn_Kwin_allSZ{k2-2,jj});
       mdx(jj) = abs(max(-log10(pp1)));
       
       clear hh pp ci tst hh1 pp1 ci1 tst1
    end
    mmx = max(mx);
    mmdx = max(mdx);
    
    G1 = figure('Color', 'w', 'Name', ['FNC_clustersMNGroupDiff_allone_Subj_k' num2str(k2) '_corr_scrubbed_cnew_orig'], 'Position', CPOSn);
    for jj = 1:k2
        [hh pp ci tst] = ttest2(mean_Kwin_allHC{k2-2,jj} ,mean_Kwin_allSZ{k2-2,jj});
        
     %%
        subplot(3,k2,jj);
        cmat1 = vec2mat(squeeze(mean(mean_Kwin_allHC{k2-2,jj})));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I),cmat1(mnOrd_final_ord,mnOrd_final_ord) );
        %imagesc(vec2mat(squeeze(mean(HCclustmeans(jj,list_k_HC{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('HC: (N = %d)',nS_Kwin_allHC{k2-2}(jj)));
        set(gca, 'clim', CLIM); axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        axis square
        
        subplot(3,k2,jj+k2);
        cmat2 = vec2mat(squeeze(mean(mean_Kwin_allSZ{k2-2,jj})));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), cmat2(mnOrd_final_ord,mnOrd_final_ord) );
        %imagesc(vec2mat(squeeze(mean(SZclustmeans(jj,list_k_SZ{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('SZ: (N = %d)',nS_Kwin_allSZ{k2-2}(jj)));
        set(gca, 'clim', CLIM); axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')

        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        axis square
             
        subplot(3,k2,jj+2*k2);
        dmat = vec2mat(log10(pp).*sign(tst.tstat));
        %imagesc(vec2mat(-log10(pp).*sign(tst.tstat))); 
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), dmat(mnOrd_final_ord,mnOrd_final_ord));
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('SZ - HC; maxP=%1.2g', max(-log10(pp))));
        CLIM1 = [-mmx mmx];
        set(gca, 'clim', CLIM1);axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
         set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        set(c(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(c(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(c(2*(length(RSN_I)+1)-cMOD), 'top')
        A = gca;        
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T)
        axis square
        p_masked = fnctb_fdr(pp,0.01);
        disp(p_masked)
        if p_masked == 0
            fdrlim = 0;
        else
            fdrlim =  -log10(p_masked);%-log10(truep{cc}(p_masked));
        end
        yt = get(C, 'YTick');
        if fdrlim > 0
            set(C, 'YTick', sort([yt, -fdrlim fdrlim]))
        end
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G1, 'Name') '.pdf']), G1);
  %  close(G1)
    
    G2 = figure('Color', 'w', 'Name', ['FNC_clustersMDNGroupDiff_allone_Subj_k' num2str(k2) '_corr_scrubbed_cnew_orig'], 'Position', CPOSn);
    for jj = 1:k2
        [hh pp ci tst] = ttest2(mdn_Kwin_allHC{k2-2,jj} ,mdn_Kwin_allSZ{k2-2,jj});
        
     %%
        subplot(3,k2,jj);
        cdmat1 = vec2mat(squeeze(mean(mdn_Kwin_allHC{k2-2,jj})));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), cdmat1(mnOrd_final_ord,mnOrd_final_ord));
        %imagesc(vec2mat(squeeze(mean(HCclustmeans(jj,list_k_HC{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('HC: (N = %d)',nS_Kwin_allHC{k2-2}(jj)));
        axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        cc = c(find(strcmp(get(c, 'Type'),'line')));
        set(cc(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(cc(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(cc(2*(length(RSN_I)+1)-cMOD), 'top')
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        set(gca, 'clim', CLIM); 
        
        
        subplot(3,k2,jj+k2);
        cdmat2 = vec2mat(squeeze(mean(mdn_Kwin_allSZ{k2-2,jj})));
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I), cdmat2(mnOrd_final_ord,mnOrd_final_ord));
        %imagesc(vec2mat(squeeze(mean(SZclustmeans(jj,list_k_SZ{jj},:),2))), [-.5,.5]); 
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('SZ: (N = %d)',nS_Kwin_allSZ{k2-2}(jj)));
        set(gca, 'clim', CLIM); axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        cc = c(find(strcmp(get(c, 'Type'),'line')));
        set(cc(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(cc(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(cc(2*(length(RSN_I)+1)-cMOD), 'top')
        
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T1)
        
             
        subplot(3,k2,jj+2*k2);
        ddmat = vec2mat(log10(pp).*sign(tst.tstat));
        
        %imagesc(vec2mat(-log10(pp).*sign(tst.tstat))); 
        simtb_pcolor(1:length(RSN_I), 1:length(RSN_I),ddmat(mnOrd_final_ord,mnOrd_final_ord));
        axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('SZ - HC; maxP=%1.2g', max(-log10(pp))));
        CLIM1 = [-mmdx mmdx];
        set(gca, 'clim', CLIM1);axis ij;
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
        cc = c(find(strcmp(get(c, 'Type'),'line')));
        set(cc(length(RSN_I)-cMOD+1), 'Color', 'k')
        set(cc(2*(length(RSN_I)+1)-cMOD), 'Color', 'k')
        uistack(cc(2*(length(RSN_I)+1)-cMOD), 'top')
        p_masked = fnctb_fdr(pp,0.01);
        disp(p_masked)
        if p_masked == 0
            fdrlim = 0;
        else
            fdrlim =  -log10(p_masked);%-log10(truep{cc}(p_masked));
        end
        
        A = gca;
        I = c(find(strcmp(get(c, 'Type'),'image')));
        C = colorbar; 
        set(get(C, 'YLabel'), 'String', T)
        axis square
        yt = get(C, 'YTick');
        if fdrlim > 0
            set(C, 'YTick', sort([yt, -fdrlim fdrlim]))
        end
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G2, 'Name') '.pdf']), G2);
    
   % close(G2)
end
%%
%%  subcortical (SC) antagonism  vs sensorimotor (SM) hyperconnectivity
% SM : all sensory ICNs (auditory, visual, motor)
load(fullfile(fileparts(inname), 'mn_mdn_FNC_clusters_allone_ByDiag_kall_corr_scrubbed_cnew'))
k2 = 5;
sc_to_sm_corr_HC_allstates = cell(1,k2);
scnot_to_sm_corr_HC_allstates = cell(1,k2);
within_sc_corr_HC_allstates = cell(1,k2);
within_sm_corr_HC_allstates = cell(1,k2);
thal_to_sm_corr_HC_allstates = cell(1,k2);
sc_to_sm_corr_SZ_allstates = cell(1,k2);
scnot_to_sm_corr_SZ_allstates = cell(1,k2);
within_sc_corr_SZ_allstates = cell(1,k2);
within_sm_corr_SZ_allstates = cell(1,k2);
thal_to_sm_corr_SZ_allstates = cell(1,k2);

% r_SC2SMcorr_wSMcorr_HC_obs = zeros(1,k2); p_SC2SMcorr_wSMcorr_HC_obs = ones(1,k2);
% r_THAL2SMcorr_wSMcorr_HC_obs = zeros(1,k2); p_THAL2SMcorr_wSMcorr_HC_obs = ones(1,k2);
% r_wSCcorr_wSMcorr_HC_obs  = zeros(1,k2); p_wSCcorr_wSMcorr_HC_obs = ones(1,k2);
% 
% r_SC2SMcorr_wSMcorr_SZ_obs = zeros(1,k2); p_SC2SMcorr_wSMcorr_SZ_obs = ones(1,k2);
% r_THAL2SMcorr_wSMcorr_SZ_obs = zeros(1,k2); p_THAL2SMcorr_wSMcorr_SZ_obs = ones(1,k2);
% r_wSCcorr_wSMcorr_SZ_obs  = zeros(1,k2); p_wSCcorr_wSMcorr_SZ_obs = ones(1,k2);
% 
% nBOOT = 10000;
% r_SC2SMcorr_wSMcorr_HC_boot = zeros(nBOOT,k2);
% r_THAL2SMcorr_wSMcorr_HC_boot = zeros(nBOOT,k2);
% r_wSCcorr_wSMcorr_HC_boot = zeros(nBOOT,k2);
% r_SC2SMcorr_wSMcorr_SZ_boot = zeros(nBOOT,k2);
% r_THAL2SMcorr_wSMcorr_SZ_boot = zeros(nBOOT,k2);
% r_wSCcorr_wSMcorr_SZ_boot = zeros(nBOOT,k2);


for state = 1:k2,
    
    mdnHC_state = vec2mat(mdn_Kwin_allHC{k2-2,state});
    mdnHC_state = mdnHC_state(:,mnOrd_final_ord,mnOrd_final_ord);
    %figure;imagesc(squeeze(mean(mdnHC_state)),[-0.5 0.5])

    mdnSZ_state = vec2mat(mdn_Kwin_allSZ{k2-2,state});
    mdnSZ_state = mdnSZ_state(:,mnOrd_final_ord,mnOrd_final_ord);
    
    sc_to_sm_corr_HC = zeros(nS_Kwin_allHC{k2-2}(state),1);
    scnot_to_sm_corr_HC = zeros(nS_Kwin_allHC{k2-2}(state),1);
    within_sc_corr_HC = sc_to_sm_corr_HC;within_sm_corr_HC = sc_to_sm_corr_HC;
    thal_to_sm_corr_HC = sc_to_sm_corr_HC;
    for ii = 1:nS_Kwin_allHC{k2-2}(state)
        sc_to_sm_corr_HC(ii) = sum(sum(squeeze(mdnHC_state(ii,1:5,6:24))))/numel(squeeze(mdnHC_state(ii,1:5,6:24))); 
        scnot_to_sm_corr_HC(ii) = sum(sum(squeeze(mdnHC_state(ii,1:4,6:24))))/numel(squeeze(mdnHC_state(ii,1:4,6:24))); 
        within_sc_corr_HC(ii) = 0.5*(sum(sum(squeeze(mdnHC_state(ii,1:5,1:5)))))/nchoosek(5,2);
        within_sm_corr_HC(ii) = 0.5*(sum(sum(squeeze(mdnHC_state(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
        thal_to_sm_corr_HC(ii) = (sum(squeeze(mdnHC_state(ii,5,6:24))))/numel(6:24); 
    end
    
    sc_to_sm_corr_SZ = zeros(nS_Kwin_allSZ{k2-2}(state),1);
    scnot_to_sm_corr_SZ = zeros(nS_Kwin_allSZ{k2-2}(state),1);
    within_sc_corr_SZ = sc_to_sm_corr_SZ;within_sm_corr_SZ = sc_to_sm_corr_SZ;
    thal_to_sm_corr_SZ = sc_to_sm_corr_SZ;
    for ii = 1:nS_Kwin_allSZ{k2-2}(state)
        sc_to_sm_corr_SZ(ii) = sum(sum(squeeze(mdnSZ_state(ii,1:5,6:24))))/numel(squeeze(mdnSZ_state(ii,1:5,6:24))); 
        scnot_to_sm_corr_SZ(ii) = sum(sum(squeeze(mdnSZ_state(ii,1:4,6:24))))/numel(squeeze(mdnSZ_state(ii,1:4,6:24))); 
        within_sc_corr_SZ(ii) = 0.5*(sum(sum(squeeze(mdnSZ_state(ii,1:5,1:5)))))/nchoosek(5,2);
        within_sm_corr_SZ(ii) = 0.5*(sum(sum(squeeze(mdnSZ_state(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
        thal_to_sm_corr_SZ(ii) = (sum(squeeze(mdnSZ_state(ii,5,6:24))))/numel(6:24); 
    end
    
    sc_to_sm_corr_HC_allstates{state} = sc_to_sm_corr_HC;
    scnot_to_sm_corr_HC_allstates{state} = scnot_to_sm_corr_HC;
    within_sc_corr_HC_allstates{state} = within_sc_corr_HC;
    within_sm_corr_HC_allstates{state} = within_sm_corr_HC;
    thal_to_sm_corr_HC_allstates{state} = thal_to_sm_corr_HC;
    
    sc_to_sm_corr_SZ_allstates{state} = sc_to_sm_corr_SZ;
    scnot_to_sm_corr_SZ_allstates{state} = scnot_to_sm_corr_SZ;
    within_sc_corr_SZ_allstates{state} = within_sc_corr_SZ;
    within_sm_corr_SZ_allstates{state} = within_sm_corr_SZ;
    thal_to_sm_corr_SZ_allstates{state} = thal_to_sm_corr_SZ;
    
%     [r_SC2SMcorr_wSMcorr_HC_obs(state) p_SC2SMcorr_wSMcorr_HC_obs(state)] =  corr(sc_to_sm_corr_HC_allstates{state}, within_sm_corr_HC_allstates{state});
%     [r_THAL2SMcorr_wSMcorr_HC_obs(state) p_THAL2SMcorr_wSMcorr_HC_obs(state)] =  corr(thal_to_sm_corr_HC_allstates{state}, within_sm_corr_HC_allstates{state});
%     [r_wSCcorr_wSMcorr_HC_obs(state) p_wSCcorr_wSMcorr_HC_obs(state)] =  corr(within_sc_corr_HC_allstates{state}, within_sm_corr_HC_allstates{state});
%     
%     [r_SC2SMcorr_wSMcorr_SZ_obs(state) p_SC2SMcorr_wSMcorr_SZ_obs(state)] =  corr(sc_to_sm_corr_SZ_allstates{state}, within_sm_corr_SZ_allstates{state});
%     [r_THAL2SMcorr_wSMcorr_SZ_obs(state) p_THAL2SMcorr_wSMcorr_SZ_obs(state)] =  corr(thal_to_sm_corr_SZ_allstates{state}, within_sm_corr_SZ_allstates{state});
%     [r_wSCcorr_wSMcorr_SZ_obs(state) p_wSCcorr_wSMcorr_SZ_obs(state)] =  corr(within_sc_corr_SZ_allstates{state}, within_sm_corr_SZ_allstates{state});
% 
%     for bb = 1:nBOOT,
%         if mod(bb,100) == 1
%             disp(['now working on ' num2str(bb) ' of ' num2str(nBOOT) ' : state ' num2str(state)])
%         end
%         bootHCind = ceil(rand(1,nS_Kwin_allHC{k2-2}(state))*nS_Kwin_allHC{k2-2}(state));        
%         bootSZind = ceil(rand(1,nS_Kwin_allSZ{k2-2}(state))*nS_Kwin_allSZ{k2-2}(state));
%         
%         mdnHC_state = vec2mat(mdn_Kwin_allHC{k2-2,state});
%         mdnHC_state = mdnHC_state(:,mnOrd_final_ord,mnOrd_final_ord);
%         mdnHC_state = mdnHC_state(bootHCind,:,:);
%         mdnSZ_state = vec2mat(mdn_Kwin_allSZ{k2-2,state});
%         mdnSZ_state = mdnSZ_state(:,mnOrd_final_ord,mnOrd_final_ord);
%         mdnSZ_state = mdnSZ_state(bootSZind,:,:);
%         
%         for ii = 1:nS_Kwin_allHC{k2-2}(state)
%             sc_to_sm_corr_HC(ii) = sum(sum(squeeze(mdnHC_state(ii,1:5,6:24))))/numel(squeeze(mdnHC_state(ii,1:5,6:24))); 
%             within_sc_corr_HC(ii) = 0.5*(sum(sum(squeeze(mdnHC_state(ii,1:5,1:5)))))/nchoosek(5,2);
%             within_sm_corr_HC(ii) = 0.5*(sum(sum(squeeze(mdnHC_state(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
%             thal_to_sm_corr_HC(ii) = (sum(squeeze(mdnHC_state(ii,5,6:24))))/numel(6:24); 
%         end
%         
%         for ii = 1:nS_Kwin_allSZ{k2-2}(state)
%             sc_to_sm_corr_SZ(ii) = sum(sum(squeeze(mdnSZ_state(ii,1:5,6:24))))/numel(squeeze(mdnSZ_state(ii,1:5,6:24))); 
%             within_sc_corr_SZ(ii) = 0.5*(sum(sum(squeeze(mdnSZ_state(ii,1:5,1:5)))))/nchoosek(5,2);
%             within_sm_corr_SZ(ii) = 0.5*(sum(sum(squeeze(mdnSZ_state(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
%             thal_to_sm_corr_SZ(ii) = (sum(squeeze(mdnSZ_state(ii,5,6:24))))/numel(6:24); 
%         end
%         
%         r_SC2SMcorr_wSMcorr_HC_boot(bb,state) = corr(sc_to_sm_corr_HC,within_sm_corr_HC);
%         r_THAL2SMcorr_wSMcorr_HC_boot(bb,state) = corr(thal_to_sm_corr_HC,within_sm_corr_HC);
%         r_wSCcorr_wSMcorr_HC_boot(bb,state) = corr(within_sc_corr_HC,within_sm_corr_HC);
%         
%         r_SC2SMcorr_wSMcorr_SZ_boot(bb,state) = corr(sc_to_sm_corr_SZ,within_sm_corr_SZ);
%         r_THAL2SMcorr_wSMcorr_SZ_boot(bb,state) = corr(thal_to_sm_corr_SZ,within_sm_corr_SZ);
%         r_wSCcorr_wSMcorr_SZ_boot(bb,state) = corr(within_sc_corr_SZ,within_sm_corr_SZ);
%     end
    disp(['Done state ' num2str(state)])
end

%save(fullfile(fileparts(inname),  'dfnc_SC_SM_relationships_k5_corr'),  'r_SC*', 'r_THAL*', 'r_wSC*', 'p_SC*', 'p_THAL*', 'p_wSC*', '*allstates')
%%
matched_citystate = [3 1 2 5 4];
matchOrd = [2 3 1 5 4];

F = figure;set(F,'color',[1 1 1],'Name','r_SC2SMcorr_wSMcorr_k5_corr','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    boxplot([r_SC2SMcorr_wSMcorr_HC_boot(:,matchOrd(st));r_SC2SMcorr_wSMcorr_SZ_boot(:,matchOrd(st))],[ones(nBOOT,1);2*ones(nBOOT,1)],'notch','on')
    title(['state ' num2str(st)])
end

F1 = figure;set(F1,'color',[1 1 1],'Name','r_wSCcorr_wSMcorr_k5_corr','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    boxplot([r_wSCcorr_wSMcorr_HC_boot(:,matchOrd(st));r_wSCcorr_wSMcorr_SZ_boot(:,matchOrd(st))],[ones(nBOOT,1);2*ones(nBOOT,1)],'notch','on')
    title(['state ' num2str(st)])
end

F2 = figure;set(F2,'color',[1 1 1],'Name','r_THAL2SMcorr_wSMcorr_k5_corr','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    boxplot([r_THAL2SMcorr_wSMcorr_HC_boot(:,matchOrd(st));r_THAL2SMcorr_wSMcorr_SZ_boot(:,matchOrd(st))],[ones(nBOOT,1);2*ones(nBOOT,1)],'notch','on')
    title(['state ' num2str(st)])
end

F3 = figure;set(F3,'color',[1 1 1],'Name','scatter_SC2SMcorr_wSMcorr_k5_corr','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(sc_to_sm_corr_HC_allstates{matchOrd(st)},within_sm_corr_HC_allstates{matchOrd(st)},'ko');hold;
    plot(sc_to_sm_corr_SZ_allstates{matchOrd(st)},within_sm_corr_SZ_allstates{matchOrd(st)},'ro');
    xlim([-0.4 0.4])
    title(['state ' num2str(st)])
end


F4 = figure;set(F4,'color',[1 1 1],'Name','scatter_wSCcorr_wSMcorr_k5_corr','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(within_sc_corr_HC_allstates{matchOrd(st)},within_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(within_sc_corr_SZ_allstates{matchOrd(st)},within_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    title(['state ' num2str(st)])
end

F5 = figure;set(F5,'color',[1 1 1],'Name','scatter_THAL2SMcorr_wSMcorr_k5_corr','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(thal_to_sm_corr_HC_allstates{matchOrd(st)},within_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(thal_to_sm_corr_SZ_allstates{matchOrd(st)},within_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    title(['state ' num2str(st)])
end

stats_SC2SMcorr_wSMcorr_st3_matchedcorr = mwwtest(r_SC2SMcorr_wSMcorr_HC_boot(:,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,1));
stats_SC2SMcorr_wSMcorr_st1_matchedcorr = mwwtest(r_SC2SMcorr_wSMcorr_HC_boot(:,2),r_SC2SMcorr_wSMcorr_SZ_boot(:,2));
stats_SC2SMcorr_wSMcorr_st2_matchedcorr = mwwtest(r_SC2SMcorr_wSMcorr_HC_boot(:,3),r_SC2SMcorr_wSMcorr_SZ_boot(:,3));
stats_SC2SMcorr_wSMcorr_st5_matchedcorr = mwwtest(r_SC2SMcorr_wSMcorr_HC_boot(:,4),r_SC2SMcorr_wSMcorr_SZ_boot(:,4));
stats_SC2SMcorr_wSMcorr_st4_matchedcorr = mwwtest(r_SC2SMcorr_wSMcorr_HC_boot(:,5),r_SC2SMcorr_wSMcorr_SZ_boot(:,5));

stats_THAL2SMcorr_wSMcorr_st3_matchedcorr = mwwtest(r_THAL2SMcorr_wSMcorr_HC_boot(:,1),r_THAL2SMcorr_wSMcorr_SZ_boot(:,1));
stats_THAL2SMcorr_wSMcorr_st1_matchedcorr = mwwtest(r_THAL2SMcorr_wSMcorr_HC_boot(:,2),r_THAL2SMcorr_wSMcorr_SZ_boot(:,2));
stats_THAL2SMcorr_wSMcorr_st2_matchedcorr = mwwtest(r_THAL2SMcorr_wSMcorr_HC_boot(:,3),r_THAL2SMcorr_wSMcorr_SZ_boot(:,3));
stats_THAL2SMcorr_wSMcorr_st5_matchedcorr = mwwtest(r_THAL2SMcorr_wSMcorr_HC_boot(:,4),r_THAL2SMcorr_wSMcorr_SZ_boot(:,4));
stats_THAL2SMcorr_wSMcorr_st4_matchedcorr = mwwtest(r_THAL2SMcorr_wSMcorr_HC_boot(:,5),r_THAL2SMcorr_wSMcorr_SZ_boot(:,5));

stats_wSCcorr_wSMcorr_st3_matchedcorr = mwwtest(r_wSCcorr_wSMcorr_HC_boot(:,1),r_wSCcorr_wSMcorr_SZ_boot(:,1));
stats_wSCcorr_wSMcorr_st1_matchedcorr = mwwtest(r_wSCcorr_wSMcorr_HC_boot(:,2),r_wSCcorr_wSMcorr_SZ_boot(:,2));
stats_wSCcorr_wSMcorr_st2_matchedcorr = mwwtest(r_wSCcorr_wSMcorr_HC_boot(:,3),r_wSCcorr_wSMcorr_SZ_boot(:,3));
stats_wSCcorr_wSMcorr_st5_matchedcorr = mwwtest(r_wSCcorr_wSMcorr_HC_boot(:,4),r_wSCcorr_wSMcorr_SZ_boot(:,4));
stats_wSCcorr_wSMcorr_st4_matchedcorr = mwwtest(r_wSCcorr_wSMcorr_HC_boot(:,5),r_wSCcorr_wSMcorr_SZ_boot(:,5));



%% mancova
load DEMO_ICA_fBIRNp3
DEM.demo_col_names = DEMO.demo_col_names;% DEMO.demo_col_names(1:4);
DEM.demo_col_var = DEMO.demo_col_var;%DEMO.demo_col_var(1:4);
DEM.demo_labels = DEMO.demo_labels;%DEMO.demo_labels(1:4);
DEM.demo_version = DEMO.demo_version;
vnames = DEM.demo_col_var;
outputdir = '/export/mialab/users/eswar/fbirn_p3/results_C100_vn/FNCdynamics/';

for k2 = 5 % 3:8
    %%
    
    for jj = 1:k2
       statehcid = any(aINDallone_HC == jj,2);hid = find(statehcid == 1);
       stateszid = any(aINDallone_SZ == jj,2);sid = find(stateszid == 1);
       %desXX = DEMO.demo_out([HCid(hid);SZid(sid)],1:4);
       desXX = DEMO.demo_out([HCid(hid);SZid(sid)],:);
       DEM.demo_out = desXX;
       
       MODEL = make_design_matrix_dyn(vnames,[],[],[],DEM);
       
       
       [DEMO1, MULT1, UNI1] = run_model_fbirn(MODEL, 0.01, [mean_Kwin_allHC{k2-2,jj} ;mean_Kwin_allSZ{k2-2,jj}], [], [], 1, 0);
       
       outname = ['MANCOVAN_dynFNCstatsFULL_allone_corr_k' num2str(k2) '_state' num2str(jj) ];
       save([outputdir filesep outname], 'DEMO1', 'MULT1', 'UNI1')
       
       for sigind = 1:length(MULT1.final_terms)
           plot_model_results_SZvsHC(DEMO1,MULT1,UNI1,'FNC',sigind);
           export_fig([FIGoutputdir filesep outname '_' MULT1.final_terms{sigind} '.pdf'],'-pdf',gcf)
           close all
       end
    end
end
%%
k2 = 5;
load(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallone','aINDallone', 'Callone','SUMDallone','Dallone')
IDXalloner = reshape(IDXallone,M,Nwin);
IDXalloner_HC = IDXalloner(HCid,:);
IDXalloner_SZ = IDXalloner(SZid,:);
stc = cell(1,2);len_state = zeros(length(unique(IDXallone)),2);
stall = unique(IDXallone);

sil_val_allonek5_corr = silhouette(Fdynflat,IDXallone,dmethod);
sil_valrCB = reshape(sil_val_allonek5_corr,M,Nwin);
sil_valrCB_SZ = sil_valrCB(SZid,:);
sil_valrCB_HC = sil_valrCB(HCid,:);



for st = 1:length(stall) 
    st_hcid = find(IDXalloner_HC == st);
    st_szid = find(IDXalloner_SZ == st);
    [hst(st) pst(st)] = ttest2(sil_valrCB_HC(st_hcid),sil_valrCB_SZ(st_szid));
    [hkst(st) pkst(st) kst(st)] = kstest2(sil_valrCB_HC(st_hcid),sil_valrCB_SZ(st_szid));
   % len_stateCB(st,:) = [length(st_hcid) length(st_szid)];
    stc{1} = sil_valrCB_HC(st_hcid);
    stc{2} = sil_valrCB_SZ(st_szid);
    figure;distributionPlot(stc);
    clear stc
end
%%  obtain time series spectra for different dyn windows
fbroot = '/export/mialab/users/eswar/fbirn_p3/';
TC = load([ fbroot 'results_C100_vn/TC/fbirnp3_rest_C100_ica_TC_scrubbed']);
TC_icn = TC.TC_scrubbed(keepID,:,RSN_I);
TC_icnz = zscore(TC_icn,[],2);TR = 2;thalID = find(strcmp(L(RSN_I),'Thalamus')==1);
for k2 = 5 %2:9
%     load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) 'scrubbed_cnew']), 'SPflat', 'IDXp','Cp');
%         
%     [IDXallone, Callone, SUMDallone, Dallone] = kmeans(Fdynflat, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
%     
% 
%     aINDallone = reshape(IDXallone, M,Nwin);    
%     
%     save(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_scrubbed_cnew']), 'IDXallone','aINDallone', 'Callone','SUMDallone','Dallone')
%     
    load(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallone','aINDallone', 'Callone','SUMDallone','Dallone');
    aINDallone_SZ = aINDallone(SZid,:);
    aINDallone_HC = aINDallone(HCid,:);
    %%
    for ii = 1:size(aINDallone,1)
        tranSumm.tranWin{ii} = [1 find(abs(diff(aINDallone(ii,:)))>= 1)+1];
        tranSumm.tranSt{ii} = aINDallone(ii,tranSumm.tranWin{ii});
        if ~ismember(tranSumm.tranWin{ii},size(aINDallone,2))
            tranSumm.tranLen{ii} = diff([1 find(abs(diff(aINDallone(ii,:)))>= 1)+1 size(aINDallone,2) ]);
        else
            disp(ii)
            tranSumm.tranLen{ii} = [diff([1 find(abs(diff(aINDallone(ii,:)))>= 1) ]) 1];
        end
    end
    %%
    longStateList = cell(size(aINDallone,1),1);cnt1 = 0;exst = zeros(1,length(keepID));    
    exst1 = zeros(k2,length(keepID));
    
    for ii = 1:size(aINDallone,1)
        longStid = find(tranSumm.tranLen{ii} > 30);
        if ~isempty(longStid)
            longStateList{ii} = tranSumm.tranSt{ii}(longStid);cnt1 = cnt1+1;exst(ii) = 1;
            
            for kk = longStateList{ii},                
                exst1(kk,ii) = exst1(kk,ii)+1;       
                
            end
            disp([ii  tranSumm.tranSt{ii}(longStid)  tranSumm.tranLen{ii}(longStid)])
        end
    end
    %%
    subj_with_longSt = find(sum(exst1,1) > 0);cnt_longSt_bySubj = sum(exst1,1);
    cnt_longWin_bySt =  sum(exst1,2);
    state1cnt = 0;state2cnt = 0;state3cnt =0; state4cnt = 0;state5cnt = 0;
    
    for ii = 1:length(subj_with_longSt)
        
        cii = subj_with_longSt(ii);
        ckid = find(exst1(:,cii) > 0);
        st = find(tranSumm.tranLen{cii} > 30);ck1 = 0;
        stlist = find(exst1(:,cii) > 0);
        for ck = 1:length(stlist)            
            nreps = exst1(ckid(ck),cii);
            mystate = stlist(ck);
            stlid = find(longStateList{cii} == mystate);
            
            for ckk = 1:nreps
                ck1 = ck1+1;
                
                disp([cii  st(ckk) mystate ckk]) 
                
                cwin = tranSumm.tranWin{cii}(st(stlid(ckk)))+1:tranSumm.tranWin{cii}(st(stlid(ckk)))+(tranSumm.tranLen{cii}(st(stlid(ckk))))-2;
                sigma = 3; % for smoothing Gaussian
                wsize = length(cwin); % for box (44 seconds)
                nT = size(TC_icnz,2);
                gw = gaussianwindow(nT,nT/2,sigma);
                b = zeros(nT,1);  
                if mod(wsize,2) ==1                    
                    wsize = wsize -1;
                end
                
                b((nT/2 - wsize/2 + 1):(nT/2+wsize/2)) = 1;
                c = conv(gw,b); c = c/max(c); c=c(nT/2+1:end-nT/2 +1);
                cshift = circshift(c,-nT/2 + wsize/2 + cwin(1));
                
                [S f] = compute_spectrum(squeeze(TC_icnz(cii,:,thalID)).*cshift',TR,1);
%                 if length(f) > 33
%                     S = S(1:2:end);
%                 end
                SS(1:length(S),ckk) = S;
                
                
            end
                
            clear ckk 
            
            if size(SS,2) > 1
                SS = mean(SS,2);
            end
            
            
            switch mystate
                case 1
                    state1cnt = state1cnt+1;                    
                    spec_state1(state1cnt,1:length(SS)) = SS;
                    
                case 2
                    state2cnt = state2cnt+1;
                    spec_state2(state2cnt,1:length(SS)) = SS;
                case 3
                    state3cnt = state3cnt+1;
                    spec_state3(state3cnt,1:length(SS)) = SS;
                case 4
                    state4cnt = state4cnt+1;
                    spec_state4(state4cnt,1:length(SS)) = SS;
                case 5
                    state5cnt = state5cnt+1;
                    spec_state5(state5cnt,1:length(SS)) = SS;
            end
                    
                    
             clear SS;       
             %   Sf = 
            
        end
    end
                
end

save spectra_thal_by_state_kmeans_corr_k5_thr30  tranSumm longStateList spec_state*
%% compare spectra by diagnosis by state
% F = figure; hold on
% set(F,'Name','ThalamusSpecra_byState')
matched_citystate = [3 1 2 5 4];
for state = 1:5
    [tf1 loc1] = ismember(find(exst1(matched_citystate(state),:) > 0),HCid);
    [tf2 loc2] = ismember(find(exst1(matched_citystate(state),:) > 0),SZid);
    disp([ matched_citystate(state) length(find(tf1 == 1)) length(find(tf2 == 1))])
    eval(['hcdata =  spec_state' num2str(matched_citystate(state)) '(tf1==1,:);']);
    eval(['szdata =  spec_state' num2str(matched_citystate(state)) '(tf2==1,:);']);
    if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8) 
        [hhh ppp] = ttest2(hcdata,szdata);
        figure;plot(sortrows([ppp' myFDR(ppp)']))
        title(['state' num2str(matched_citystate(state))])
    end
%      Fs = subplot(3,2,state);
%      plot_with_ste_area(Fs,f,hcdata);hold on
%      plot_with_ste_area(Fs,f,szdata,[],'b','c');
%      title(['state ' num2str((state))])
%      CH = get(Fs,'children');
%      legend(CH([2 5]),{'SZ' 'HC'});
end


%% spectra from all state segments


astate1cnt = 0;astate2cnt = 0;astate3cnt =0; astate4cnt = 0;astate5cnt = 0;

gw = gaussianwindow(nT,nT/2,sigma);
sigma = 3; % for smoothing Gaussian
nT = size(TC_icnz,2);

for ii = 1:length(tranSumm.tranWin)
    
    cii = ii;
    
    st = tranSumm.tranLen{cii};
    
    stlist = unique(tranSumm.tranSt{cii});
    for ck = 1:length(stlist)
        nreps = length(find(tranSumm.tranSt{cii} == stlist(ck)));
        mystate = stlist(ck);
        stlid = find(tranSumm.tranSt{cii} == mystate);
        
        for ckk = 1:nreps
            
            
            disp([cii  st(ckk) mystate ckk])
            
            cwin = tranSumm.tranWin{cii}(stlid(ckk))+1:tranSumm.tranWin{cii}(stlid(ckk))+(tranSumm.tranLen{cii}(stlid(ckk)));
            
            if length(cwin) < 22 %
            	cwin = cwin(1):cwin(1)+22;
            end                 
            
            wsize = length(cwin); % for box (44 seconds)
            
            
            b = zeros(nT,1);
            if mod(wsize,2) ==1
                wsize = wsize -1;
            end
            
            b((nT/2 - wsize/2 + 1):(nT/2+wsize/2)) = 1;
            c = conv(gw,b); c = c/max(c); c=c(nT/2+1:end-nT/2 +1);
            cshift = circshift(c,-nT/2 + wsize/2 + cwin(1));
            
            [S f] = compute_spectrum(squeeze(TC_icnz(cii,:,thalID)).*cshift',TR,1);
            %                 if length(f) > 33
            %                     S = S(1:2:end);
            %                 end
            SS(1:length(S),ckk) = S;
            
        end
        
        clear ckk
        
        if size(SS,2) > 1
            SS = mean(SS,2);
        end
        
        
        switch mystate
            case 1
                astate1cnt = astate1cnt+1;
                aspec_state1(astate1cnt,1:length(SS)) = SS;
                astate1_subjid(astate1cnt) = cii;
            case 2
                astate2cnt = astate2cnt+1;
                aspec_state2(astate2cnt,1:length(SS)) = SS;
                astate2_subjid(astate2cnt) = cii;
            case 3
                astate3cnt = astate3cnt+1;
                aspec_state3(astate3cnt,1:length(SS)) = SS;
                astate3_subjid(astate3cnt) = cii;
            case 4
                astate4cnt = astate4cnt+1;
                aspec_state4(astate4cnt,1:length(SS)) = SS;
                astate4_subjid(astate4cnt) = cii;
            case 5
                astate5cnt = astate5cnt+1;
                aspec_state5(astate5cnt,1:length(SS)) = SS;
                astate5_subjid(astate5cnt) = cii;
        end
        
        
        clear SS;
        %   Sf =
        
    end
end
%%
 F = figure; hold on
 set(F,'Name','ThalamusSpecra_byState_corr')
%load spectra_thal_by_state_kmeans_city_k5_thr30
matched_citystate = [3 1 2 5 4];
for state = 1:5
    state1 = matched_citystate(state);
    eval(['[tf1 loc1] = ismember(astate' num2str(state) '_subjid,HCid);']);
    eval(['[tf2 loc2] = ismember(astate' num2str(state) '_subjid,SZid);']);
    
    disp([ state length(find(tf1 == 1)) length(find(tf2 == 1))])
    eval(['hcdata =  aspec_state' num2str(state) '(tf1==1,:);']);
    eval(['szdata =  aspec_state' num2str(state) '(tf2==1,:);']);
     if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8) 
         [hhh(state,:) ppp(state,:)] = ttest2(hcdata,szdata);
         %figure;plot(sortrows([ppp' myFDR(ppp)']))
         %title(['state ' num2str(state)])
         
     end
     Fs = subplot(3,2,state1);
     plot_with_ste_area(Fs,f,hcdata);hold on
     plot_with_ste_area(Fs,f,szdata,[],'b','c');
     title(['state ' num2str(state1)]);
     CH = get(Fs,'children');
     legend(CH([2 5]),{'SZ' 'HC'});
end

%%
for state = 1:5
    state1 = matched_citystate(state);
    figure;plot(sortrows([ppp(state1,:)' myFDR(ppp(state1,:))']))
    title(['state ' num2str(state1)])
end
%%
%%  spectra for all segments

gw = gaussianwindow(nT,nT/2,sigma);
sigma = 3; % for smoothing Gaussian
nT = size(TC_icnz,2);
for rsn =  [1 5 7 16  23 24 29 47]
    astate1cnt = 0;astate2cnt = 0;astate3cnt =0; astate4cnt = 0;astate5cnt = 0;


    for ii = 1:length(tranSumm.tranWin)
        
        cii = ii;
        
        st = tranSumm.tranLen{cii};
        
        stlist = unique(tranSumm.tranSt{cii});
        for ck = 1:length(stlist)
            nreps = length(find(tranSumm.tranSt{cii} == stlist(ck)));
            mystate = stlist(ck);
            stlid = find(tranSumm.tranSt{cii} == mystate);
            
            for ckk = 1:nreps
                
                
                disp([cii  st(ckk) mystate ckk])
                
                cwin = tranSumm.tranWin{cii}(stlid(ckk))+1:tranSumm.tranWin{cii}(stlid(ckk))+(tranSumm.tranLen{cii}(stlid(ckk)));
                
                if length(cwin) < 22 %
                    cwin = cwin(1):cwin(1)+22;
                end
                
                wsize = length(cwin); % for box (44 seconds)
                
                
                b = zeros(nT,1);
                if mod(wsize,2) ==1
                    wsize = wsize -1;
                end
                
                b((nT/2 - wsize/2 + 1):(nT/2+wsize/2)) = 1;
                c = conv(gw,b); c = c/max(c); c=c(nT/2+1:end-nT/2 +1);
                cshift = circshift(c,-nT/2 + wsize/2 + cwin(1));
                
                [S f] = compute_spectrum(squeeze(TC_icnz(cii,:,rsn)).*cshift',TR,1);
                %                 if length(f) > 33
                %                     S = S(1:2:end);
                %                 end
                SS(1:length(S),ckk) = S;
                
            end
            
            clear ckk
            
            if size(SS,2) > 1
                SS = mean(SS,2);
            end
            
            
            switch mystate
                case 1
                    astate1cnt = astate1cnt+1;
                    aspec_state1(astate1cnt,1:length(SS)) = SS;
                    astate1_subjid(astate1cnt) = cii;
                case 2
                    astate2cnt = astate2cnt+1;
                    aspec_state2(astate2cnt,1:length(SS)) = SS;
                    astate2_subjid(astate2cnt) = cii;
                case 3
                    astate3cnt = astate3cnt+1;
                    aspec_state3(astate3cnt,1:length(SS)) = SS;
                    astate3_subjid(astate3cnt) = cii;
                case 4
                    astate4cnt = astate4cnt+1;
                    aspec_state4(astate4cnt,1:length(SS)) = SS;
                    astate4_subjid(astate4cnt) = cii;
                case 5
                    astate5cnt = astate5cnt+1;
                    aspec_state5(astate5cnt,1:length(SS)) = SS;
                    astate5_subjid(astate5cnt) = cii;
            end
            
            
            clear SS;
            %   Sf =
            
        end
    end
    outputdir = '/export/mialab/users/eswar/fbirn_p3/results_C100_vn/FNC/';
    outn = [fileparts(outputdir) 'dynamics' filesep 'spectra_k5_corr_by_state_ICN_' L{RSN_I(rsn)}];
    save(outn,'astate*','aspec*')

    CPOS2 = [ 180   448   820   612];
    Fa = figure;set(Fa,'Color',[1 1 1], 'Position', CPOS2); hold on
    set(Fa,'Name',[L{RSN_I(rsn)}  'Specra_byState_corr'])
    %load spectra_thal_by_state_kmeans_city_k5_thr30
    for state = 1:5
        state1 = matched_citystate(state);
        eval(['[tf1 loc1] = ismember(astate' num2str(state) '_subjid,HCid);']);
        eval(['[tf2 loc2] = ismember(astate' num2str(state) '_subjid,SZid);']);
        
        disp([ state length(find(tf1 == 1)) length(find(tf2 == 1))])
        eval(['hcdata =  aspec_state' num2str(state) '(tf1==1,:);']);
        eval(['szdata =  aspec_state' num2str(state) '(tf2==1,:);']);
        if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8)
            [hhh(state,:) ppp(state,:)] = ttest2(hcdata,szdata);
            %figure;plot(sortrows([ppp' myFDR(ppp)']))
            %title(['state ' num2str(state)])
            
        end
        Fs = subplot(3,2,state1);
        plot_with_ste_area(Fs,f,hcdata);hold on
        plot_with_ste_area(Fs,f,szdata,[],'b','c');
        title(['state ' num2str(state1)]);
        set(Fs,'FontSize',14)
        ylabel('amplitude A.U.')
        xlabel('frequency')
        ylim([0.2 1.2])
        CH = get(Fs,'children');
        legend(CH([2 5]),{'SZ' 'HC'});
        
    end
    export_fig(fullfile(FIGoutputdir, [get(Fa, 'Name') '_allwin.pdf']), Fa);

    CPOS3 = [  180         598        1056         462];
    F = figure; set(F,'Color',[1 1 1], 'Position', CPOS3);
    set(F,'Name',[L{RSN_I(rsn)}  'Specra_byStateGroup_corr'])
    %load spectra_thal_by_state_kmeans_city_k5_thr30
    clrf = [0 0 0;0 0 1;1 0 0; 0.6 1 0.2; 0.3333    0.0196    0.3686];
    clrb = [.5 .5 .5; .5373 .8118 .9412;1 1 0;0 1 0;.929  .2235  .98];
    clrf = clrf([2 3 1 5 4],:);clrb = clrb([2 3 1 5 4],:);
    for state1 = 1:5
        state = matched_citystate(state1);
        eval(['[tf1 loc1] = ismember(astate' num2str(state) '_subjid,HCid);']);
        eval(['[tf2 loc2] = ismember(astate' num2str(state) '_subjid,SZid);']);

        disp([ state length(find(tf1 == 1)) length(find(tf2 == 1))])
        eval(['hcdata =  aspec_state' num2str(state) '(tf1==1,:);']);
        eval(['szdata =  aspec_state' num2str(state) '(tf2==1,:);']);
        if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8)
            [hhh(state,:) ppp(state,:)] = ttest2(hcdata,szdata);
            %figure;plot(sortrows([ppp' myFDR(ppp)']))
            %title(['state ' num2str(state)])

        end
        Fs1 = subplot(1,2,1);hold on
        plot_with_ste_area(Fs1,f,hcdata,[],clrf(state,:),clrb(state,:));
        ylim([0.2 1.2])
        if state1 == 5            
            title(['HC: ' L{RSN_I(rsn)} ' Mean Power Spectra']);
            CH1 = get(Fs1,'children');
            legend(CH1([2:3:14]),{'state5' 'state4' 'state3' 'state2' 'state1'});
            set(Fs1,'FontSize',14)
            ylabel('amplitude A.U.')
            xlabel('frequency')
        end

        Fs2 = subplot(1,2,2); hold on
        plot_with_ste_area(Fs2,f,szdata,[],clrf(state,:),clrb(state,:));
        ylim([0.2 1.2])
        if state1 == 5
            title(['SZ: ' L{RSN_I(rsn)} ' Mean Power Spectra']);
            CH2 = get(Fs2,'children');
            legend(CH2([2:3:14]),{'state5' 'state4' 'state3' 'state2' 'state1'});
            set(Fs2,'FontSize',14)
            ylabel('amplitude A.U.')
            xlabel('frequency')
        end

        % CH = get(Fs,'children');
        % legend(CH([2 5]),{'SZ' 'HC'});
    end
    export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '_allwin.pdf']), F);
    clear astate* aspec*
end

%%  spectra for all dynamic FNC segments
k2 = 5;
load(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'IDXallone','aINDallone', 'Callone','SUMDallone','Dallone');
aINDallone_SZ = aINDallone(SZid,:);
aINDallone_HC = aINDallone(HCid,:);
load('/export/mialab/users/eswar/fbirn_p3/results_C100_vn/FNCdynamics/Spec_dynamicSWTukey22_nowrap_notaper'); % loads dynSpec_nt_nw_HC [subj * nWin * nC* freq]
%load('/export/mialab/users/eswar/fbirn_p3/results_C100_vn/FNCdynamics/Spec_dynamicSWTukey22_nowrap_taper'); % loads dynSpec_t_nw_HC [subj * nWin * nC* freq]

dynSpec_nt_nw_HC = dynSpec_nt_nw(HCid,:,:,:,:);
dynSpec_nt_nw_SZ = dynSpec_nt_nw(SZid,:,:,:,:);
%dynSpec_nt_nw_HC = reshape(dynSpec_nt_nw_HC,[numel(aINDallone_HC) length(RSN_I) length(f)]);
%dynSpec_nt_nw_SZ = reshape(dynSpec_nt_nw_SZ,[numel(aINDallone_SZ) length(RSN_I) length(f)]);

%dynSpec_t_nw_HC = dynSpec_t_nw(HCid,:,:,:,:);
%dynSpec_t_nw_SZ = dynSpec_t_nw(SZid,:,:,:,:);
%dynSpec_t_nw_HC = reshape(dynSpec_t_nw_HC,[numel(aINDallone_HC) length(RSN_I) length(f)]);
%dynSpec_t_nw_SZ = reshape(dynSpec_t_nw_SZ,[numel(aINDallone_SZ) length(RSN_I) length(f)]);


CPOS2 = [ 180   448   820   612];
CPOS3 = [  180         598        1056         462];
clrf = [0 0 0;0 0 1;1 0 0; 0.6 1 0.2; 0.3333    0.0196    0.3686];
clrb = [.5 .5 .5; .5373 .8118 .9412;1 1 0;0 1 0;.929  .2235  .98];
f = 0:0.25/128:0.25;
lf = find(f> 0.024 & f < 0.08);
dm = load('DEMO_ICA_fBIRNp3');
desX = dm.DEMO.demo_out;
labs = dm.DEMO.demo_col_names;

desX(:,strcmp(labs,'age')) = bsxfun(@minus, desX(:,strcmp(labs,'age')) , mean(desX(:,strcmp(labs,'age'))));
desX(:,strcmp(labs,'meanFD')) = bsxfun(@minus, desX(:,strcmp(labs,'meanFD')) , mean(desX(:,strcmp(labs,'meanFD'))));
desX_HC = desX(HCid,[1 2 4:end]);desX_SZ = desX(SZid,[1 2 4:end]);
load('DEMO_ICA_fBIRNp3')
ICATC = load('../results_C100_vn/TC/ICATC_metrics.mat');
fd_dv_orig = load('FD_DVARS_fbirn_p3_ICA');
goodindHC = (intersect(find(DEMO.demo_out(:,4) < 0.2),find(DEMO.demo_out(:,3) ==0) ));
goodindSZ = (intersect(find(DEMO.demo_out(:,4) < 0.2),find(DEMO.demo_out(:,3) ==1) ));

desX_HCltd = desX(HCid,[1 2 4]);
desX_SZltd = desX(SZid,[1 2 4]);
%%
 %FdynZ  = zscore(permute(Fdyn,[3 1 2]));
 %FdynZ = permute(FdynZ,[  2 3 1]);
 FdynflatHC = Fdyn(HCid,:,:);
 FdynflatSZ = Fdyn(SZid,:,:);


for rsn = 4%[1  5 7 16  23 24 29 30 47]
   
    
%     Fa = figure;set(Fa,'Color',[1 1 1], 'Position', CPOS2); hold on
%     set(Fa,'Name',[L{RSN_I(rsn)}  'Specra_byState_city'])
%     
%     F = figure; set(F,'Color',[1 1 1], 'Position', CPOS3);
%     set(F,'Name',[L{RSN_I(rsn)}  'Specra_byStateGroup_city'])
   %%
    t_obs_hc = zeros(1,5);p_obs_hc = t_obs_hc;
    p_obs_sz = t_obs_hc;t_obs_sz = t_obs_hc;
    b_obs_hc = t_obs_hc;b_obs_sz = t_obs_hc;
    t_obs_hc_goodind = zeros(1,5);p_obs_hc_goodind = t_obs_hc;
    p_obs_sz_goodind = t_obs_hc;t_obs_sz_goodind = t_obs_hc;
    b_obs_hc_goodind = t_obs_hc;b_obs_sz_goodind = t_obs_hc;
  %  b_boot_state_hc = zeros(k2,1000);b_boot_state_sz = zeros(k2,1000);
    
    for state = 1:k2     
        
            
        statehcid = any(aINDallone_HC == state,2);hid = find(statehcid == 1);
        dynSpec_nt_nw_HC_st = zeros(sum(statehcid),numel(f));normHC = dynSpec_nt_nw_HC_st;
        for ss = 1:sum(statehcid)
            dynSpec_nt_nw_HC_st(ss,:) = squeeze(mean(dynSpec_nt_nw_HC(hid(ss),aINDallone_HC(hid(ss),:)==state,rsn,:),2));
            normHC(ss,:) = dynSpec_nt_nw_HC_st(ss,:)/max(dynSpec_nt_nw_HC_st(ss,:));
        end
        stateszid = any(aINDallone_SZ == state,2);sid = find(stateszid == 1);
        dynSpec_nt_nw_SZ_st = zeros(sum(stateszid),numel(f));normSZ = dynSpec_nt_nw_SZ_st;
        for ss = 1:sum(stateszid)
            dynSpec_nt_nw_SZ_st(ss,:) = squeeze(mean(dynSpec_nt_nw_SZ(sid(ss),aINDallone_SZ(sid(ss),:)==state,rsn,:),2));
            normSZ(ss,:) = dynSpec_nt_nw_SZ_st(ss,:)/max(dynSpec_nt_nw_SZ_st(ss,:));
        end                  
            
        
        hcdata = normHC;   %dynSpec_t_nw_HC_st;
        szdata = normSZ;  %dynSpec_t_nw_SZ_st;        
        
        %%%%%%%%%
        %if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8)
        %    [hhsdyn(state,:) ppsdyn(state,:)] = ttest2(hcdata,szdata);
            %figure;plot(sortrows([ppp' myFDR(ppp)']))
            %title(['state ' num2str(state)])
            
        %end
        %%%%%%%%
        
%         figure(Fa);
%         Fs = subplot(3,2,state);
%         plot_with_ste_area(Fs,f,hcdata);hold on
%         plot_with_ste_area(Fs,f,szdata,[],'b','c');
%         title(['state ' num2str(state)]);
%         set(Fs,'FontSize',14)
%         ylabel('amplitude A.U.')
%         xlabel('frequency')
%         ylim([0 1.0])
%         CH = get(Fs,'children');
%         legend(CH([2 5]),{'SZ' 'HC'});
%         
%         
%         figure(F);
%         Fs1 = subplot(1,2,1);hold on
%         plot_with_ste_area(Fs1,f,hcdata,[],clrf(state,:),clrb(state,:));
%         ylim([0  1.0])
%         if state == 5            
%             title(['HC: ' L{RSN_I(rsn)} ' Mean Power Spectra']);
%             CH1 = get(Fs1,'children');
%             legend(CH1([2:3:14]),{'state5' 'state4' 'state3' 'state2' 'state1'});
%             set(Fs1,'FontSize',14)
%             ylabel('amplitude A.U.')
%             xlabel('frequency')
%         end
% 
%         Fs2 = subplot(1,2,2); hold on
%         plot_with_ste_area(Fs2,f,szdata,[],clrf(state,:),clrb(state,:));
%         ylim([0 1.0])
%         if state == 5
%             title(['SZ: ' L{RSN_I(rsn)} ' Mean Power Spectra']);
%             CH2 = get(Fs2,'children');
%             legend(CH2([2:3:14]),{'state5' 'state4' 'state3' 'state2' 'state1'});
%             set(Fs2,'FontSize',14)
%             ylabel('amplitude A.U.')
%             xlabel('frequency')
%         end
%     end
%     
%     export_fig(fullfile(FIGoutputdir, [get(Fa, 'Name') '_allwin_norm.pdf']), Fa);
%     close(Fa) 
%     
%       
%     
%     export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '_allwin_norm.pdf']), F);
%     close(F)
%     
    
    
    tempHC = FdynflatHC(hid,:,:);    
    tempmat1 = zeros([ length(hid)         136   length(RSN_I)  length(RSN_I) ]);
    for ss = 1:length(hid), tempmat1(ss,:,:,:) = vec2mat(squeeze(tempHC(ss,:,:)));end
    tempmat_thalrhc = squeeze(mean(tempmat1(:,:,rsn,6:24),4));
    tempmat_SMrhc = squeeze(sum(sum(tempmat1(:,:,6:24,6:24),4),3))/342;
    
    sindhc = aINDallone_HC(hid,:) == state;tempmat_thalrhc_mean = zeros(length(hid),1);
    tempmat_SMrhc_mean= zeros(length(hid),1);
    for ss = 1:length(hid),
        tempmat_thalrhc_mean(ss,:) = squeeze(mean(tempmat_thalrhc(ss,sindhc(ss,:)),2));
        tempmat_SMrhc_mean(ss,:) = squeeze(mean(tempmat_SMrhc(ss,sindhc(ss,:)),2));
    end
      
    %j = sum(normHC(:,lf),2)./sum(normHC(:,hf),2);
    
    tempSZ = FdynflatSZ(sid,:,:);
    tempmat2 = zeros([ length(sid)         136   length(RSN_I)  length(RSN_I) ]);
    for ss = 1:length(sid), tempmat2(ss,:,:,:) = vec2mat(squeeze(tempSZ(ss,:,:)));end
    tempmat_thalrsz = squeeze(mean(tempmat2(:,:,rsn,6:24),4));
    tempmat_SMrsz = squeeze(sum(sum(tempmat2(:,:,6:24,6:24),4),3))/342;
    
    sindsz = aINDallone_SZ(sid,:) == state;tempmat_thalrsz_mean = zeros(length(sid),1);
    tempmat_SMrsz_mean= zeros(length(sid),1);
    for ss = 1:length(sid),
        tempmat_thalrsz_mean(ss,:) = squeeze(mean(tempmat_thalrsz(ss,sindsz(ss,:)),2));
        tempmat_SMrsz_mean(ss,:) = squeeze(mean(tempmat_SMrsz(ss,sindsz(ss,:)),2));
    end
%j1 = sum(normSZ(:,lf),2)./sum(normSZ(:,hf),2);
    figure; plot([tempmat_thalrhc_mean;],[tempmat_SMrhc_mean;],'ko');
    hold on;scatter([tempmat_thalrsz_mean],[tempmat_SMrsz_mean],'ro');
    title(['State : ' num2str(state)])
    
    spectHC_lf_sum = sum(dynSpec_nt_nw_HC_st(:,lf),2);    
    spectSZ_lf_sum = sum(dynSpec_nt_nw_SZ_st(:,lf),2);
    
%      Fc = figure;set(Fc,'Color',[1 1 1],'name',[L{RSN_I(rsn)} 'LFspect_' L{RSN_I(rsn)} 'Corr_state'  num2str(state) ], 'Position', [ 7   351   854   658]); 
%      plot(spectHC_lf_sum, tempmat_thalrhc_mean, 'b.'); hold on; lsline; plot(spectSZ_lf_sum, tempmat_thalrsz_mean, 'r.'); lsline
%      xlabel('LF power');ylabel([ L{RSN_I(rsn)} ' correlation with sensory ICNs'])
%      set(gca,'FontSize',14)
%      export_fig(fullfile(FIGoutputdir, [get(Fc, 'Name') '_allwin_norm.pdf']), Fc);
%      close(Fc)
    
     zind = find(sum(desX_HCltd(hid,:)) == 0);
     okind1 = setdiff(1:size(desX_HCltd,2),zind); 
     rstathc = regstats(tempmat_thalrhc_mean,[desX_HCltd(hid,okind1) spectHC_lf_sum]); 
     zind = find(sum(desX_SZltd(sid,:)) == 0);
     okind2 = setdiff(1:size(desX_SZltd,2),zind); 
     rstatsz = regstats(tempmat_thalrsz_mean,[desX_SZltd(sid,okind2) spectSZ_lf_sum]);
     
     
     
     %remhc = (desX_HC(hid,okind1))* rstathc.tstat.beta(2:end-1);
     %remsz = (desX_SZ(sid,okind2))* rstatsz.tstat.beta(2:end-1);

%     zind = find(sum(desX_HC(hid,:)) == 0);
%     okind1 = setdiff(1:size(desX_HC,2),zind); 
%     rstathc = regstats(tempmat_thalrhc_mean,[desX_HC(hid,okind1) spectHC_lf_sum]); 
%     zind = find(sum(desX_SZ(sid,:)) == 0);
%     okind2 = setdiff(1:size(desX_SZ,2),zind); 
%     rstatsz = regstats(tempmat_thalrsz_mean,[desX_SZ(sid,okind2) spectSZ_lf_sum]);
%     
%     remhc = (desX_HC(hid,okind1))* rstathc.tstat.beta(2:end-1);
%     remsz = (desX_SZ(sid,okind2))* rstatsz.tstat.beta(2:end-1);
%     tempmat_thalrhc_meanortho = tempmat_thalrhc_mean - remhc;
%     tempmat_thalrsz_meanortho = tempmat_thalrsz_mean - remsz;
    [intv1 inti1] = intersect(HCid(hid),goodindHC);
    [intv2 inti2] = intersect(SZid(sid),goodindSZ);

    rstathc_goodind = regstats(tempmat_thalrhc_mean(inti1),[desX_HCltd(hid(inti1),okind1) spectHC_lf_sum(inti1)]);      
    rstatsz_goodind  = regstats(tempmat_thalrsz_mean(inti2),[desX_SZltd(sid(inti2),okind2) spectSZ_lf_sum(inti2)]);
    rstatdiff = regstats([tempmat_thalrhc_mean;tempmat_thalrsz_mean],[desX_HCltd(hid,okind1) zeros(length(hid),1) spectHC_lf_sum zeros(length(hid),1);desX_SZltd(sid,okind2) ones(length(sid),1) spectSZ_lf_sum  spectSZ_lf_sum]); % diag * lfpower interaction
    rstatdiff_goodind = regstats([tempmat_thalrhc_mean(inti1);tempmat_thalrsz_mean(inti2)],[desX_HCltd(hid(inti1),okind1) zeros(length(inti1),1) spectHC_lf_sum(inti1) zeros(length(inti1),1) ;desX_SZltd(sid(inti2),okind2) ones(length(inti2),1) spectSZ_lf_sum(inti2)  spectSZ_lf_sum(inti2)]); 

    
    FF = figure;set(FF,'Name',['corr_'  L{RSN_I(rsn)} 'toSensoryICN_vs_' L{RSN_I(rsn)} 'LFamp_state' num2str(state) '_corrZ'],'Color',[1 1 1]);
    scatter(spectHC_lf_sum,tempmat_thalrhc_mean,'k');hold on;
    HH1 = lsline;  set(HH1,'Color','k','LineWidth',1)
    
    scatter(spectSZ_lf_sum,tempmat_thalrsz_mean,'r');
    HH2 = lsline;set(HH2,'Color','r','LineWidth',1)
    scatter(spectHC_lf_sum(inti1),tempmat_thalrhc_mean(inti1),'c','filled');
    scatter(spectSZ_lf_sum(inti2),tempmat_thalrsz_mean(inti2),'m','filled');
    
    export_fig(fullfile(FIGoutputdir,[get(FF,'Name') '.pdf']),'-pdf',FF)
%     
     t_obs_hc(state) = rstathc.tstat.t(end);p_obs_hc(state) = rstathc.tstat.pval(end);
     t_obs_sz(state) = rstatsz.tstat.t(end);p_obs_sz(state) = rstatsz.tstat.pval(end);
     b_obs_hc(state) = rstathc.tstat.beta(end);b_obs_sz(state) = rstatsz.tstat.beta(end);
     t_obs_hc_goodind(state) = rstathc_goodind.tstat.t(end);p_obs_hc_goodind(state) = rstathc_goodind.tstat.pval(end);
     t_obs_sz_goodind(state) = rstatsz_goodind.tstat.t(end);p_obs_sz_goodind(state) = rstatsz_goodind.tstat.pval(end);
     b_obs_hc_goodind(state) = rstathc_goodind.tstat.beta(end);b_obs_sz_goodind(state) = rstatsz_goodind.tstat.beta(end);
     t_obs_diff(state) = rstatdiff.tstat.t(end);p_obs_diff(state) = rstatdiff.tstat.pval(end);b_obs_diff(state) = rstatdiff.tstat.beta(end);
     t_obs_diff_goodind(state) = rstatdiff_goodind.tstat.t(end);p_obs_diff_goodind(state) = rstatdiff_goodind.tstat.pval(end);b_obs_diff_goodind(state) = rstatdiff_goodind.tstat.beta(end);

%     
%     for bb = 1:1000
%         pickME1 = ceil(rand(1,length(hid))*length(hid));
%         zind = find(sum(desX_HC(hid(pickME1),:)) == 0);
%         okind = setdiff(1:size(desX_HC,2),zind);        
%         brstathc = regstats(tempmat_thalrhc_mean(pickME1),[desX_HC(hid(pickME1),okind) spectHC_lf_sum(pickME1)]);
%         b_boot_state_hc(state,bb) = brstathc.tstat.beta(end);
%         
%         pickME2 = ceil(rand(1,length(sid))*length(sid));
%         zind = find(sum(desX_SZ(sid(pickME2),:)) == 0);
%         okind = setdiff(1:size(desX_SZ,2),zind);
%         
%         brstatsz = regstats(tempmat_thalrsz_mean(pickME2),[desX_SZ(sid(pickME2),okind) spectSZ_lf_sum(pickME2)]);
%         b_boot_state_sz(state,bb) = brstatsz.tstat.beta(end);
%     end

    clear hcdata szdata
    end
    save(['corr_' L{RSN_I(rsn)} 'toSensoryICN_vs_' L{RSN_I(rsn)} 'LFamp_mFD_Diag_allone_corr'],    'p_obs_*', 't_obs_*', 'b_obs_*');
end
%save corr_thaltoSensoryICN_vs_thalLFamp  b_boot_state_*  p_obs_* t_obs_* b_obs_*

%% 
matchOrd = [2 3 1 5 4];
for rsn = [4 5]%[1  5 7 16  23 24 29 30 47]
    %rsn = 30;
    cRSN = load(['corr_' L{RSN_I(rsn)} 'toSensoryICN_vs_' L{RSN_I(rsn)} 'LFamp_mFD_Diag_allone_corr']);
    disp(['---------------------- ' num2str(rsn) ' : ' L{RSN_I(rsn)}   ' -------------------------'])
    [cRSN.t_obs_hc(matchOrd);cRSN.t_obs_sz(matchOrd);cRSN.t_obs_diff(matchOrd)]
    disp(['-----------------------------------------------------------'])
    [cRSN.p_obs_hc(matchOrd);cRSN.p_obs_sz(matchOrd);cRSN.p_obs_diff(matchOrd)]
    disp(['-----------------------------------------------------------'])
    [cRSN.b_obs_hc(matchOrd);cRSN.b_obs_sz(matchOrd);cRSN.b_obs_diff(matchOrd)]
    disp(['%%%%%%%---------------------------------------------------------%%%%%%%'])
    
%     F = figure;set(F,'Color',[1 1 1],'Name',['beta_spect'  L{RSN_I(rsn)}  '_vs_r'  L{RSN_I(rsn)} '_allone_corr' ]);
%     errorbar(1:5,cRSN.b_obs_hc(matchOrd),cRSN.b_obs_hc(matchOrd)./cRSN.t_obs_hc(matchOrd),'k');xlim([0.5 5.5]); set(gca,'XTick',[1:5])
%     hold on;errorbar(1:5,cRSN.b_obs_sz(matchOrd),cRSN.b_obs_sz(matchOrd)./cRSN.t_obs_sz(matchOrd),'r');xlim([0.5 5.5]); set(gca,'XTick',[1:5])
%     plot(1:5,[0 0 0 0 0],'k--')
%     export_fig(fullfile(FIGoutputdir,[get(F,'Name') '.pdf']),'-pdf',F)
%     
%     F1 = figure;set(F1,'Color',[1 1 1],'Name',['beta_spect'  L{RSN_I(rsn)}  '_vs_r' L{RSN_I(rsn)}   '_goodind_allone_corr']);
%     errorbar(1:5,cRSN.b_obs_hc_goodind(matchOrd),cRSN.b_obs_hc_goodind(matchOrd)./cRSN.t_obs_hc_goodind(matchOrd),'k');xlim([0.5 5.5]); set(gca,'XTick',[1:5])
%     hold on;errorbar(1:5,cRSN.b_obs_sz_goodind(matchOrd),cRSN.b_obs_sz_goodind(matchOrd)./cRSN.t_obs_sz_goodind(matchOrd),'r');xlim([0.5 5.5]); set(gca,'XTick',[1:5])
%     plot(1:5,[0 0 0 0 0],'k--')
%     export_fig(fullfile(FIGoutputdir,[get(F1,'Name') '.pdf']),'-pdf',F1)
end
%%
for rsn = [  16  23  47] %[1 5 7 16  23 24 29 47]
    
    
    
    
    Fa = figure;set(Fa,'Color',[1 1 1], 'Position', CPOS2); hold on
    set(Fa,'Name',[L{RSN_I(rsn)}  'Specra_byState_corr'])
    %load spectra_thal_by_state_kmeans_city_k5_thr30
    
    for state = 1:5
        state1 = matched_citystate(state);
        
        statehcid = find(aINDallone_HC(:) == state);
        stateszid = find(aINDallone_SZ(:) == state);        
        
        hcdata =  squeeze(dynSpec_nt_nw_HC(statehcid,rsn,:));
        szdata =  squeeze(dynSpec_nt_nw_SZ(stateszid,rsn,:));        
        
        %if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8)
        %    [hhsdyn(state,:) ppsdyn(state,:)] = ttest2(hcdata,szdata);
            %figure;plot(sortrows([ppp' myFDR(ppp)']))
            %title(['state ' num2str(state)])
            
        %end
        Fs = subplot(3,2,state1);
        plot_with_ste_area(Fs,f,hcdata);hold on
        plot_with_ste_area(Fs,f,szdata,[],'b','c');
        title(['state ' num2str(state1)]);
        set(Fs,'FontSize',14)
        ylabel('amplitude A.U.')
        xlabel('frequency')
        ylim([0 9])
        CH = get(Fs,'children');
        legend(CH([2 5]),{'SZ' 'HC'});
        
    end
    
    export_fig(fullfile(FIGoutputdir, [get(Fa, 'Name') '_allwin.pdf']), Fa);
    close(Fa) 
    
   
    F = figure; set(F,'Color',[1 1 1], 'Position', CPOS3);
    set(F,'Name',[L{RSN_I(rsn)}  'Specra_byStateGroup_corr'])
    %load spectra_thal_by_state_kmeans_city_k5_thr30
    clrf = [0 0 0;0 0 1;1 0 0; 0.6 1 0.2; 0.3333    0.0196    0.3686];
    clrb = [.5 .5 .5; .5373 .8118 .9412;1 1 0;0 1 0;.929  .2235  .98];
    clrf = clrf([2 3 1 5 4],:);clrb = clrb([2 3 1 5 4],:);
    for state1 = 1:5
        state = matched_citystate(state1);
        statehcid = find(aINDallone_HC(:) == state);
        stateszid = find(aINDallone_SZ(:) == state);        
        
        hcdata =  squeeze(dynSpec_nt_nw_HC(statehcid,rsn,:));
        szdata =  squeeze(dynSpec_nt_nw_SZ(stateszid,rsn,:));
        
        %if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8)
        %    [hhh(state,:) ppp(state,:)] = ttest2(hcdata,szdata);
        %    %figure;plot(sortrows([ppp' myFDR(ppp)']))
        %    %title(['state ' num2str(state)])
        %end
        
        Fs1 = subplot(1,2,1);hold on
        plot_with_ste_area(Fs1,f,hcdata,[],clrf(state,:),clrb(state,:));
        ylim([0  9])
        if state1 == 5            
            title(['HC: ' L{RSN_I(rsn)} ' Mean Power Spectra']);
            CH1 = get(Fs1,'children');
            legend(CH1([2:3:14]),{'state5' 'state4' 'state3' 'state2' 'state1'});
            set(Fs1,'FontSize',14)
            ylabel('amplitude A.U.')
            xlabel('frequency')
        end

        Fs2 = subplot(1,2,2); hold on
        plot_with_ste_area(Fs2,f,szdata,[],clrf(state,:),clrb(state,:));
        ylim([0 9])
        if state1 == 5
            title(['SZ: ' L{RSN_I(rsn)} ' Mean Power Spectra']);
            CH2 = get(Fs2,'children');
            legend(CH2([2:3:14]),{'state5' 'state4' 'state3' 'state2' 'state1'});
            set(Fs2,'FontSize',14)
            ylabel('amplitude A.U.')
            xlabel('frequency')
        end

        % CH = get(Fs,'children');
        % legend(CH([2 5]),{'SZ' 'HC'});
    end
    
    export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '_allwin.pdf']), F);
    close(F)
    clear hcdata szdata
end
%%  tapered

for rsn = [1 5 7 16  23 24 29 47]   
    
    Fa = figure;set(Fa,'Color',[1 1 1], 'Position', CPOS2); hold on
    set(Fa,'Name',[L{RSN_I(rsn)}  'Specra_byState_corr'])
    %load spectra_thal_by_state_kmeans_city_k5_thr30
    
    for state = 1:5
        state1 = matched_citystate(state);
        
        statehcid = find(aINDallone_HC(:) == state);
        stateszid = find(aINDallone_SZ(:) == state);        
        
        hcdata =  squeeze(dynSpec_t_nw_HC(statehcid,rsn,:));
        szdata =  squeeze(dynSpec_t_nw_SZ(stateszid,rsn,:));        
        
        %if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8)
        %    [hhsdyn(state,:) ppsdyn(state,:)] = ttest2(hcdata,szdata);
            %figure;plot(sortrows([ppp' myFDR(ppp)']))
            %title(['state ' num2str(state)])
            
        %end
        Fs = subplot(3,2,state1);
        plot_with_ste_area(Fs,f,hcdata);hold on
        plot_with_ste_area(Fs,f,szdata,[],'b','c');
        title(['state ' num2str(state1)]);
        set(Fs,'FontSize',14)
        ylabel('amplitude A.U.')
        xlabel('frequency')
        ylim([0 0.6])
        CH = get(Fs,'children');
        legend(CH([2 5]),{'SZ' 'HC'});
        
    end
    
    export_fig(fullfile(FIGoutputdir, [get(Fa, 'Name') '_allwint.pdf']), Fa);
    %close(Fa) 
    
   
    F = figure; set(F,'Color',[1 1 1], 'Position', CPOS3);
    set(F,'Name',[L{RSN_I(rsn)}  'Specra_byStateGroup_corr'])
    %load spectra_thal_by_state_kmeans_city_k5_thr30
    clrf = [0 0 0;0 0 1;1 0 0; 0.6 1 0.2; 0.3333    0.0196    0.3686];
    clrb = [.5 .5 .5; .5373 .8118 .9412;1 1 0;0 1 0;.929  .2235  .98];
    clrf = clrf([2 3 1 5 4],:);clrb = clrb([2 3 1 5 4],:);
    for state1 = 1:5
        state = matched_citystate(state1);
        statehcid = find(aINDallone_HC(:) == state);
        stateszid = find(aINDallone_SZ(:) == state);        
        
        hcdata =  squeeze(dynSpec_t_nw_HC(statehcid,rsn,:));
        szdata =  squeeze(dynSpec_t_nw_SZ(stateszid,rsn,:));
        
        %if (length(find(tf1 == 1)) > 8) && (length(find(tf2 == 1)) > 8)
        %    [hhh(state,:) ppp(state,:)] = ttest2(hcdata,szdata);
        %    %figure;plot(sortrows([ppp' myFDR(ppp)']))
        %    %title(['state ' num2str(state)])
        %end
        
        Fs1 = subplot(1,2,1);hold on
        plot_with_ste_area(Fs1,f,hcdata,[],clrf(state,:),clrb(state,:));
        ylim([0  0.6])
        if state1 == 5            
            title(['HC: ' L{RSN_I(rsn)} ' Mean Power Spectra']);
            CH1 = get(Fs1,'children');
            legend(CH1([2:3:14]),{'state5' 'state4' 'state3' 'state2' 'state1'});
            set(Fs1,'FontSize',14)
            ylabel('amplitude A.U.')
            xlabel('frequency')
        end

        Fs2 = subplot(1,2,2); hold on
        plot_with_ste_area(Fs2,f,szdata,[],clrf(state,:),clrb(state,:));
        ylim([0 0.6])
        if state1 == 5
            title(['SZ: ' L{RSN_I(rsn)} ' Mean Power Spectra']);
            CH2 = get(Fs2,'children');
            legend(CH2([2:3:14]),{'state5' 'state4' 'state3' 'state2' 'state1'});
            set(Fs2,'FontSize',14)
            ylabel('amplitude A.U.')
            xlabel('frequency')
        end

        % CH = get(Fs,'children');
        % legend(CH([2 5]),{'SZ' 'HC'});
    end
    
    export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '_allwint.pdf']), F);
    %close(F)
    clear hcdata szdata
end
%%  kmedoids

for k2 = 2:9
   
    [label, energy, index] = kmedoids(Fdynflat' ,k2);
   
    save(fullfile(fileparts(inname), ['FNC_kmedclusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew']), 'label','energy', 'index')
    
    label = reshape(label,M,Nwin);
    label_HC = label(HCid,:);
    label_SZ = label(SZid,:);
    
    G = figure('Color', 'w', 'Name', ['FNC_kmedclusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew'], 'Position', CPOS);
    for ii = 1:k2,
        subplot(2,k2,ii);
        kmedh = vec2mat(median(FdynflatHC(label_HC == ii ,:)),1);
        imagesc(kmedh(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(sum(label_HC == ii)),  round(100*sum(sum(label_HC == ii))/length(label_HC(:)))));
        
        subplot(2,k2,ii+k2);
        kmeds = vec2mat(median(FdynflatSZ(label_SZ == ii ,:)),1);
        imagesc(kmeds(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
        set(gca, 'XTick', [], 'YTick', [])
        title(sprintf('%d (%d%%)', sum(sum(label_SZ == ii)),  round(100*sum(sum(label_SZ == ii))/length(label_SZ(:)))));
    end
    IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);
end
%%

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
addpath(genpath('/export/mialab/hcp/rest/analysis/CODE_dyn/classification/distributionPlot/'))
k2 = 5;
all1_k5 = load(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_corr_scrubbed_cnew']));
aFR = zeros(M,k2);
aTM = zeros(M,k2,k2);
aMDT = zeros(M,k2);
aNT = zeros(M,1);

for ii = 1:M
    [FRii, TMii, MDTii, NTii] = statevector_stats(all1_k5.aINDallone(ii,:)', k2);
    aFR(ii,:) = FRii;
    aTM(ii,:,:) = TMii;
    aMDT(ii,:) = MDTii;
    aNT(ii) = NTii;
end

G = figure('Color', 'w', 'Name', 'Frequency');
edges = linspace(0,1,20);
for ii = 1:k2
    subplot(k2,1,ii);
    N=histc(aFR(HCid, ii), edges);
    bar(edges,N, 'histc');
    box off;
    set(gca, 'xlim', [0,1], 'ylim', [0 200]);
end
H = figure('Color', 'w', 'Name', 'Frequency');
edges = linspace(0,1,20);
for ii = 1:k2
    subplot(k2,1,ii);
    N=histc(aFR(SZid, ii), edges);
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

%%
matchOrd = [2 3 1 5 4];

aFR = aFR(:,matchOrd);
aTM = aTM(:,matchOrd,matchOrd);
aMDT = aMDT(:,matchOrd);

I = plot_statevector_stats(k2, aFR, aTM, aMDT, aNT);
set(I, 'Name', 'Average_STATESTATS_corr')
IM = export_fig(fullfile(FIGoutputdir, [get(I, 'Name') '.pdf']), I);
I1 = plot_statevector_stats(k2, aFR(HCid,:), aTM(HCid,:,:), aMDT(HCid,:), aNT(HCid,:));
set(I1, 'Name', 'Average_STATESTATS_HC_corr')
IM = export_fig(fullfile(FIGoutputdir, [get(I1, 'Name') '.pdf']), I1);
I2 = plot_statevector_stats(k2, aFR(SZid,:), aTM(SZid,:,:), aMDT(SZid,:), aNT(SZid,:));
set(I2, 'Name', 'Average_STATESTATS_SZ_corr')
IM = export_fig(fullfile(FIGoutputdir, [get(I2, 'Name') '.pdf']), I2);
 A = squeeze(mean(aTM));
 statvec = compute_statprob(A);
%%
[hMDT pMDT ciMDT tstMDT] = ttest2(aMDT(HCid,:),aMDT(SZid,:));
[hFR pFR ciFR tstFR] = ttest2(aFR(HCid,:),aFR(SZid,:));
I3 = plot_statevector_stats(k2, aFR, aTM, aMDT, aNT,ica_ursi_group_keepID);
set(I3, 'Name', 'Average_STATESTATS_HCSZ_corr')
IM = export_fig(fullfile(FIGoutputdir, [get(I3, 'Name') '.pdf']), I3);
%%

for site = 1:7
    
    sii = eval(['ica_ursi_site' num2str(site) '_ID' ]);
    %hcsid = find(strcmp(grp,['HC_s' num2str(site)])==1);
    %szsid = find(strcmp(grp,['SZ_s' num2str(site)])==1);
  %  subplot(4,2,site)
    plot_statevector_stats(k2, aFR(sii,:), aTM(sii,:,:), aMDT(sii,:), aNT(sii),ica_ursi_group_keepID(sii));
    set(gcf, 'Name', ['Average_STATESTATS_HCSZ_site' num2str(site) '_corr_newcmap'])
    C = colorbar;
    ylm = get(C,'YLim');
    set(C,'YLim', [floor(ylm(1)) 0]);
    set(C,'YTick', floor(ylm(1)):1:0,'YTickLabel',10.^(floor(ylm(1)):1:0))
    export_fig(fullfile(FIGoutputdir, [get(gcf, 'Name') '.pdf']), gcf);
end
%close all
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





% %% age difference - stat prob vector
% 
% minage = 17;
% maxage = 24;
% 
% Af = squeeze(mean(aTM(a < minage,:,:)));
% stvecYOUNG = compute_statprob(Af);
% 
% Am = squeeze(mean(aTM(a > maxage,:,:)));
% stvecOLD = compute_statprob(Am);
% gdist = sqrt(sum((stvecYOUNG-stvecOLD).^2));
% adist = abs(stvecYOUNG-stvecOLD);
% %gtdist = sqrt(sum(sum((Af-Am).^2))); 
% 
% aTMr = aTM((a>maxage | a<minage),:,:);
% ar = a(a>maxage | a<minage);
% 
% nboot = 100;
% stvecOLDboot = zeros(nboot,k2);
% stvecYOUNGboot = zeros(nboot,k2);
% for ii = 1:nboot
%     pickme = ceil(length(ar)*rand(1,length(ar)));
%     ab = ar(pickme);
%     aTMb = aTMr(pickme,:,:);
%     Afboot = squeeze(mean(aTMb(ab < minage,:,:)));
%     stvecfboot = compute_statprob(Afboot);
%     stvecYOUNGboot(ii,:) = stvecfboot;
%     Amboot = squeeze(mean(aTMb(ab > maxage,:,:)));
%     stvecmboot = compute_statprob(Amboot);
%     stvecOLDboot(ii,:) = stvecmboot;
% end
% otext = sprintf('older subjects (%d to %d, n = %d)', maxage, max(a), sum(a>maxage));
% ytext = sprintf('younger subjects (%d to %d, n = %d)', min(a), minage, sum(a<minage));
% 
% F = figure('Color', 'w', 'Name', 'StatProbVec_AGE', 'Position', [680   762   617   336]);
% Og = plot(1:k2, stvecOLDboot); set(Og, 'Color', [.7 .7 .7])
% hold on
% Ogl = plot(1:k2, stvecOLD, 'k', 'LineWidth', 2);
% Yg = plot(1:k2, stvecYOUNGboot); set(Yg, 'Color', [1 .7 .7])
% hold on
% Ygl = plot(1:k2, stvecYOUNG, 'r', 'LineWidth', 2);
% set(gca, 'TickDir', 'out'); box off
% legend([Ogl, Ygl], {otext,ytext}, 'Location', 'Northwest')
% xlabel('State (cluster index)')
% ylabel('Stationary Probability')
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% aTMr = aTM((a>maxage | a<minage),:,:);
% ar = a(a>maxage | a<minage);
% nperm = 1000;
% gdistnull = zeros(1,nboot);
% adistnull = zeros(k2,nboot);
% for ii = 1:nperm
%     aboot = ar(randperm(length(ar)));
%     Afboot = squeeze(mean(aTMr(aboot < minage,:,:)));
%     stvecfboot = compute_statprob(Afboot);
%     Amboot = squeeze(mean(aTMr(aboot > maxage,:,:)));
%     stvecmboot = compute_statprob(Amboot);
%     gdistnull(ii) = sqrt(sum((stvecfboot-stvecmboot).^2)); 
%     adistnull(:,ii) = abs(stvecfboot-stvecmboot);
%     %gtdistnull(ii) = sqrt(sum(sum((Afboot-Amboot).^2))); 
% end
% 
% F = figure('Color', 'w', 'Name', 'StatProbVec_AGE_Distance', 'Position', [100         375        490         300]);
% [N,bc] = hist(gdistnull, 75);
% B = bar(bc, N, 1); box off; 
% set(B, 'EdgeColor', 'none', 'FaceColor', 'k')
% set(gca, 'TickDir', 'out')
% xlabel('Euclidean distance between state probability vectors')
% ylabel('Frequency')
% hold on
% plot(gdist, 10, 'r*')
% TO = text(gdist, 100, 'observed'); set(TO, 'Color', 'r', 'HorizontalAlignment', 'center')
% TN = text(mean(gdistnull), max(N)+30, 'null'); set(TN, 'Color', 'k', 'HorizontalAlignment', 'center');
% %axis tight;
% %drawaxismargins(gca, .1,0)
% title(sprintf('Z(obs) = %0.1f, P = %0.3f', (gdist-mean(gdistnull))/std(gdistnull), sum(gdistnull>gdist)/nperm))
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% 
% 
% %% distribution of ages
% F = figure('Color', 'w', 'Name', 'SubjectAge', 'Position', [680   762   350   250]);
% hist(a, min(a):max(a)); box off; set(gca, 'TickDir', 'out'); axis tight
% drawaxismargins(gca, .05,0)
% xlabel('Age')
% ylabel('Frequency')
% ax = axis;
% hold on
% plot([minage, minage], ax(3:4), 'r--')
% plot([maxage, maxage], ax(3:4), 'r--')
% IM = export_fig(fullfile(FIGoutputdir, [get(F, 'Name') '.pdf']), F);
% 
% %% gender difference/age differences for statvector not found -- come back to this
% 
% g_o = loadlabel('gender');
% a = loadlabel('age');
% 
% aTMr = aTM;
% g = g_o;
% Af = squeeze(mean(aTMr(g == 1,:,:)));
% stvecf = compute_statprob(Af);
% 
% Am = squeeze(mean(aTMr(g == -1,:,:)));
% stvecm = compute_statprob(Am);
% %gdist =  sum(stvecm.*log(stvecm./stvecf));
% gdist = sqrt(sum((stvecf-stvecm).^2));
% 
% nboot = 1000;
% gdistnull = zeros(1,nboot);
% for ii = 1:nboot
%     gboot = g(randperm(length(g)));
%     Afboot = squeeze(mean(aTMr(gboot == 1,:,:)));
%     stvecfboot = compute_statprob(Afboot);
%     Amboot = squeeze(mean(aTMr(gboot == -1,:,:)));
%     stvecmboot = compute_statprob(Amboot);
%     gdistnull(ii) = sqrt(sum((stvecfboot-stvecmboot).^2));
%     %gdistnull(ii) =  sum(stvecmboot.*log(stvecmboot./stvecfboot));
% end