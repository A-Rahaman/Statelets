addpath(genpath('/export/mialab/users/eswar/fbirn_p3/scripts/cluster'))
addpath(genpath('/export/mialab/users/eswar/fbirn_p3/scripts/network_estimation'))
addpath(genpath('/export/mialab/users/eswar/software/nic'))
%inname = '/export/mialab/hcp/dynamics/FNC_dynamicSWTukey18_L1ICOV_zscored';
%FIGoutputdir = '/export/mialab/hcp/dynamics/figures/L1';

FIGoutputdir = '/export/mialab/users/eswar/fbirn_p3/results_C100_vn/figures/dyn_L1';
inname = '/export/mialab/users/eswar/fbirn_p3/results_C100_vn/FNCdynamics/FNC_dynamicSWTukey22_L1ICOV_zscored_scr_cnew';
load orderedRSN_by_modularity
TIMEPOS = [ 680   752   497   346];
CPOS = [ 180         886        1717         212];
load(inname)
% loads FNCdynflat [Nwin, M, C*(C-1)/2]
getICA_SubjGroupID;
FNCdynflat = FNCdynflat(:,keepID,:);
[Nwin, M, nCC] = size(FNCdynflat);
fd_dv_orig = load('FD_DVARS_fbirn_p3_ICA');
DVARS_icasubj = fd_dv_orig.dvars_v_ica(keepID,:);
FD_icasubj = fd_dv_orig.FD_ica(keepID,:);
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

%% make the sliding window
nT = 158;
sigma = 3; % for smoothing Gaussian
wsize = 22; % for box (44 seconds)
gw = gaussianwindow(nT,nT/2,sigma);
b = zeros(nT,1);  b((nT/2 - wsize/2 + 1):(nT/2+wsize/2)) = 1;
c = conv(gw,b); c = c/max(c); c=c(nT/2+1:end-nT/2 +1);

A = repmat(c,1,length(RSN_I));
nTwin = length(find(c > 1e-4));
tp_win = zeros(Nwin,nTwin);
for ii=1:Nwin
	Ashift =circshift(A,-nT/2 + wsize/2 + ii);
    cshift = circshift(c,-nT/2 + wsize/2 + ii);
    csind = find(cshift > 1e-4)';
    tp_win(ii,:) = csind;        
end
%%
p=3; % level of detrending
c1 = 2.5; % cutoff
keepWin = cell(M,1);keepWin1 = cell(M,1);
keepWin2 = cell(M,1);keepWin3 = cell(M,1);
keepWin4 = cell(M,1);keepWin5 = cell(M,1);
keepWin6 = cell(M,1);
for ss = 1:M
    DVARS_subj = DVARS_icasubj(ss,4:end);
    r = length(DVARS_subj);
    b = ((1 : r)' * ones (1, p + 1)) .^ (ones (r, 1) * (0 : p));  % build the regressors
    lestimates = robustfit(b(:,2:end), DVARS_subj);
    yfit = b*lestimates;
    res = DVARS_subj' - yfit;     
    mad_res = median(abs(res - median(res))); % median absolute deviation of residuals
    sigma = mad_res* sqrt(pi/2);
    s = res/sigma;
    ind = find(abs(s) > c1)+1;
    
    FD_subj = FD_icasubj(ss,4:end);
    keepWin4{ss} = find(FD_subj < 0.3)+1;
    keepWin5{ss} = find(FD_subj < 0.4)+1;
    keepWin6{ss} = find(FD_subj < 0.5)+1;
    lestimates1 = robustfit(b(:,2:end), FD_subj);
    yfit1 = b*lestimates1;
    res1 = FD_subj' - yfit1;     
    mad_res1 = median(abs(res1 - median(res1))); % median absolute deviation of residuals
    sigma1 = mad_res1* sqrt(pi/2);
    s1 = res1/sigma1;
    ind1 = find(abs(s1) > c1)+1;
    tmp = [];tmp1 = [];tmp2 = [];tmp3 = [];
    for ii = 1:Nwin
        if ~any(ismember(ind,tp_win(ii,:)))
            tmp = [tmp ii];
        end
        if sum(ismember(ind,tp_win(ii,:))) < 2
            tmp1 = [tmp1 ii];
        end
        if ~any(ismember(ind1,tp_win(ii,:)))
            tmp2 = [tmp2 ii];
        end
        if sum(ismember(ind1,tp_win(ii,:))) < 2
            tmp3 = [tmp3 ii];
        end
    end
    keepWin{ss} = tmp;keepWin1{ss} = tmp1;
    keepWin2{ss} = tmp2;keepWin3{ss} = tmp3;
    clear tmp tmp1 tmp2 tmp3;
end
%%
kwinFDmeanall = cell(M,1); 
for ss = 1:M
    kvec = keepWin3{ss};
    if ~isempty(kvec)
        
        kwinFDmean = zeros(1,length(kvec));
        for ii = 1:length(kvec)
            kwinFDmean(ii) = mean(FD_icasubj(ss,kvec(ii)+2:kvec(ii)+24));
        end
        kwinFDmeanall{ss} = kwinFDmean;
    end
end
kwinFDmeanmeanall   = 999*ones(M,1);  
for ss =  1:M
    if ~isempty(kwinFDmeanall{ss})
        kwinFDmeanmeanall(ss) = mean(kwinFDmeanall{ss});
    end
end
    
%%
 lenall = ones(4,M)*999; 
for ss = 1:M, 
    l1 =  length(keepWin{ss}); lenall(1,ss)=l1;
    l2 =  length(keepWin1{ss});lenall(2,ss)=l2;
    l3 =  length(keepWin2{ss}); lenall(3,ss)=l3;
    l4 =  length(keepWin3{ss});lenall(4,ss)=l4;
end
%%
slenall = sum(lenall');
cslenall = cumsum(lenall')';
Fdyn_good_flat = zeros(slenall(4),size(Fdyn,3));
for ss = 1:M,
    if lenall(4,ss) > 0
        if ss == 1
            Fdyn_good_flat(1:cslenall(4,ss),:) = squeeze(Fdyn(ss, keepWin3{ss},:));
        else
            Fdyn_good_flat(cslenall(4,ss-1)+1:cslenall(4,ss),:) = squeeze(Fdyn(ss, keepWin3{ss},:));
        end
    end
end
%%
k2 = 5;
dmethod = 'city';

load(fullfile(fileparts(inname), ['FNC_group_clusters_k2_' num2str(k2) 'scrubbed_cnew']), 'SPflat', 'IDXp','Cp');
[gIDXallone, gCallone, gSUMDallone, gDallone] = kmeans(Fdyn_good_flat, k2, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', 150, 'empty', 'drop', 'Start', Cp);
save(fullfile(fileparts(inname), ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_scrubbed_cnew_goodFDWin']), 'gIDXallone', 'gCallone','gSUMDallone','gDallone')

HCid = find(ica_ursi_group_keepID == 0);
SZid = find(ica_ursi_group_keepID == 1);
%%
gIDXallone_HC = [];hcwinind = [];
for M1 = 1:length(HCid)
    cs = HCid(M1);
    if lenall(4,cs) > 0
        if cs == 1
            gIDXallone_HC = [gIDXallone_HC;gIDXallone(1:cslenall(4,cs))];
            hcwinind = [hcwinind  1:cslenall(4,cs)];
        else
            gIDXallone_HC = [gIDXallone_HC;gIDXallone(cslenall(4,cs-1)+1:cslenall(4,cs))];
            hcwinind = [hcwinind  cslenall(4,cs-1)+1:cslenall(4,cs)];
        end
    end
end
    
%%
gIDXallone_SZ = [];szwinind = [];
for M2 = 1:length(SZid)
    cs = SZid(M2);
    if lenall(4,cs) > 0
        if cs == 1
            gIDXallone_SZ = [gIDXallone_SZ;gIDXallone(1:cslenall(4,cs))];
            szwinind = [szwinind  1:cslenall(4,cs)];
        else
            gIDXallone_SZ = [gIDXallone_SZ;gIDXallone(cslenall(4,cs-1)+1:cslenall(4,cs))];
            szwinind = [szwinind  cslenall(4,cs-1)+1:cslenall(4,cs)];
        end
    end
end
%%
Fdyn_good_flat_HC = Fdyn_good_flat(hcwinind,:);
Fdyn_good_flat_SZ = Fdyn_good_flat(szwinind,:);
G = figure('Color', 'w', 'Name', ['FNC_clusters_allone_HCSZ_k' num2str(k2) '_scrubbed_cnew_goodFDwin'], 'Position', CPOS);
for ii = 1:k2    
    subplot(2,k2,ii);
    aset1 = vec2mat(median(Fdyn_good_flat_HC(gIDXallone_HC==ii,:)),1);
    imagesc(aset1(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%d (%d%%)', sum(sum(gIDXallone_HC== ii)),  round(100*sum(sum(gIDXallone_HC== ii))/length(gIDXallone_HC(:)))));
    
    subplot(2,k2,ii+k2);
    aset2 = vec2mat(median(Fdyn_good_flat_SZ(gIDXallone_SZ == ii ,:)),1);
    imagesc(aset2(mnOrd_final_ord,mnOrd_final_ord), [-.5,.5]); axis square;
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%d (%d%%)', sum(sum(gIDXallone_SZ== ii)),  round(100*sum(sum(gIDXallone_SZ== ii))/length(gIDXallone_SZ(:)))));
end
IM = export_fig(fullfile(FIGoutputdir, [get(G, 'Name') '.pdf']), G);

%%
