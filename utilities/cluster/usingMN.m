%%  subcortical (SC) antagonism  vs sensorimotor (SM) hyperconnectivity MEDIAN
% SM : all sensory ICNs (auditory, visual, motor)
load(fullfile(fileparts(inname), 'mn_mdn_FNC_clusters_alloneByDiag_kall_scrubbed_cnew_corr_resid1model'))
k2 = 5;
sc_to_sm_corr_HC_allstates = cell(1,k2);
scnot_to_sm_corr_HC_allstates = cell(1,k2);
within_sc_corr_HC_allstates = cell(1,k2);
within_sm_corr_HC_allstates = cell(1,k2);
thal_to_sm_corr_HC_allstates = cell(1,k2);
pput_to_sm_corr_HC_allstates = cell(1,k2);
pput_to_sm_corr_SZ_allstates = cell(1,k2);
sc_to_sm_corr_SZ_allstates = cell(1,k2);
scnot_to_sm_corr_SZ_allstates = cell(1,k2);
within_sc_corr_SZ_allstates = cell(1,k2);
within_sm_corr_SZ_allstates = cell(1,k2);
thal_to_sm_corr_SZ_allstates = cell(1,k2);

r_SC2SMcorr_wSMcorr_HC_obs = zeros(1,k2); p_SC2SMcorr_wSMcorr_HC_obs = ones(1,k2);
r_THAL2SMcorr_wSMcorr_HC_obs = zeros(1,k2); p_THAL2SMcorr_wSMcorr_HC_obs = ones(1,k2);
r_wSCcorr_wSMcorr_HC_obs  = zeros(1,k2); p_wSCcorr_wSMcorr_HC_obs = ones(1,k2);

r_SC2SMcorr_wSMcorr_SZ_obs = zeros(1,k2); p_SC2SMcorr_wSMcorr_SZ_obs = ones(1,k2);
r_THAL2SMcorr_wSMcorr_SZ_obs = zeros(1,k2); p_THAL2SMcorr_wSMcorr_SZ_obs = ones(1,k2);
r_wSCcorr_wSMcorr_SZ_obs  = zeros(1,k2); p_wSCcorr_wSMcorr_SZ_obs = ones(1,k2);
% 
 nBOOT = 1000;
r_SC2SMcorr_wSMcorr_HC_boot = zeros(nBOOT,k2);
r_THAL2SMcorr_wSMcorr_HC_boot = zeros(nBOOT,k2);
r_wSCcorr_wSMcorr_HC_boot = zeros(nBOOT,k2);
r_SC2SMcorr_wSMcorr_SZ_boot = zeros(nBOOT,k2);
r_THAL2SMcorr_wSMcorr_SZ_boot = zeros(nBOOT,k2);
r_wSCcorr_wSMcorr_SZ_boot = zeros(nBOOT,k2);


for state = 1:k2,
    
    mnHC_state = vec2mat(mean_Kwin_allHC{k2-2,state});
    mnHC_state = mnHC_state(:,mnOrd_final_ord,mnOrd_final_ord);
    %figure;imagesc(squeeze(mean(mnHC_state)),[-0.5 0.5])

    mnSZ_state = vec2mat(mean_Kwin_allSZ{k2-2,state});
    mnSZ_state = mnSZ_state(:,mnOrd_final_ord,mnOrd_final_ord);
    
    sc_to_sm_corr_HC = zeros(nS_Kwin_allHC{k2-2}(state),1);
    scnot_to_sm_corr_HC = zeros(nS_Kwin_allHC{k2-2}(state),1);
    within_sc_corr_HC = sc_to_sm_corr_HC;within_sm_corr_HC = sc_to_sm_corr_HC;
    thal_to_sm_corr_HC = sc_to_sm_corr_HC;
    pput_to_sm_corr_HC = zeros(nS_Kwin_allHC{k2-2}(state),1);
    for ii = 1:nS_Kwin_allHC{k2-2}(state)
        sc_to_sm_corr_HC(ii) = sum(sum(squeeze(mnHC_state(ii,1:5,6:24))))/numel(squeeze(mnHC_state(ii,1:5,6:24))); 
        scnot_to_sm_corr_HC(ii) = sum(sum(squeeze(mnHC_state(ii,1:4,6:24))))/numel(squeeze(mnHC_state(ii,1:4,6:24))); 
        within_sc_corr_HC(ii) = 0.5*(sum(sum(squeeze(mnHC_state(ii,1:5,1:5)))))/nchoosek(5,2);
        within_sm_corr_HC(ii) = 0.5*(sum(sum(squeeze(mnHC_state(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
        thal_to_sm_corr_HC(ii) = (sum(squeeze(mnHC_state(ii,5,6:24))))/numel(6:24); 
        pput_to_sm_corr_HC(ii) = (sum(squeeze(mnHC_state(ii,4,6:24))))/numel(6:24);
    end
    
    sc_to_sm_corr_SZ = zeros(nS_Kwin_allSZ{k2-2}(state),1);
    scnot_to_sm_corr_SZ = zeros(nS_Kwin_allSZ{k2-2}(state),1);
    within_sc_corr_SZ = sc_to_sm_corr_SZ;within_sm_corr_SZ = sc_to_sm_corr_SZ;
    thal_to_sm_corr_SZ = sc_to_sm_corr_SZ;
    pput_to_sm_corr_SZ = zeros(nS_Kwin_allSZ{k2-2}(state),1);
    for ii = 1:nS_Kwin_allSZ{k2-2}(state)
        sc_to_sm_corr_SZ(ii) = sum(sum(squeeze(mnSZ_state(ii,1:5,6:24))))/numel(squeeze(mnSZ_state(ii,1:5,6:24))); 
        scnot_to_sm_corr_SZ(ii) = sum(sum(squeeze(mnSZ_state(ii,1:4,6:24))))/numel(squeeze(mnSZ_state(ii,1:4,6:24))); 
        within_sc_corr_SZ(ii) = 0.5*(sum(sum(squeeze(mnSZ_state(ii,1:5,1:5)))))/nchoosek(5,2);
        within_sm_corr_SZ(ii) = 0.5*(sum(sum(squeeze(mnSZ_state(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
        thal_to_sm_corr_SZ(ii) = (sum(squeeze(mnSZ_state(ii,5,6:24))))/numel(6:24); 
        pput_to_sm_corr_SZ(ii) = (sum(squeeze(mnSZ_state(ii,4,6:24))))/numel(6:24);
    end
    
    sc_to_sm_corr_HC_allstates{state} = sc_to_sm_corr_HC;
    scnot_to_sm_corr_HC_allstates{state} = scnot_to_sm_corr_HC;
    within_sc_corr_HC_allstates{state} = within_sc_corr_HC;
    within_sm_corr_HC_allstates{state} = within_sm_corr_HC;
    thal_to_sm_corr_HC_allstates{state} = thal_to_sm_corr_HC;
    pput_to_sm_corr_HC_allstates{state} = pput_to_sm_corr_HC;

    sc_to_sm_corr_SZ_allstates{state} = sc_to_sm_corr_SZ;
    scnot_to_sm_corr_SZ_allstates{state} = scnot_to_sm_corr_SZ;
    within_sc_corr_SZ_allstates{state} = within_sc_corr_SZ;
    within_sm_corr_SZ_allstates{state} = within_sm_corr_SZ;
    thal_to_sm_corr_SZ_allstates{state} = thal_to_sm_corr_SZ;
    pput_to_sm_corr_SZ_allstates{state} = pput_to_sm_corr_SZ;

    [r_SC2SMcorr_wSMcorr_HC_obs(state) p_SC2SMcorr_wSMcorr_HC_obs(state)] =  corr(sc_to_sm_corr_HC_allstates{state}, within_sm_corr_HC_allstates{state});
    [r_THAL2SMcorr_wSMcorr_HC_obs(state) p_THAL2SMcorr_wSMcorr_HC_obs(state)] =  corr(thal_to_sm_corr_HC_allstates{state}, within_sm_corr_HC_allstates{state});
    [r_wSCcorr_wSMcorr_HC_obs(state) p_wSCcorr_wSMcorr_HC_obs(state)] =  corr(within_sc_corr_HC_allstates{state}, within_sm_corr_HC_allstates{state});
    
    [r_SC2SMcorr_wSMcorr_SZ_obs(state) p_SC2SMcorr_wSMcorr_SZ_obs(state)] =  corr(sc_to_sm_corr_SZ_allstates{state}, within_sm_corr_SZ_allstates{state});
    [r_THAL2SMcorr_wSMcorr_SZ_obs(state) p_THAL2SMcorr_wSMcorr_SZ_obs(state)] =  corr(thal_to_sm_corr_SZ_allstates{state}, within_sm_corr_SZ_allstates{state});
    [r_wSCcorr_wSMcorr_SZ_obs(state) p_wSCcorr_wSMcorr_SZ_obs(state)] =  corr(within_sc_corr_SZ_allstates{state}, within_sm_corr_SZ_allstates{state});
% 
    for bb = 1:nBOOT,
        if mod(bb,100) == 1
            disp(['now working on ' num2str(bb) ' of ' num2str(nBOOT) ' : state ' num2str(state)])
        end
        bootHCind = ceil(rand(1,nS_Kwin_allHC{k2-2}(state))*nS_Kwin_allHC{k2-2}(state));        
        bootSZind = ceil(rand(1,nS_Kwin_allSZ{k2-2}(state))*nS_Kwin_allSZ{k2-2}(state));
        
        mnHC_state = vec2mat(mean_Kwin_allHC{k2-2,state});
        mnHC_state = mnHC_state(:,mnOrd_final_ord,mnOrd_final_ord);
        mnHC_state = mnHC_state(bootHCind,:,:);
        mnSZ_state = vec2mat(mean_Kwin_allSZ{k2-2,state});
        mnSZ_state = mnSZ_state(:,mnOrd_final_ord,mnOrd_final_ord);
        mnSZ_state = mnSZ_state(bootSZind,:,:);
        
        for ii = 1:nS_Kwin_allHC{k2-2}(state)
            sc_to_sm_corr_HC(ii) = sum(sum(squeeze(mnHC_state(ii,1:5,6:24))))/numel(squeeze(mnHC_state(ii,1:5,6:24))); 
            within_sc_corr_HC(ii) = 0.5*(sum(sum(squeeze(mnHC_state(ii,1:5,1:5)))))/nchoosek(5,2);
            within_sm_corr_HC(ii) = 0.5*(sum(sum(squeeze(mnHC_state(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
            thal_to_sm_corr_HC(ii) = (sum(squeeze(mnHC_state(ii,5,6:24))))/numel(6:24); 
        end
        
        for ii = 1:nS_Kwin_allSZ{k2-2}(state)
            sc_to_sm_corr_SZ(ii) = sum(sum(squeeze(mnSZ_state(ii,1:5,6:24))))/numel(squeeze(mnSZ_state(ii,1:5,6:24))); 
            within_sc_corr_SZ(ii) = 0.5*(sum(sum(squeeze(mnSZ_state(ii,1:5,1:5)))))/nchoosek(5,2);
            within_sm_corr_SZ(ii) = 0.5*(sum(sum(squeeze(mnSZ_state(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
            thal_to_sm_corr_SZ(ii) = (sum(squeeze(mnSZ_state(ii,5,6:24))))/numel(6:24); 
        end

        
        r_SC2SMcorr_wSMcorr_HC_boot(bb,state) = corr(sc_to_sm_corr_HC,within_sm_corr_HC);
        r_THAL2SMcorr_wSMcorr_HC_boot(bb,state) = corr(thal_to_sm_corr_HC,within_sm_corr_HC);
        r_wSCcorr_wSMcorr_HC_boot(bb,state) = corr(within_sc_corr_HC,within_sm_corr_HC);
        
        r_SC2SMcorr_wSMcorr_SZ_boot(bb,state) = corr(sc_to_sm_corr_SZ,within_sm_corr_SZ);
        r_THAL2SMcorr_wSMcorr_SZ_boot(bb,state) = corr(thal_to_sm_corr_SZ,within_sm_corr_SZ);
        r_wSCcorr_wSMcorr_SZ_boot(bb,state) = corr(within_sc_corr_SZ,within_sm_corr_SZ);
    end
    
    disp(['Done state ' num2str(state)])
end

save(fullfile(fileparts(inname),  'dfnc_SC_SM_relationships_k5_corr_recid1model_mn'),  'r_SC*', 'r_THAL*', 'r_wSC*', 'p_SC*', 'p_THAL*', 'p_wSC*', '*allstates')
%%
load(fullfile(fileparts(inname), 'mn_mdn_FNC_clusters_alloneByDiag_kall_scrubbed_cnew_corr_resid1model'))
X1_HC_SC2SMcorr_recid1corr_boot = zeros(nBOOT,k2);
Y_SZ_thal2SMcorr_recid1corr_boot = zeros(nBOOT,k2);
Y_HC_thal2SMcorr_recid1corr_boot = zeros(nBOOT,k2);
X2_HC_wSMcorr_recid1corr_boot = zeros(nBOOT,k2);
Y_SZ_pput2SMcorr_recid1corr_boot = zeros(nBOOT,k2);
Y_HC_pput2SMcorr_recid1corr_boot = zeros(nBOOT,k2);


for bb = 1:nBOOT,
    sc_to_sm_corr_HC_allstatesb = cell(1,k2);
    scnot_to_sm_corr_HC_allstatesb = cell(1,k2);
    within_sc_corr_HC_allstatesb = cell(1,k2);
    within_sm_corr_HC_allstatesb = cell(1,k2);
    thal_to_sm_corr_HC_allstatesb = cell(1,k2);
    sc_to_sm_corr_SZ_allstatesb = cell(1,k2);
    scnot_to_sm_corr_SZ_allstatesb = cell(1,k2);
    within_sc_corr_SZ_allstatesb = cell(1,k2);
    within_sm_corr_SZ_allstatesb = cell(1,k2);
    thal_to_sm_corr_SZ_allstatesb = cell(1,k2);
    pput_to_sm_corr_HC_allstatesb = cell(1,k2);
    pput_to_sm_corr_SZ_allstatesb = cell(1,k2);
    


    for state = 1:k2
        if mod(bb,100) == 1
            disp(['now working on ' num2str(bb) ' of ' num2str(nBOOT) ' : state ' num2str(state)])
        end
        bootHCind = ceil(rand(1,nS_Kwin_allHC{k2-2}(state))*nS_Kwin_allHC{k2-2}(state));        
        bootSZind = ceil(rand(1,nS_Kwin_allSZ{k2-2}(state))*nS_Kwin_allSZ{k2-2}(state));
        
        sc_to_sm_corr_HCb = zeros(nS_Kwin_allHC{k2-2}(state),1);
        scnot_to_sm_corr_HCb = sc_to_sm_corr_HCb;
        within_sc_corr_HCb = sc_to_sm_corr_HCb;
        within_sm_corr_HCb = sc_to_sm_corr_HCb;
        thal_to_sm_corr_HCb = sc_to_sm_corr_HCb;
        pput_to_sm_corr_HCb = sc_to_sm_corr_HCb;
    
        sc_to_sm_corr_SZb = zeros(nS_Kwin_allSZ{k2-2}(state),1);
        scnot_to_sm_corr_SZb = sc_to_sm_corr_SZb;
        within_sc_corr_SZb = sc_to_sm_corr_SZb;
        within_sm_corr_SZb = sc_to_sm_corr_SZb;
        thal_to_sm_corr_SZb = sc_to_sm_corr_SZb;
        pput_to_sm_corr_SZb = sc_to_sm_corr_SZb;
        
        mnHC_stateb = vec2mat(mean_Kwin_allHC{k2-2,state});
        mnHC_stateb = mnHC_stateb(:,mnOrd_final_ord,mnOrd_final_ord);
        mnHC_stateb = mnHC_stateb(bootHCind,:,:);
        mnSZ_stateb = vec2mat(mean_Kwin_allSZ{k2-2,state});
        mnSZ_stateb = mnSZ_stateb(:,mnOrd_final_ord,mnOrd_final_ord);
        mnSZ_stateb = mnSZ_stateb(bootSZind,:,:);
        
        for ii = 1:nS_Kwin_allHC{k2-2}(state)
            sc_to_sm_corr_HCb(ii) = sum(sum(squeeze(mnHC_stateb(ii,1:5,6:24))))/numel(squeeze(mnHC_stateb(ii,1:5,6:24))); 
            scnot_to_sm_corr_HCb(ii) = sum(sum(squeeze(mnHC_stateb(ii,1:4,6:24))))/numel(squeeze(mnHC_stateb(ii,1:4,6:24))); 
            within_sc_corr_HCb(ii) = 0.5*(sum(sum(squeeze(mnHC_stateb(ii,1:5,1:5)))))/nchoosek(5,2);
            within_sm_corr_HCb(ii) = 0.5*(sum(sum(squeeze(mnHC_stateb(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
            thal_to_sm_corr_HCb(ii) = (sum(squeeze(mnHC_stateb(ii,5,6:24))))/numel(6:24); 
            pput_to_sm_corr_HCb(ii) = (sum(squeeze(mnHC_stateb(ii,4,6:24))))/numel(6:24); 
        end
        
        for ii = 1:nS_Kwin_allSZ{k2-2}(state)
            sc_to_sm_corr_SZb(ii) = sum(sum(squeeze(mnSZ_stateb(ii,1:5,6:24))))/numel(squeeze(mnSZ_stateb(ii,1:5,6:24))); 
            scnot_to_sm_corr_SZb(ii) = sum(sum(squeeze(mnSZ_stateb(ii,1:4,6:24))))/numel(squeeze(mnSZ_stateb(ii,1:4,6:24))); 
            within_sc_corr_SZb(ii) = 0.5*(sum(sum(squeeze(mnSZ_stateb(ii,1:5,1:5)))))/nchoosek(5,2);
            within_sm_corr_SZb(ii) = 0.5*(sum(sum(squeeze(mnSZ_stateb(ii,6:24,6:24)))))/nchoosek(numel(6:24),2);
            thal_to_sm_corr_SZb(ii) = (sum(squeeze(mnSZ_stateb(ii,5,6:24))))/numel(6:24); 
            pput_to_sm_corr_SZb(ii) = (sum(squeeze(mnSZ_stateb(ii,4,6:24))))/numel(6:24);
        end
        sc_to_sm_corr_HC_allstatesb{state} = sc_to_sm_corr_HCb;
        scnot_to_sm_corr_HC_allstatesb{state} = scnot_to_sm_corr_HCb;
        within_sc_corr_HC_allstatesb{state} = within_sc_corr_HCb;
        within_sm_corr_HC_allstatesb{state} = within_sm_corr_HCb;
        thal_to_sm_corr_HC_allstatesb{state} = thal_to_sm_corr_HCb;
        pput_to_sm_corr_HC_allstatesb{state} = pput_to_sm_corr_HCb;
        
        sc_to_sm_corr_SZ_allstatesb{state} = sc_to_sm_corr_SZb;
        scnot_to_sm_corr_SZ_allstatesb{state} = scnot_to_sm_corr_SZb;
        within_sc_corr_SZ_allstatesb{state} = within_sc_corr_SZb;
        within_sm_corr_SZ_allstatesb{state} = within_sm_corr_SZb;
        thal_to_sm_corr_SZ_allstatesb{state} = thal_to_sm_corr_SZb;
        pput_to_sm_corr_SZ_allstatesb{state} = pput_to_sm_corr_SZb;
    end
    X1_HC_SC2SMcorr_recid1corr_boot(bb,:) = cellfun(@mean,sc_to_sm_corr_HC_allstatesb);
    Y_SZ_thal2SMcorr_recid1corr_boot(bb,:) =  cellfun(@mean, thal_to_sm_corr_SZ_allstatesb);
    Y_HC_thal2SMcorr_recid1corr_boot(bb,:)  =  cellfun(@mean, thal_to_sm_corr_HC_allstatesb);
    X2_HC_wSMcorr_recid1corr_boot(bb,:) = cellfun(@mean,within_sm_corr_HC_allstatesb);
    Y_SZ_pput2SMcorr_recid1corr_boot(bb,:) =  cellfun(@mean, pput_to_sm_corr_SZ_allstatesb);
    Y_HC_pput2SMcorr_recid1corr_boot(bb,:)  =  cellfun(@mean, pput_to_sm_corr_HC_allstatesb);
    clear *_allstatesb
end
save(fullfile(fileparts(inname),  'dfnc_SC_SM_relationships_k5_corr_recid1model_mn_boot'),'*_boot')
%%
load(fullfile(fileparts(inname),  'dfnc_SC_SM_relationships_k5_corr_recid1model_mn'))
X1_HC_SC2SMcorr_recid1corr = cellfun(@mean,sc_to_sm_corr_HC_allstates);
Y_SZ_thal2SMcorr_recid1corr =  cellfun(@mean, thal_to_sm_corr_SZ_allstates);
Y_HC_thal2SMcorr_recid1corr  =  cellfun(@mean, thal_to_sm_corr_HC_allstates);
Y_SZ_pput2SMcorr_recid1corr =  cellfun(@mean, pput_to_sm_corr_SZ_allstates);
Y_HC_pput2SMcorr_recid1corr  =  cellfun(@mean, pput_to_sm_corr_HC_allstates);
X2_HC_wSMcorr_recid1corr = cellfun(@mean,within_sm_corr_HC_allstates);
lab =  {'S5', 'S3', 'S4', 'S2', 'S1'};

F = figure; set(F,'color',[1 1 1]);
plot(X1_HC_SC2SMcorr_recid1corr, Y_SZ_thal2SMcorr_recid1corr-Y_HC_thal2SMcorr_recid1corr, 'go'); text(X1_HC_SC2SMcorr_recid1corr, Y_SZ_thal2SMcorr_recid1corr-Y_HC_thal2SMcorr_recid1corr, lab)
xlabel('HC: mean subcortical to sensory FC antagonism')
ylabel('SZ - HC: mean thalamus to sensory FC')
set(gca,'FontSize',14)
export_fig(fullfile(FIGoutputdir,'thal2SMgroupdiff_vs_scsmantagonism_mn.pdf'),'-pdf')

F = figure; set(F,'color',[1 1 1]);
plot(X2_HC_wSMcorr_recid1corr, Y_SZ_thal2SMcorr_recid1corr-Y_HC_thal2SMcorr_recid1corr, 'go'); text(X2_HC_wSMcorr_recid1corr, Y_SZ_thal2SMcorr_recid1corr-Y_HC_thal2SMcorr_recid1corr, lab)
xlabel('HC: mean within sensory FC')
ylabel('SZ - HC: mean thalamus to sensory FC')
set(gca,'FontSize',14)
export_fig(fullfile(FIGoutputdir,'thal2SMgroupdiff_vs_withinSM-FC_mn.pdf'),'-pdf')

F = figure; set(F,'color',[1 1 1]);
plot(X1_HC_SC2SMcorr_recid1corr, Y_SZ_pput2SMcorr_recid1corr-Y_HC_pput2SMcorr_recid1corr, 'go'); text(X1_HC_SC2SMcorr_recid1corr, Y_SZ_pput2SMcorr_recid1corr-Y_HC_pput2SMcorr_recid1corr, lab)
xlabel('HC: mean subcortical to sensory FC antagonism')
ylabel('SZ - HC: mean putamen to sensory FC')
set(gca,'FontSize',14)
export_fig(fullfile(FIGoutputdir,'pput2SMgroupdiff_vs_scsmantagonism_mn.pdf'),'-pdf')

F = figure; set(F,'color',[1 1 1]);
plot(X2_HC_wSMcorr_recid1corr, Y_SZ_pput2SMcorr_recid1corr-Y_HC_pput2SMcorr_recid1corr, 'go'); text(X2_HC_wSMcorr_recid1corr, Y_SZ_pput2SMcorr_recid1corr-Y_HC_pput2SMcorr_recid1corr, lab)
xlabel('HC: mean within sensory FC')
ylabel('SZ - HC: mean putamen to sensory FC')
set(gca,'FontSize',14)
export_fig(fullfile(FIGoutputdir,'pput2SMgroupdiff_vs_withinSM-FC_mn.pdf'),'-pdf')

FF = figure;errorbar(mean(X1_HC_SC2SMcorr_recid1corr_boot),mean(Y_SZ_thal2SMcorr_recid1corr_boot - Y_HC_thal2SMcorr_recid1corr_boot),std(Y_SZ_thal2SMcorr_recid1corr_boot - Y_HC_thal2SMcorr_recid1corr_boot),'xr','linestyle','none');
hold on;errorbar(mean(X1_HC_SC2SMcorr_recid1corr_boot),mean(Y_SZ_thal2SMcorr_recid1corr_boot - Y_HC_thal2SMcorr_recid1corr_boot),std(X1_HC_SC2SMcorr_recid1corr_boot),'xk','linestyle','none');
text(mean(X1_HC_SC2SMcorr_recid1corr_boot), mean(Y_SZ_thal2SMcorr_recid1corr_boot-Y_HC_thal2SMcorr_recid1corr_boot), lab)
xlabel('HC:subcortical sensory antagonism','FontSize',14)
ylabel('SZ - HC:thalamus to sensory connectivity','FontSize',14)
set(gca,'FontSize',12)
set(FF,'color',[1 1 1]);
export_fig(fullfile(FIGoutputdir,'dFNC_thalamus2SM-subcortex2SM-relationships_mn.pdf'),'-pdf')

FF1 = figure;errorbar(mean(X2_HC_wSMcorr_recid1corr_boot),mean(Y_SZ_thal2SMcorr_recid1corr_boot - Y_HC_thal2SMcorr_recid1corr_boot),std(Y_SZ_thal2SMcorr_recid1corr_boot - Y_HC_thal2SMcorr_recid1corr_boot),'xr','linestyle','none');
hold on;errorbar(mean(X2_HC_wSMcorr_recid1corr_boot),mean(Y_SZ_thal2SMcorr_recid1corr_boot - Y_HC_thal2SMcorr_recid1corr_boot),std(X2_HC_wSMcorr_recid1corr_boot),'xk','linestyle','none');
text(mean(X2_HC_wSMcorr_recid1corr_boot), mean(Y_SZ_thal2SMcorr_recid1corr_boot-Y_HC_thal2SMcorr_recid1corr_boot), lab)
xlabel('HC:mean within sensory coherence','FontSize',14)
ylabel('SZ - HC:thalamus to sensory connectivity','FontSize',14)
set(gca,'FontSize',12)
set(FF1,'color',[1 1 1]);
export_fig(fullfile(FIGoutputdir,'dFNC_thalamus2SM-wSensory-relationships_mn.pdf'),'-pdf')

FF2 = figure;errorbar(mean(X1_HC_SC2SMcorr_recid1corr_boot),mean(Y_SZ_pput2SMcorr_recid1corr_boot - Y_HC_pput2SMcorr_recid1corr_boot),std(Y_SZ_pput2SMcorr_recid1corr_boot - Y_HC_pput2SMcorr_recid1corr_boot),'xr','linestyle','none');
hold on;errorbar(mean(X1_HC_SC2SMcorr_recid1corr_boot),mean(Y_SZ_pput2SMcorr_recid1corr_boot - Y_HC_pput2SMcorr_recid1corr_boot),std(X1_HC_SC2SMcorr_recid1corr_boot),'xk','linestyle','none');
text(mean(X1_HC_SC2SMcorr_recid1corr_boot), mean(Y_SZ_pput2SMcorr_recid1corr_boot-Y_HC_pput2SMcorr_recid1corr_boot), lab)
xlabel('HC:subcortical sensory antagonism','FontSize',14)
ylabel('SZ - HC:putamen to sensory connectivity','FontSize',14)
set(gca,'FontSize',12)
set(FF2,'color',[1 1 1]);
export_fig(fullfile(FIGoutputdir,'dFNC_putamen2SM-subcortex2SM-relationships_mn.pdf'),'-pdf')

FF3= figure;errorbar(mean(X2_HC_wSMcorr_recid1corr_boot),mean(Y_SZ_pput2SMcorr_recid1corr_boot - Y_HC_pput2SMcorr_recid1corr_boot),std(Y_SZ_pput2SMcorr_recid1corr_boot - Y_HC_pput2SMcorr_recid1corr_boot),'xr','linestyle','none');
hold on;errorbar(mean(X2_HC_wSMcorr_recid1corr_boot),mean(Y_SZ_pput2SMcorr_recid1corr_boot - Y_HC_pput2SMcorr_recid1corr_boot),std(X2_HC_wSMcorr_recid1corr_boot),'xk','linestyle','none');
text(mean(X2_HC_wSMcorr_recid1corr_boot), mean(Y_SZ_pput2SMcorr_recid1corr_boot-Y_HC_pput2SMcorr_recid1corr_boot), lab)
xlabel('HC:mean within sensory coherence','FontSize',14)
ylabel('SZ - HC:putamen to sensory connectivity','FontSize',14)
set(gca,'FontSize',12)
set(FF3,'color',[1 1 1]);
export_fig(fullfile(FIGoutputdir,'dFNC_putamen2SM-wSensory-relationships_mn.pdf'),'-pdf')

%%
matched_citystate = [5 3 4 2 1];
matchOrd = [5 4 2 3 1];

% F = figure;set(F,'color',[1 1 1],'Name','r_SC2SMcorr_wSMcorr_k5_corr','Position',[309         659        1006         206]);
% for st = 1:k2,
%     subplot(1,5,st);
%     % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
%     % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
%     boxplot([r_SC2SMcorr_wSMcorr_HC_boot(:,matchOrd(st));r_SC2SMcorr_wSMcorr_SZ_boot(:,matchOrd(st))],[ones(nBOOT,1);2*ones(nBOOT,1)],'notch','on')
%     title(['state ' num2str(st)])
% end
% 
% F1 = figure;set(F1,'color',[1 1 1],'Name','r_wSCcorr_wSMcorr_k5_recid1corr','Position',[309         659        1006         206]);
% for st = 1:k2,
%     subplot(1,5,st);
%     % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
%     % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
%     boxplot([r_wSCcorr_wSMcorr_HC_boot(:,matchOrd(st));r_wSCcorr_wSMcorr_SZ_boot(:,matchOrd(st))],[ones(nBOOT,1);2*ones(nBOOT,1)],'notch','on')
%     title(['state ' num2str(st)])
% end
% 
% F2 = figure;set(F2,'color',[1 1 1],'Name','r_THAL2SMcorr_wSMcorr_k5_recid1corr','Position',[309         659        1006         206]);
% for st = 1:k2,
%     subplot(1,5,st);
%     % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
%     % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
%     boxplot([r_THAL2SMcorr_wSMcorr_HC_boot(:,matchOrd(st));r_THAL2SMcorr_wSMcorr_SZ_boot(:,matchOrd(st))],[ones(nBOOT,1);2*ones(nBOOT,1)],'notch','on')
%     title(['state ' num2str(st)])
% end

F3 = figure;set(F3,'color',[1 1 1],'Name','scatter_SC2SMcorr_wSMcorr_k5_recid1corr_mn','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(sc_to_sm_corr_HC_allstates{matchOrd(st)},within_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(sc_to_sm_corr_SZ_allstates{matchOrd(st)},within_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    xlim([-0.5 0.5]);ylim([-0.5 0.5]);axis square
    title(['state ' num2str(st)])
end
export_fig(fullfile(FIGoutputdir,[get(F3,'Name') '.pdf']),'-pdf');

F4 = figure;set(F4,'color',[1 1 1],'Name','scatter_wSCcorr_wSMcorr_k5_recid1corr_mn','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(within_sc_corr_HC_allstates{matchOrd(st)},within_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(within_sc_corr_SZ_allstates{matchOrd(st)},within_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    xlim([-0.5 0.5]);ylim([-0.5 0.5]);axis square
    title(['state ' num2str(st)])
end
export_fig(fullfile(FIGoutputdir,[get(F4,'Name') '.pdf']),'-pdf');

F5 = figure;set(F5,'color',[1 1 1],'Name','scatter_THAL2SMcorr_wSMcorr_k5_recid1corr_mn','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(thal_to_sm_corr_HC_allstates{matchOrd(st)},within_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(thal_to_sm_corr_SZ_allstates{matchOrd(st)},within_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    xlim([-0.5 0.5]);ylim([-0.5 0.5]);axis square
    title(['state ' num2str(st)])
end
export_fig(fullfile(FIGoutputdir,[get(F5,'Name') '.pdf']),'-pdf');
[r_SC2SMcorr_wSMcorr_HC_obs(matchOrd);r_SC2SMcorr_wSMcorr_SZ_obs(matchOrd)]
[p_SC2SMcorr_wSMcorr_HC_obs(matchOrd);p_SC2SMcorr_wSMcorr_SZ_obs(matchOrd)]
[r_THAL2SMcorr_wSMcorr_HC_obs(matchOrd);r_THAL2SMcorr_wSMcorr_SZ_obs(matchOrd)]
[p_THAL2SMcorr_wSMcorr_HC_obs(matchOrd);p_THAL2SMcorr_wSMcorr_SZ_obs(matchOrd)]
[r_wSCcorr_wSMcorr_HC_obs(matchOrd);r_wSCcorr_wSMcorr_SZ_obs(matchOrd)]
[p_wSCcorr_wSMcorr_HC_obs(matchOrd);p_wSCcorr_wSMcorr_SZ_obs(matchOrd)]

F3 = figure;set(F3,'color',[1 1 1],'Name','scatter_thal2SM_wSMcorr_k5_recid1corr_mn','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(thal_to_sm_corr_HC_allstates{matchOrd(st)},within_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(thal_to_sm_corr_SZ_allstates{matchOrd(st)},within_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    xlim([-0.5 0.5]);ylim([-0.5 0.5]);axis square
    title(['state ' num2str(st)])
end
F4 = figure;set(F4,'color',[1 1 1],'Name','scatter_thal2SM_SC2SMcorr_k5_recid1corr_mn','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(thal_to_sm_corr_HC_allstates{matchOrd(st)},sc_to_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(thal_to_sm_corr_SZ_allstates{matchOrd(st)},sc_to_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    xlim([-0.5 0.5]);ylim([-0.5 0.5]);axis square
    title(['state ' num2str(st)])
end
F13 = figure;set(F13,'color',[1 1 1],'Name','scatter_pput2SM_wSMcorr_k5_recid1corr_mn','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(pput_to_sm_corr_HC_allstates{matchOrd(st)},within_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(pput_to_sm_corr_SZ_allstates{matchOrd(st)},within_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    xlim([-0.5 0.5]);ylim([-0.5 0.5]);axis square
    title(['state ' num2str(st)])
end
F14 = figure;set(F14,'color',[1 1 1],'Name','scatter_pput2SM_SC2SMcorr_k5_recid1corr_mn','Position',[309         659        1006         206]);
for st = 1:k2,
    subplot(1,5,st);
    % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
    % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
    plot(pput_to_sm_corr_HC_allstates{matchOrd(st)},sc_to_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
    plot(pput_to_sm_corr_SZ_allstates{matchOrd(st)},sc_to_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    xlim([-0.5 0.5]);ylim([-0.5 0.5]);axis square
    title(['state ' num2str(st)])
end

F15 = figure;set(F15,'color',[1 1 1],'Name','scatter_pput2SMcorr_thal2SMcorr_k5_recid1corr_mn','Position',[309         659        1006         206]);
for st = 1:k2,
%     subplot(1,5,st);
%     % scatter(ones(nBOOT,1),r_SC2SMcorr_wSMcorr_HC_boot(:,st));hold;
%     % scatter(1.1*ones(nBOOT,1),r_SC2SMcorr_wSMcorr_SZ_boot(:,st),'r');xlim([0.9 1.2]);
%     plot(thal_to_sm_corr_HC_allstates{matchOrd(st)},pput_to_sm_corr_HC_allstates{matchOrd(st)},'ko');lsline;hold;
%     plot(thal_to_sm_corr_SZ_allstates{matchOrd(st)},pput_to_sm_corr_SZ_allstates{matchOrd(st)},'ro');lsline;
    [rrh pph] = corr(thal_to_sm_corr_HC_allstates{matchOrd(st)},pput_to_sm_corr_HC_allstates{matchOrd(st)});  
    [rrs pps] = corr(thal_to_sm_corr_SZ_allstates{matchOrd(st)},pput_to_sm_corr_SZ_allstates{matchOrd(st)});
    disp([rrh rrs]);disp([pph pps])
    
%     xlim([-0.5 0.5]);ylim([-0.5 0.5]);axis square
%     title(['state ' num2str(st)])
end