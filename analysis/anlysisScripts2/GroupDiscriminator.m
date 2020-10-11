
%realPeaks = peaks_i_j;
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis'))
clear
clc
load allPairshapes_HC.mat
load allPairshapes_SZ.mat
load allPairshapesPD_HC.mat
load allPairshapesPD_SZ.mat
load Global_SZ_peaks_real.mat
SZ_peaks = realPeaks;
clear realPeaks
load Global_HC_peaks_real.mat
HC_peaks=realPeaks;
clear realPeaks

SZP = p_SZ(SZ_peaks);
HCP = p_HC(HC_peaks);
[sSZP,sSZPIDX] = sort(SZP,'descend');
SZ10 = SZ_peaks(sSZPIDX(1:48));
[sHCP,sHCPIDX] = sort(HCP,'descend');
HC10 = HC_peaks(sHCPIDX(1:48));
%% get the pairs with average PD. since one pair can show multiple shapes 
%Pair_wise_PD = zeros(1081,1);

for i =1:size(SZ_allPairshapes,1)
sz_shapes_pairs(i) = SZ_allPairshapes{i,4};
end

for i =1:size(HC_allPairshapes,1)
hc_shapes_pairs(i) = HC_allPairshapes{i,4};
end

for i = 1:1081
    % HC
    hcpairsPD = find(hc_shapes_pairs==i);
    [~,idx]  = max(p_HC(hcpairsPD)); 
    hc_shape_p_i(i)= hcpairsPD(idx);
    hcPDs(i) = sum(p_HC(hcpairsPD));
    % SZ
    szpairsPD = find(sz_shapes_pairs==i);
    [~,idx1]  = max(p_SZ(szpairsPD)); 
    sz_shape_p_i(i)= szpairsPD(idx1);
    szPDs(i) = sum(p_SZ(szpairsPD));
end
[~,HC_Pair_wise_PD] = sort(hcPDs,'descend');
[~,SZ_Pair_wise_PD] = sort(szPDs,'descend');

topHCP = HC_Pair_wise_PD(1:10)
topSZP = SZ_Pair_wise_PD(1:10)
tHCP = hc_shape_p_i(topHCP);
tSZP = sz_shape_p_i(topSZP);
% plotting
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/data'))
load Comp2DomMapping.mat
load Paircoordinates.mat
load Name_of_components.mat
figure()
for i =1:length(tSZP)
    subplot(5,2,i)
    plot(SZ_allPairshapes{(tSZP(i)),1},'r','LineWidth',2);
    hold on
     
    plot(HC_allPairshapes{(tHCP(i)),1},'b','LineWidth',2);
    plot(zeros(length(SZ_allPairshapes{(tSZP(i)),1}),1),'k','LineWidth',2)
    %title(['S: ' num2str(SZ_allPairshapes{tSZP(i),4}) '           H: ' num2str(HC_allPairshapes{tHCP(i),4})])
    HC_comps = coord(HC_allPairshapes{tHCP(i),4},:);
    SZ_comps = coord(SZ_allPairshapes{tSZP(i),4},:);
    title({['S: ' components_47{SZ_comps(1)} ' , ' components_47{SZ_comps(2)}];[ 'H: '  components_47{HC_comps(1)} ' , ' components_47{HC_comps(2)}]})
    
%     HCDom(i,1) = Comp2Dom(HC_comps(1));
%     HCDom(i,2) = Comp2Dom(HC_comps(2));
%     SZ_comps = coord(SZ_allPairshapes{tSZP(i),4},:);
%     SZDom(i,1) = Comp2Dom(SZ_comps(1));
%     SZDom(i,2) = Comp2Dom(SZ_comps(2));
    hold on
    set(gca,'fontsize', 10);
    box off
    axis tight
    axis off 
end
% hold on
% for i =1:length(tHCP)
%     subplot(2,10,i+10)
%     plot(HC_allPairshapes{(tHCP(i)),1},'b','LineWidth',2);
%     hold on
%     plot(zeros(length(HC_allPairshapes{(tHCP(i)),1}),1),'k','LineWidth',2) 
%     set(gca,'fontsize', 12);
%     box off
%     axis tight
%     axis off 
% end

%% Filtering
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/shapeletAlgorithmScripts'))
K_simscores = []; 
idss        = [];
K_simscores(1) = SZP(1); 
idss(1)        = SZ10(1);
%emd = 190;
pals(1) =1; 
for j = 2:length(SZ10)
        flag = 1;
        kx = SZ_allPairshapes{SZ10(j),2}; 
        % already normalized
        %kx = kx-min(kx);
        %kx = kx/sum(kx);
        
        for l = 1:length(idss)
            ky = SZ_allPairshapes{SZ10(l),2}; 
            %ky = ky-min(ky);
            %ky = ky/sum(ky);
            %dEMD(kx,ky)
            if(corr(kx,ky,'Type','Spearman')> 0.7)
            fprintf(" EMD between%d and idss %d is %f\n",j,l,corr(kx,ky,'Type','Spearman'));        
            %if(dEMD(kx,ky) <= emd)
            %fprintf(" EMD between%d and idss %d is %f\n",j,l,dEMD(kx,ky));    
            flag=0;
            break;
            end
        end
        if(flag~=0)
            K_simscores(end+1) = SZP(j); 
            idss(end+1)        = SZ10(j);
            pals(end+1)        = j;
        end
        
end

K_simscores = []; 
idss_hc        = [];
K_simscores(1) = HCP(1); 
idss_hc(1)        = HC10(1);
%emd = 190;
pals(1) =1; 
for j = 2:length(HC10)
        flag = 1;
        kx = HC_allPairshapes{SZ10(j),2}; 
        % already normalized
        %kx = kx-min(kx);
        %kx = kx/sum(kx);
        
        for l = 1:length(idss_hc)
            ky = HC_allPairshapes{HC10(l),2}; 
            %ky = ky-min(ky);
            %ky = ky/sum(ky);
            %dEMD(kx,ky)
            if(corr(kx,ky,'Type','Spearman')> 0.7)
            fprintf(" EMD between%d and idss %d is %f\n",j,l,corr(kx,ky,'Type','Spearman'));        
            %if(dEMD(kx,ky) <= emd)
            %fprintf(" EMD between%d and idss %d is %f\n",j,l,dEMD(kx,ky));    
            flag=0;
            break;
            end
        end
        if(flag~=0)
            K_simscores(end+1) = HCP(j); 
            idss_hc(end+1)        = HC10(j);
            pals(end+1)        = j;
        end
        
end


%%
figure()
for i =1:length(HC10)
    subplot(6,8,i)
    plot(SZ_allPairshapes{(SZ10(i)),1},'r','LineWidth',2);
    hold on
    plot(HC_allPairshapes{(HC10(i)),1},'b','LineWidth',2);
    plot(zeros(length(SZ_allPairshapes{(SZ10(i)),1}),1),'k','LineWidth',2)
    hold on
    set(gca,'fontsize', 12);
    box off
    axis tight
    axis off 
end

%% Plot Pot 
figure()
for i = 6:10%length(SZ10)
    subplot(1,5,i-5)
    plot(SZ_allPairshapes{SZ10(i),1},'r','LineWidth',2);
    hold on
    plot(HC_allPairshapes{HC10(i),1},'b','LineWidth',2);
    %maxl = max(length(tempSZ));
    plot(zeros(length(SZ_allPairshapes{SZ10(i),1}),1),'k','LineWidth',2)
    %num2str(HC_allPairshapes{HC10(i),4})
    %title(['S: ' num2str(SZ_allPairshapes{SZ10(i),4}) ' H: ' num2str(HC_allPairshapes{HC10(i),4})])
    set(gca,'fontsize', 12);
    box off
    axis tight
    axis off 
end

% Top 12 diversified 
figure()
for i = 1:length(idss)%HC10)
    subplot(3,4,i)
    plot(HC_allPairshapes{(idss_hc(i)),1},'b','LineWidth',2);
    hold on
     plot(SZ_allPairshapes{(idss(i)),1},'r','LineWidth',2);
    %maxl = max(length(tempSZ));
    plot(zeros(length(HC_allPairshapes{idss_hc(i),1}),1),'k','LineWidth',2)
    %title(['H' num2str(HC_allPairshapes{idss_hc(i),4}) ' S' num2str(SZ_allPairshapes{idss(i),4})])
    set(gca,'fontsize', 12);
    box off
    axis tight
    axis off 
end
legend('HC','SZ')
%% % ******** Probably not required if it is working correctly **********  
clear
clc
dirr = '/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/EMD';
emd = 3;
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/probabilityDensity.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/subWISEPeaks.mat
% Impose diversity
%count = 1;
for h = 1:size(peaksIJ,2)
  load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/shapeletsSubWISE/',['sub_',num2str(h,'%03.f')]));
  kk = peaksIJ{h};
  sims_ii = pd{h};
  sims_i  = sims_ii(kk);      % get p(i) for those peaks 
  %shapes = zeros(length(kk),50);
  [vals,inds] = sort(sims_i,'descend');
  K_simscores = []; 
  idss        = [];
  K_simscores(1) = vals(1); 
  idss(1)        = kk(inds(1));
  % diversity ------------------------
    for j = 2:length(inds)
        flag = 1;
        kx = allshapes_s{kk(inds(j))}; 
        kx = kx-min(kx);
        kx = kx/sum(kx);
        
        for l = 1:length(idss)
            ky = allshapes_s{idss(l)}; 
            ky = ky-min(ky);
            ky = ky/sum(ky);
            %dEMD(kx,ky)
            %if(corr(kx,ky,'Type','Spearman')> 0.5)
            if(dEMD(kx,ky) <= emd)
            fprintf(" EMD between%d and idss %d is %f\n",j,l,dEMD(kx,ky));    
            flag=0;
            break;
            end
        end
        if(flag~=0)
            K_simscores(end+1) = vals(j); 
            idss(end+1)        = kk(inds(j));
        end
    end
    % -------------------------------------
    shapesIdx{h} = idss;
end
%% Random plotting. TOp 10 before global tSNE 

[~,HC10] = sort(p_HC,'descend');
[~,SZ10] = sort(p_SZ,'descend');


figure()
k =1;
for i = 1:length(SZ10(1:10))
    subplot(2,5,k)
    k = k+1;
    SZ_t10Sub_Pair(i,1) = SZ_allPairshapes{SZ10(i),3};
    SZ_t10Sub_Pair(i,2) = SZ_allPairshapes{SZ10(i),4};
    HC_t10Sub_Pair(i,1) = HC_allPairshapes{HC10(i),3};
    HC_t10Sub_Pair(i,2) = HC_allPairshapes{HC10(i),4};
    
    plot(SZ_allPairshapes{SZ10(i),1},'r','LineWidth',2);
    hold on
    plot(HC_allPairshapes{HC10(i),1},'b','LineWidth',2);
    %maxl = max(length(tempSZ));
    plot(zeros(length(SZ_allPairshapes{SZ10(i),1}),1),'k','LineWidth',2)
    %num2str(HC_allPairshapes{HC10(i),4})
    %title(['S: ' num2str(SZ_allPairshapes{SZ10(i),4}) ' H: ' num2str(HC_allPairshapes{HC10(i),4})])
    set(gca,'fontsize', 12);
    box off
    axis tight
    axis off 
end
%% Using Subject and Pair info get the full signal
clear
clc
%load SZ_top10Sub_Pair.mat
load HC_top10Sub_Pair.mat
SZ_t10Sub_Pair = HC_t10Sub_Pair;
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/data';      % server 
load (fullfile(dataDir,'dFNCs.mat'));
load (fullfile(dataDir,'candidates_22to50.mat'));
figure()
k=1;
for i=1:size(SZ_t10Sub_Pair)
    tempSig = squeeze(rFdyn(SZ_t10Sub_Pair(i,1),:,SZ_t10Sub_Pair(i,2)));
   % load(fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/SZ_pairwise_shapelet', ['pair_' num2str(SZ_t10Sub_Pair(i,2),'%04.f')]))
    load(fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/HC_pairwise_shapelet', ['pair_' num2str(SZ_t10Sub_Pair(i,2),'%04.f')]))
    for j =1:size(shapes,2)
        if(shapes(j).isLeader==1 && shapes(j).whichsubject==SZ_t10Sub_Pair(i,1))
            occur = shapes(j).occurences;
            occr{i} = shapes(j).occurences;
            candsize(i) = length(shapes(j).real_length);
            exloc = shapes(j).locExactshapelet;
            break;
        end
    end
    cands = candidates{candsize(i)-22+1};
% %     subplot(10,2,k)
% %     plot(tempSig(cands{exloc}))
% %     set(gca,'fontsize', 12);
% %     box off
% %     axis tight
% %     axis off 
%     k=k+1;
    subplot(10,1,k)
    plot(tempSig,'k')
    hold on
    for p = 1:length(occur)
        pl = nan(1,136);
        pl(cands{occur(p)}) = tempSig(cands{occur(p)});
        plot(pl,'r','LineWidth',2)
        %plot(pl,'r')
        set(gca,'fontsize', 12);
    box off
    axis tight
    axis off 
    end
    k=k+1;
end

%%