clear
clc
cd('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load SZ_top10Sub_Pair.mat
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/data';      % server 
load (fullfile(dataDir,'dFNCs.mat'));
load (fullfile(dataDir,'SZsubjectIdx.mat'));
load allPairshapes_SZ.mat
load allPairshapesPD_SZ.mat
[~,SZ10] = sort(p_SZ,'descend');
figDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation/SZConnPropagation';
% Do the central convolution: zero padded, keep the data point in the middle and take half left and half right. 
% Will return identical length like the original signal    
% conv(x,y,'same')
%% SZ
tpWISE_edge_counts = zeros(1081,136); % for Domain graph
tpWISE_edge_weights = zeros(1081,136); % for Domain graph
allCSZ = zeros(length(SZIdx),1081,136);
SZ_dominant = SZ_allPairshapes{SZ10(1),1};
SZ_dominant = flip(SZ_dominant); % flipping so that while conv function flips it will get to the original shape 
pairs_across_subs = [];
msf = 0;
sz_paircount = zeros(1081,1);
for sz = 1:length(SZIdx)
sz
SZSub = SZIdx(sz);%SUBs(i);%SZ_t10Sub_Pair(1,1);
%SZSub = SZ_t10Sub_Pair(1,1);
SZ_subdata = squeeze(rFdyn(SZSub,:,:));
SZ_subdata = SZ_subdata';
cSZ_subdata = conv2(SZ_subdata,SZ_dominant,'same');

% figure()
% H2 =imagesc(cSZ_subdata);
% colorbar
% title('original')
% alphaMask2 = cSZ_subdata ~= 0;
% set(H2, 'AlphaData', alphaMask2)
%meanCONVSZ(sz) = mean(cSZ_subdata(:));
% Threshodling
th = 0.8;
CSZ = abs(cSZ_subdata)>=th;
allCSZ(sz,:,:) = CSZ;

%CSZ = cSZ_subdata>=1;
thcSZ = cSZ_subdata.*CSZ;

% figure()
% H =imagesc(thcSZ);
% colorbar
% title('thresholded')
% alphaMask = thcSZ ~= 0;
% set( H, 'AlphaData', alphaMask )
% Filtering 
fthcSZ = [];
fCSZ = [];
sPairs = [];
firstone = [];
k =1;
for i = 1:size(thcSZ,1)
    if(nnz(CSZ(i,:))>1)% if the appearence is >1 across timepoints 
    fthcSZ(k,:) = thcSZ(i,:);
    fCSZ(k,:) = CSZ(i,:);
    temp1 = find(fCSZ(k,:)==1);
    firstone(k) = min(temp1);
    sPairs(k) = i;
    k = k+1;
    end 
end

%{
 figure()
 H1 =imagesc(fthcSZ);
 colorbar
 title('Filtered')
 alphaMask1 = fthcSZ ~= 0;
 set(H1, 'AlphaData', alphaMask1)
%}
% which comes earlier/sorting
[appeared_ts,sfirstoneIDX] = sort(firstone);
sfthcSZ = fthcSZ(sfirstoneIDX,:);
SortedPairs_Links_sz{sz} = sfthcSZ;
msf = max(max(abs(sfthcSZ(:))),msf);
SortedPairs = sPairs(sfirstoneIDX);
SortedPairs_for_sz{sz} = SortedPairs;
sz_paircount(SortedPairs) = sz_paircount(SortedPairs)+1;
% Sorted pairs according to their appearence [Pairs Timestamps]
sPairs_Tiemstamps = [SortedPairs',appeared_ts'];
sPairs_TS_sz{sz} = sPairs_Tiemstamps;
pairs_across_subs = union(pairs_across_subs,SortedPairs);
% Timepointwise counting for domain graph
for t = 1:size(sfthcSZ,2)
pairs_tidx = find(sfthcSZ(:,t)~=0);
pairs_t = SortedPairs(pairs_tidx);
    for p =1:length(pairs_t)
        tpWISE_edge_counts(pairs_t(p),t) = tpWISE_edge_counts(pairs_t(p),t)+1; 
        tpWISE_edge_weights(pairs_t(p),t) = tpWISE_edge_weights(pairs_t(p),t)+SZ_subdata(pairs_t(p),t);
    end
end

%{
figure()
Hf =imagesc(sfthcSZ);
caxis([-3.5840 3.5840])
colorbar
title(['SZ ' num2str(SZSub)])
alphaMaskf = sfthcSZ ~= 0;
set(Hf, 'AlphaData', alphaMaskf)
fName = ['SZLinkstrength',num2str(sz),'.png'];
saveas(Hf,fullfile(figDir,fName));
close
%saveas('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation/SZ',gcf)
%}
end
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation/StrongLinks_SZ1.mat','SortedPairs_Links_sz');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation/sLinkedPairs_SZ1.mat','SortedPairs_for_sz');
%% HC Initialize
clear
clc
cd('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load HC_top10Sub_Pair.mat
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/data';      % server 
load (fullfile(dataDir,'dFNCs.mat'));
load (fullfile(dataDir,'HCsubjectIdx.mat'));
load allPairshapes_HC.mat
load allPairshapesPD_HC.mat
[~,HC10] = sort(p_HC,'descend');
figDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation/HCConnPropagation';
%% HC
tpWISE_edge_counts = zeros(1081,136); % for Domain graph
tpWISE_edge_weights = zeros(1081,136); % for Domain graph
allCHC = zeros(length(CtIdx),1081,136);
HC_dominant = HC_allPairshapes{HC10(1),1};
HC_dominant = flip(HC_dominant); % flipping so that while conv function flips it will get to the original shape 
pairs_across_subs = [];
hc_paircount = zeros(1081,1);
msf = 0;
for hc =1:size(CtIdx,1)
hc
HCSub = CtIdx(hc);%HC_t10Sub_Pair(1,1);
%HCSub = HC_t10Sub_Pair(1,1);
HC_subdata = squeeze(rFdyn(HCSub,:,:));
HC_subdata = HC_subdata';
cHC_subdata = conv2(HC_subdata,HC_dominant,'same');

% figure()
% H2 =imagesc(cHC_subdata);
% colorbar
% alphaMask2 = cHC_subdata ~= 0;
% set(H2, 'AlphaData', alphaMask2)
meanCONVHC(hc) = mean(cHC_subdata(:));
% Threshodling
%CHC = cHC_subdata>=0.6;
th = 0.8;
CHC = abs(cHC_subdata) >= th;
allCHC(hc,:,:) = CHC;
thcHC = cHC_subdata.*CHC;

% figure()
% H =imagesc(thcHC);
% colorbar
% alphaMask = thcHC ~= 0;
% set( H, 'AlphaData', alphaMask )

% Filtering 
fthcHC = [];
fCHC = [];
sPairs = [];
peak1 = [];
k = 1;
%maxnumberofpeaks =0;
%allPeaks = zeros(577,4);
for i = 1:size(thcHC,1)
    if(nnz(CHC(i,:))>1)
    %i
    fthcHC(k,:) = thcHC(i,:);
    fCHC(k,:) = CHC(i,:);
    temp1 = find(fCHC(k,:)==1);
    stemp1 = sort(temp1);
    peak1(k) = stemp1(1);
    sPairs(k) = i;
    k = k+1;
    end
end

[appeared_ts,sfirstoneIDX] = sort(peak1);
sfthcHC = fthcHC(sfirstoneIDX,:);
SortedPairs_Links_hc{hc} = sfthcHC;
msf = max(max(abs(sfthcHC(:))),msf);
SortedPairs = sPairs(sfirstoneIDX);
SortedPairs_for_hc{hc} = SortedPairs;

% Sorted pairs according to their appearence [Pairs Timestamps]
sPairs_Tiemstamps = [SortedPairs',appeared_ts'];
sPairs_TS_hc{hc} = sPairs_Tiemstamps;
hc_paircount(SortedPairs) = hc_paircount(SortedPairs)+1;
pairs_across_subs = union(pairs_across_subs,SortedPairs);
% Timepointwise counting for domain graph
for t = 1:size(sfthcHC,2)
pairs_tidx = find(sfthcHC(:,t)~=0);
pairs_t = SortedPairs(pairs_tidx);
    for p =1:length(pairs_t)
        tpWISE_edge_counts(pairs_t(p),t) = tpWISE_edge_counts(pairs_t(p),t)+1; 
        tpWISE_edge_weights(pairs_t(p),t) = tpWISE_edge_weights(pairs_t(p),t)+HC_subdata(pairs_t(p),t);
    end
end

%{
figure()
Hf =imagesc(sfthcHC);
caxis([-3.5840 3.5840])
colorbar
title(['HC ' num2str(HCSub)])
alphaMaskf = sfthcHC ~= 0;
set(Hf, 'AlphaData', alphaMaskf)
fName = ['HCLinkstrength',num2str(hc),'.png'];
saveas(Hf,fullfile(figDir,fName));
close
%}
end
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation/StrongLinks_HC1.mat','SortedPairs_Links_hc');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation/sLinkedPairs_HC1.mat','SortedPairs_for_hc');
%% Get the common link propagation order for both groups 

% HC
clear
clc
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation'))
load sLinkedPairs_HC.mat
%load StrongLinks_HC.mat

% Find the counts of each pair across the subjects
HC_Pair_occureces = zeros(1081,1);  
for i = 1:size(SortedPairs_for_hc,2)
    sub_i_pairs = SortedPairs_for_hc{1,i};
    for j=1:size(sub_i_pairs,2)
        HC_Pair_occureces(sub_i_pairs(j),1) = HC_Pair_occureces(sub_i_pairs(j),1)+ 1;
    end
end
[sHC_Pairs(:,2),sHC_Pairs(:,1)] = sort(HC_Pair_occureces,'descend');

% SZ
clear
clc
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation'))
load sLinkedPairs_SZ.mat
%load StrongLinks_SZ.mat
SZ_Pair_occureces = zeros(1081,1);  
for i = 1:size(SortedPairs_for_sz,2)
    sub_i_pairs = SortedPairs_for_sz{1,i};
    for j=1:size(sub_i_pairs,2)
       SZ_Pair_occureces(sub_i_pairs(j),1) = SZ_Pair_occureces(sub_i_pairs(j),1)+ 1;
    end
end
[sSZ_Pairs(:,2),sSZ_Pairs(:,1)] = sort(SZ_Pair_occureces,'descend');

% Plotting SZ & HC occurences 
load sHCpairs_appeareces.mat
load sSZpairs_appeareces.mat

%figure()
%plot(sSZ_Pairs(:,2),'r','LineWidth',2)
%hold on 
%plot(sHC_Pairs(:,2),'b','LineWidth',2)
%% Get the graph ready for each subject for both group 
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/data/Paircoordinates.mat
load sLinkedPairs_HC.mat

for i = 1:size(SortedPairs_for_hc,2)
tshownHCpairs = SortedPairs_for_hc{i};     
HCE = coord(tshownHCpairs,:);
allHCE{i} = HCE;
%HCE(:,3) = tshownHCpairs(:,2);
end


%tshownHCpairs = sHC_Pairs(1:100,:);
%tshownSZpairs = sSZ_Pairs(1:100,:);
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/data/Paircoordinates.mat
load sLinkedPairs_SZ.mat
for i = 1:size(SortedPairs_for_sz,2)
tshownSZpairs = SortedPairs_for_sz{i};     
SZE = coord(tshownSZpairs,:);
allSZE{i} = SZE;
%HCE(:,3) = tshownHCpairs(:,2);
end
%SZE = coord(tshownSZpairs(:,1),:);
%SZE(:,3) = tshownSZpairs(:,2);
%% How the link graph changes 
clear
clc
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/data';     
load (fullfile(dataDir,'SZsubjectIdx.mat'));
load allCSZ.mat

%% SZ
allSZPropagation = zeros(length(SZIdx),136);
for sz = 1:length(SZIdx)
CSZ = squeeze(allCSZ(sz,:,:));
sCSZ(1) = sum(CSZ(:,1));
uSZpairs = find(CSZ(:,1)==1);
for i=2:size(CSZ,2)
    iSZpairs = find(CSZ(:,i)==1);
    riSZpairs = setdiff(iSZpairs,uSZpairs);
    sCSZ(i) = sum(CSZ(riSZpairs,i));
    uSZpairs = union(uSZpairs,riSZpairs);
end
allSZPropagation(sz,:) = sCSZ;
end

%% HC
clear
clc
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/data';     
load (fullfile(dataDir,'HCsubjectIdx.mat'));
load allCHC.mat
allHCPropagation = zeros(length(CtIdx),136);
for hc = 1:length(CtIdx)
CHC = squeeze(allCHC(hc,:,:));
sCHC(1) = sum(CHC(:,1));
uHCpairs = find(CHC(:,1)==1);
for i = 2:size(CHC,2)
    iHCpairs = find(CHC(:,i)==1);
    riHCpairs = setdiff(iHCpairs,uHCpairs);
    sCHC(i) = sum(CHC(riHCpairs,i));
    uHCpairs = union(uHCpairs,riHCpairs);
end
allHCPropagation(hc,:) = sCHC;
end

% figure()
% plot(sCHC,'b','LineWidth',2);
% hold on
% plot(sCSZ,'r','LineWidth',2);

%% end of link graph changes 

for i =1:size(SZ_allPairshapes,1)
    subj(i)  = SZ_allPairshapes{i,3};
end
%% Domain graph creation for both SZ and HC 
% edges types: inter domians, within domain  
% The weight of edges between two domains corresponds to teh strength of the overall connectivity produced by all the the pairs within 
% these two domians 
% (number of active conenctions * their weights)/all possible connections between two given domains 

% HC
clear
clc
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation'));
load HCedgecounts.mat
load HCedgeweights.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/Paircoordinates.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/Comp2DomMapping.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/IC2DOM_int.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/DomainEdges.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/DOmainEdges_string.mat
e_src  = edges(:,1);
e_dest = edges(:,2);
max_edge_weight = 0;
load edgeWISEpossibleconnections.mat
for i = 1:size(tpWISE_edge_counts,2)
    i
    pairs_i = tpWISE_edge_counts(:,i);
    [midx] = find(pairs_i>=15);
    mcidx{i} = midx; % number of hyper active edges  
    e_weight = zeros(28,1);
    e_activeconnections = zeros(28,1);
    flag = 0;
    for j=1:length(midx)
        flag =1;
        comp1 = coord(midx(j),1); 
        comp2 = coord(midx(j),2); 
        % ........get the edge id......
        for h =1:length(e_src)
            if((e_src(h) == IC2DOM(1,comp1)&& e_dest(h)==IC2DOM(1,comp2))||(e_src(h) == IC2DOM(1,comp2)&& e_dest(h)==IC2DOM(1,comp1)))
                id = h; % which edge it is?
                break;
            end
        end
        %.................................
        e_weight(id)= e_weight(id)+ (tpWISE_edge_weights(midx(j),i)/tpWISE_edge_counts(midx(j),i)); % Multiple pair can correspond to same edge 
        %.. because it is domainwise. This is average across the subjects
        e_activeconnections(id) = e_activeconnections(id)+1; % count the number of pairs/connections weighting that domain-domain edge
        
        %me_weight = max(abs(e_weight(:)));
    end
    for e = 1:size(e_weight,1)
        if(e_activeconnections(e)~=0)
        % compute the mean activity weight by an inter domian edge
        % weighted sum of all active connections/ possible connections between the domians  
        %e_weight(e) = e_weight(e)/e_possible_connections(e); 
        e_weight(e) = e_weight(e)/e_activeconnections(e); 
        end
    end
    me_weight = max(abs(e_weight(:)));
    if(flag==1)
        e_connections{i} = e_activeconnections;
        e_weights{i} = e_weight;
      
    else
      e_weights{i} = 0;
      e_connections{i} = 0;
    end
    max_edge_weight = max(max_edge_weight,me_weight);
end
max_edge_weight
%% SZ
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/data/Paircoordinates.mat
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis/LinkPropagation'));
load SZedgecounts.mat
load SZedgeweights.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/Comp2DomMapping.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/IC2DOM_int.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/DomainEdges.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/DOmainEdges_string.mat
load edgeWISEpossibleconnections.mat
e_src  = edges(:,1);
e_dest = edges(:,2);
max_edge_weight = 0;
for i = 1:size(tpWISE_edge_counts,2)
    i
    pairs_i = tpWISE_edge_counts(:,i);
    [midx] = find(pairs_i>=15);
    mcidx{i} = midx;
    e_weight = zeros(28,1);
    e_activeconnections = zeros(28,1);
    flag = 0;
    length(midx);
    for j=1:length(midx)
        
        flag =1;
        comp1 = coord(midx(j),1); 
        comp2 = coord(midx(j),2); 
        % ........geth the edge id......
        for h =1:length(e_src)
            if((e_src(h) == IC2DOM(1,comp1)&& e_dest(h)==IC2DOM(1,comp2))||(e_src(h) == IC2DOM(1,comp2)&& e_dest(h)==IC2DOM(1,comp1)))
                id = h;
                break;
            end
        end
        %.................................
        e_weight(h)= e_weight(h)+ (tpWISE_edge_weights(midx(j),i)/tpWISE_edge_counts(midx(j),i)); % Multiple pair can correspond to same edge 
        %.. because it is domainwise
        e_activeconnections(id) = e_activeconnections(id)+1; % count the number of pairs/connection weighting that domain-domain edge
         
    end
    for e = 1:size(e_weight,1)
        if(e_activeconnections(e)~=0)
        % compute the mean activity weight by an inter domian edge
        % weighted sum of all active connections/ possible connections between the domians  
        %e_weight(e) = e_weight(e)/e_possible_connections(e); 
        e_weight(e) = e_weight(e)/e_activeconnections(e); 
        end 
    end
    me_weight = max(abs(e_weight(:)));
    if(flag==1)
      e_weights{i} = e_weight;
      e_connections{i} = e_activeconnections;
    else
      e_weights{i} = 0;
      e_connections{i} = 0;
    end
    max_edge_weight = max(max_edge_weight,me_weight);
end
max_edge_weight
%% supporting data creation for the domain graph edges
%{
for d = 1:size(Comp2Dom,2)
        if(Comp2Dom(d)==1)
        IC2DOM(d) = "SC";
        elseif(Comp2Dom(d)==2)
        IC2DOM(d) = "AUD";
        elseif(Comp2Dom(d)==3)
        IC2DOM(d) = "VIS";
        elseif(Comp2Dom(d)==4)
        IC2DOM(d) = "SM";
        elseif(Comp2Dom(d)==5)
        IC2DOM(d) = "CC";
        elseif(Comp2Dom(d)==6)
        IC2DOM(d) = "DM";
        elseif(Comp2Dom(d)==7) 
        IC2DOM(d) = "CB";
        end
end

for d = 1:size(Comp2Dom,2)
        if(Comp2Dom(d)=="SC")
        IC2DOM(d) = 1;
        elseif(Comp2Dom(d)=="AUD")
        IC2DOM(d) = 2;
        elseif(Comp2Dom(d)=="VIS")
        IC2DOM(d) = 3;
        elseif(Comp2Dom(d)=="SM")
        IC2DOM(d) = 4;
        elseif(Comp2Dom(d)=="CC")
        IC2DOM(d) = 5;
        elseif(Comp2Dom(d)=="DM")
        IC2DOM(d) = 6;
        elseif(Comp2Dom(d)=="CB") 
        IC2DOM(d) = 7;
        end
end
for d = 1:size(edges,2)
    for f = 1:size(edges,1)
        if(edges(f,d)==1)
        Edges(f,d) = "SC";
        elseif(edges(f,d)==2)
        Edges(f,d) = "AUD";
        elseif(edges(f,d)==3)
        Edges(f,d) = "VIS";
        elseif(edges(f,d)==4)
        Edges(f,d) = "SM";
        elseif(edges(f,d)==5)
        Edges(f,d) = "CC";
        elseif(edges(f,d)==6)
        Edges(f,d) = "DM";
        elseif(edges(f,d)==7) 
        Edges(f,d) = "CB";
        end
    end
end

v = 1:7;
clear edges
edges = nchoosek(v,2);
k=1;
for l = size(edges,1)+1:size(edges,1)+7
edges(l,:) = [k,k];
k=k+1;
end

%}
%% Groupwise average graph based on time decay, Note: across the subjects it brings up all the pairs for both SZ and HC

clear
clc
load sPairs_Timestamps_HC.mat
load sPairs_Timestamps_SZ.mat
SZ_edge_weights = zeros(1081,1);
HC_edge_weights = zeros(1081,1);
SZ_edge_counts = zeros(1081,1);
HC_edge_counts = zeros(1081,1);
%for i =1:1081
% HC
    for j=1:size(sPairs_TS_hc,2)
        connections = sPairs_TS_hc{j};
        for k =1:size(connections,1)
            HC_edge_counts(connections(k,1),1) = HC_edge_counts(connections(k,1),1)+1; % Number of subjects
            HC_edge_weights(connections(k,1),1) = HC_edge_weights(connections(k,1),1)+ 1/connections(k,2);
        end
    end
%end
HC_edge_weights = HC_edge_weights./163;

% SZ

    for j=1:size(sPairs_TS_sz,2)
        connections1 = sPairs_TS_sz{j};
        for k =1:size(connections1,1)
            SZ_edge_counts(connections1(k,1),1) = SZ_edge_counts(connections1(k,1),1)+1; % Number of subjects
            SZ_edge_weights(connections1(k,1),1) = SZ_edge_weights(connections1(k,1),1)+ 1/connections1(k,2); % Tie decay equation
        end
    end

SZ_edge_weights = SZ_edge_weights./151;   

% thresholding

clear
clc
load HC_avggraph_edges.mat
load SZ_avggraph_edges.mat

%tHC_edge_weights = zeros(1081,1);
%tSZ_edge_weights = zeros(1081,1);
th = (mean(HC_edge_weights)+mean(SZ_edge_weights))/2;
tHC_edge_weights = HC_edge_weights>th;
tSZ_edge_weights = SZ_edge_weights>th;
%shc = find(HC_edge_weights>=th);
%ssz = find(SZ_edge_weights>=th)

%% 
clear
clc


%%







