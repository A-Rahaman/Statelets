%% post hoc analysis after subjectwise tSNE
% to get the dominant shapes in a given subject
clear
clc
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/subWISEtSNE';
load /Users/mrahaman1/Documents/Statelet_V2/fResults/probabilityDensity.mat
%% Map the coordinates a 2D matrix and weight each of the cell with probability density
pMats = cell(314,1);
for i = 1:314
    i
    Ps = pd{1,i};
    Ysub = load(fullfile(dataDir,['sub_',num2str(i),'.txt']));
    X = Ysub(:,1);
    Y = Ysub(:,2);
    pMat = zeros(1+ceil((max(X)-min(X))),1+ceil((max(Y)-min(Y))));
    %pMatfreq = zeros(1+ceil((max(X)-min(X))),1+ceil((max(Y)-min(Y))));
    for j = 1:length(X)
        x = 1+round(X(j)+(-1)*min(X));
        y = 1+round(Y(j)+(-1)*min(Y));
        pMat(x,y) = pMat(x,y)+Ps(j);
        %pMatfreq(x,y)
    end
    pMats{i} = pMat;
end
save('/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/pMats.mat','pMats');
%% Find the peaks
%[pks,locs] = findpeaks(bpMat,'MinPeakDistance',6);
%[pks,locs] = findpeaks(bpMat);
%[cent, cm] = FastPeakFind(pMat);
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/pMats.mat
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/anlysisScripts/FastPeakFind'));
peaklocs = cell(314,1);

for i = 1:length(pMats)
    i
pMat = pMats{i,1};
bpMat = imgaussfilt(pMat,4.5);
%bpMat = bpMat/max(bpMat(:));
[cent, cm] = FastPeakFind(bpMat);
%ii = imregionalmax(bpMat);
%[xloc,yloc,~] = find(ii);        % coordinates
[xloc,yloc] = find(cm==1);
locs = [xloc yloc];
peaklocs{i} = locs;
end
save('/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/Peaklocs.mat','peaklocs');
% imagesc(bpMat)
% hold on
% for fn=1:length(x)
%     text(y(fn),x(fn),'*','FontSize',25)
% end
% hold off
%mesh(0.01*cm+bpMat)
%mesh(ii+bpMat)
%% Take 10 closest neigbour and get the max of them
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/Peaklocs.mat
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/subWISEtSNE';
load /Users/mrahaman1/Documents/Statelet_V2/fResults/probabilityDensity.mat
for i = 1:length(peaklocs)
    i
    coord = peaklocs{i};
    Ysub = load(fullfile(dataDir,['sub_',num2str(i),'.txt']));
    P = pd{i};
    
    peaks_i_j = []; % Crucial   
    for j=1:length(coord)
        cx = coord(j,:);
        cx(1) = cx(1)+min(Ysub(:,1)); % map to tSNE plot coordiante x
        cx(2) = cx(2)+min(Ysub(:,2)); % map to tSNE plot coordiante y
        
        % Find the closest match ---------------------------------
        d = zeros(length(Ysub),1); % Crucial
        for k=1:length(Ysub)
            cy = Ysub(k,:);
            d(k) = norm(cx-cy);
        end
        [~,clId] = min(d);
        clCoord = Ysub(clId,:);  
        realPeaks(j) = clId;
        
        % --------------------------------------------------------
        % distance from that real point in tSNE to get 10 nearest
        dd = zeros(length(Ysub),1); % Crucial
        for k=1:length(Ysub)    
            cy = Ysub(k,:);
            dd(k) = norm(clCoord-cy);
        end
        [sortedDD,idx] = sort(dd);
        I =  idx(1:10);           
        PI = P(I); 
        % ---------------------------------------------------------
        % Find one with the optimal probability within that 10 
        [~,mPI] = max(PI);
        peak_j = I(mPI);
        peaks_i_j(j) = peak_j; 
        % ---------------------------------------------------------
    end
    peaksIJ{i} = peaks_i_j; 
    % Filter some of the peaks which are very close to each other. keep maxp always
%     peaks_i = peaksIJ{i};
%     for l = 1:length(peaks_i)
%         
%     end
    % --------------------------------------------------------------
end
save('/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/subWISEPeaks.mat','peaksIJ');
% column 1 = real_length, col 2: extrapolated, col 3: pair, col 4: subject class
%{
ss = peaksIJ{1,1};
for i=1:length(ss)
pll(i,:)= allshapes_s{ss(i),2};
end
%}
%% filter out the redundancy 

% ******** Probably not required if it is working correctly **********  
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

save(fullfile(dirr,['subWISEPeaks_',num2str(emd),'.mat']),'shapesIdx');
%}
%% Grab the real length shapelets per subject.
clear
clc
dirr = '/Users/mrahaman1/Documents/Statelet_V2/fResults/subWISE_real_shapes';
for k =1:314    
    k
    load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/shapeletsPairWISE',['sub_',num2str(k,'%03.f')]));
    c=1;
    clear subwise_real_shapes
    for l =1:size(pairWISEshapes,1)
        p = 1;
        while(p<= size(pairWISEshapes,2)&& ~isempty(pairWISEshapes{l,p}) )
            subwise_real_shapes{c} = pairWISEshapes{l,p};
            c = c+1;
            p=p+1;
        end
    end
    save(fullfile(dirr,['subWISE_real_shapes_',num2str(k),'.mat']),'subwise_real_shapes');
end
%% Getting the most dominant shapes 
clear
clc
dirr = '/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD';
%emd = 3;
%load(fullfile(dirr,['subWISEPeaks_',num2str(emd),'.mat']));
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/EMD/subWISEPeaks.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/subWISEPeaks.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/subWISEPeaks_CORRfiltered.mat
%peaksIJ = shapesIdx;  
%subWISEshapes = cell(314,1);
% load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/data/EMDresults/allshapesof_a_sub_EMD',['sub_',num2str(1,'%03.f')]));
% tt = peaksIJ{1};
% allSUBshapes(1,50) = allshapes_s{tt(1)};
count = 1;
for h =1:size(peaksIJ,2)
  %load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/shapeletsSubWISE',['sub_',num2str(h,'%03.f')]));
  load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/subWISEshapes',['sub_',num2str(h,'%03.f'),'.mat']));
  kk = peaksIJ{h};
  %shapes = zeros(length(kk),50);
  for pp = 1:length(kk)
      %shapes(pp,:) = allshapes_s{kk(pp),2};
      %allSUBshapes(count,:) = allshapes_s{kk(pp),2};
      %{
      allSUBshapes{count,1} = allshapes_s{kk(pp),1}; % for grabing the real length 
      allSUBshapes{count,2} = allshapes_s{kk(pp),2}; % extrapolated
      allSUBshapes{count,3} = allshapes_s{kk(pp),3}; % pair
      allSUBshapes{count,4} = allshapes_s{kk(pp),4}; % subejct class
      count = count+1;
      %}
      mostFreq{count,1} = allSUBshapes{kk(pp),1}; % for grabing the real length 
      mostFreq{count,2} = allSUBshapes{kk(pp),2}; % extrapolated
      mostFreq{count,3} = allSUBshapes{kk(pp),3}; % pair
      mostFreq{count,4} = allSUBshapes{kk(pp),4}; % subejct class
      mshapes{count} = allSUBshapes{kk(pp),1};
      count = count+1;
      
  end
  %subWISEshapes{h} = shapes;
end
%save(fullfile(dirr,['subWISEshapes_',num2str(emd),'.mat']),'subWISEshapes');
save(fullfile(dirr,['allSUBshapes' '.mat']),'allSUBshapes');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/subWISEshapes.mat','subWISEshapes');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes.mat','allSUBshapes');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/subWISEshapes.mat','subWISEshapes');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/allSUBshapes.mat','allSUBshapes');

%% get EMD for all the shapes from all the subjects 
clear
clc
for i=1:length(allSUBshapes)
    allSUBshapesEXT(i,:) = allSUBshapes{i,2}; 
end
%dirr = '/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/EMD';
%emd = 3;
%load(fullfile(dirr,['allSUBshapes_',num2str(emd),'.mat']));
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/EMD/allSUBshapes.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/HCshapes.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/allSUBshapes.mat
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/shapeletAlgorithmScripts'));
%allSUBshapes = HCshapes; % while testing groupwise;
for j = 1:size(allSUBshapesEXT,1)
  x = allSUBshapesEXT(j,:);  
  
  for k = 1:size(allSUBshapesEXT,1)
    %if j~=k
    y = allSUBshapesEXT(k,:);
    x = x-min(x);
    x = x/sum(x);
    y = y-min(y);
    y = y/sum(y);
    simscore(j,k) = dEMD(x,y);
       
    %end 
  end
end
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBEMD.mat','simscore');
%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/allSUBEMD.mat','simscore');
%save(fullfile(dirr,['allSUBEMD_',num2str(emd),'.mat']),'simscore');
%% shape to subejct's id mapping 
clear
clc
dirr = '/Users/mrahaman1/Documents/Statelet_V2/fResults';
%emd = 3;
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/subWISEshapes.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/subWISEshapes.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/EMD/subWISEshapes.mat
load(fullfile(dirr,['subWISEshapes' '.mat']));
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/SZsubjectIdx.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/HCsubjectIdx.mat
shape_to_sub = [];
count = 1;
shape_to_sub(1:size(subWISEshapes{1},1)) = 1;
shape_to_subcl(1:size(subWISEshapes{1})) = 0; % because subject one is SZ
for i = 2:size(subWISEshapes,1)
    pp = subWISEshapes{i};
    shape_to_sub(end+1:size(shape_to_sub,2)+size(pp,1)) = i;
    if(~isempty(intersect(SZIdx,i)))
        shape_to_subcl(end+1:size(shape_to_subcl,2)+size(pp,1)) = 0; 
    else
        shape_to_subcl(end+1:size(shape_to_subcl,2)+size(pp,1)) = 1; 
    end
    
end
%

%save(fullfile(dirr,['allshapestosubmapping_',num2str(emd),'.mat']),'shape_to_sub');

%save('/Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/allshapestosubmapping.mat','shape_to_sub');
%% k-means on EMD based dominant shapes 
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/allSUBshapes.mat
data = allSUBshapes;
k = 10;
clear C IDX SUMD D  
%Worked really bad
%[IDX, C, SUMD, D]= kmeans(data, k,'Distance','sqeuclidean','MaxIter',10000,'Replicates',100);

[IDX, C, SUMD, D]= kmeans(data, k,'Distance','correlation','MaxIter',10000,'Replicates',100);
%IDXX = reshape(IDX,314,5); % Return back to the original data dimension 

load /Users/mrahaman1/Documents/Statelet_V2/fResults/allshapestosubmapping.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/afterFiltering_dominants/CORR/allshapestosubmapping.mat
Statelets = cell(k,2);
for i=1:k
    IDZ = find(IDX==i);
    subs = shape_to_sub(IDZ);
    Statelets{i,1} = unique(subs);
    Statelets{i,2} = IDZ;
end

%{
for i = 1:k
[s,sh]= find(IDXX==i);    
pairs = reps_COMM{s(1),sh(1)};
for p =1:length(s)
    pairs = union(pairs,reps_COMM{s(p),sh(p)});
    %pairs = intersect(pairs,reps_COMM{s(p),sh(p)});
end
Statelets{i,1} = s;     % Subejcts
Statelets{i,2} = sh;    % Shapelets 
Statelets{i,3} = pairs; % Pairs 
end
%}
%% testing k-means
for i =1:size(C,1)
subplot(3,3,i)
plot(C(i,:))
end
%% testing

%{
load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/allshapesof_a_sub_EMD',['sub_',num2str(i,'%03.f')]));
kk = peaksIJ{245};
for i =1:length(kk)
    %subplot(ceil(sqrt(length(kk)),sqrt(length(kk)),i));
    subplot(4,4,i);
    plot(allshapes_s{kk(i)}); 
end
% ------------ best k
pkkj = P(kk);
[~,spkkj] = sort(pkkj,'descend');
lt = kk(spkkj(1:5));
lt = kk;
for i =1:length(lt)
    subplot(4,4,i);
    plot(allshapes_s{lt(i)}); 
end
%}
%% Distances
% dists = zeros(length(realPeaks),length(realPeaks)-1);
% for ii =1:length(realPeaks)
%    cy = Ysub(realPeaks(ii),:);
%    count =1;
%    for jj =1:length(realPeaks)
%        if(ii~=jj)
%        cc = Ysub(realPeaks(jj),:);    
%        dists(ii,count) = norm(cc-cy);
%        count = count+1;
%        end
%    end
%    [~,mdist] = min(dists(ii,:));
% end

%% subject level shapes to subject level pair.
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/subWISEPeaks.mat
for k =1:314
    k
    load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/shapeletsSubWISE',['sub_',num2str(k,'%03.f')]));
    load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/shapeletsPairWISE',['sub_',num2str(k,'%03.f')]));
    shapetopair_subk = zeros(size(allshapes_s,1),1);
    c=1;
    for l =1:size(pairWISEshapes,1)
        p = 1;
        while(p<= size(pairWISEshapes,2)&& ~isempty(pairWISEshapes{l,p}) )
            
            shapetopair_subk(c,1) = l;
            c = c+1;
            p=p+1;
        end
    end
    save(fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/shape2pair_mapping_subjectWISE',['sub__',num2str(k),'.mat']),'shapetopair_subk');
    
end
%% garb the real length of the dominant shapes. without extrapolation 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/subWISEPeaks.mat

for h =1:314
    kk = peaksIJ{h};
    
end
%%
%% Mapping the shapes to the pair for all shapes 
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allshapestosubmapping.mat
shape_to_pair = zeros(size(shape_to_sub,2),1);
load /Users/mrahaman1/Documents/Statelet_V2/fResults/subWISEPeaks.mat
%for i =1:size(shape_to_sub,2)
for i =1:314
    load (fullfile('/Users/mrahaman1/Documents/Statelet_V2/fResults/shape2pair_mapping_subjectWISE',['sub__',num2str(i)]));
    shapes_sub_ids = find(shape_to_sub==i);
    shape_to_pair(shapes_sub_ids) = shapetopair_subk(peaksIJ{i});
    %peakIDs = 
end
%% Creating histogram plot for SZ - HC pairwise summary 
clear
clc
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/shapestosubjectclassMapping.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/allSUBshapes.mat
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/allshapestoPairmapping.mat
for i=1:length(allSUBshapes)
shape_to_subcl(i) = allSUBshapes{i,4};
shape_to_pair(i)  = allSUBshapes{i,3};
end
%SZshapes = find(shape_to_subcl==0); % SZ shapes
%HCshapes = find(shape_to_subcl==1); % HC shapes
SZshapes_Pairs = shape_to_pair(shape_to_subcl==0);
HCshapes_Pairs = shape_to_pair(shape_to_subcl==1);
SZ_i = zeros(1081,1);
HC_i = zeros(1081,1);
for i =1:1081
    SZ_i(i) = length(find(SZshapes_Pairs==i));
    HC_i(i) = length(find(HCshapes_Pairs==i));
end
SZ_HC = [SZ_i HC_i];
barh(SZ_HC)
axis tight
%% Sorting red with blue and sorting blue with red 
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/SZPaircount.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/HCPaircount.mat

% sort SZ
[~,sortedSZ_idx] = sort(SZ_i,'descend');
sortedSZs = SZ_i(sortedSZ_idx);        % SZ sorted hight to low counts  
sortedHCs = HC_i(sortedSZ_idx);        % HC sorted accrodingly
sortedSZ_HC = [sortedSZs sortedHCs];

% sort HC
[~,sortedHC_idx] = sort(HC_i,'descend');
sortedHCs1 = HC_i(sortedHC_idx);        % HC sorted hight to low counts  
sortedSZs1 = SZ_i(sortedHC_idx);        % SZ sorted accrodingly 
sortedHC_SZ = [sortedSZs1 sortedHCs1];  

figure()
subplot(1,2,1)
h = barh(sortedSZ_HC);
axis tight
h(1).FaceColor = 'r'; % color
h(2).FaceColor = 'b';%[0.3 0.78 0.5];% rgb
subplot(1,2,2)
h1 = barh(sortedHC_SZ);
h1(1).FaceColor = 'r'; % color
h1(2).FaceColor = 'b';%[0.3 0.78 0.5];% rgb
axis tight

%% heat map; just plot the counts
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
load /Users/mrahaman1/Documents/mrahaman/d_NBiC/Data/idxUpper.mat % Upper 1081 idx (IND1)
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = 'hot';
%cm = 'jet';
% FDR correction of tvals then keeping only the significant tvals    
FDRtvals = [];
load Peakcounts_HC.mat
load Peakcounts_SZ.mat
%% HC
HC_i = counts_HC;
datt1 = zeros(47);    
%datt1(idx) = log(HC_i+eps); 
datt1(idx) = (HC_i);  
datt21 = rot90(fliplr(datt1));
Fdatt1 = zeros(47);
%Fdatt1(idx) = log(HC_i+eps);
Fdatt1(idx) = (HC_i);

FNC1 = datt21+Fdatt1; 
%FNC = datt2+datt;
%pvals = tvalMat{db3,1}; % pval
%CLIM = [-max(abs(HC_i(:))) max(abs(HC_i(:)))];
CLIM = [0 max((Fdatt1(:)))];
T = 'HC pairwise activity';
%fName = ['t-val',num2str(db3) '.fig'];
%[F,A,C,I] = plot_FNC(rot90(fliplr(datt)), CLIM, LABEL, RSN_I,'',T,MOD,cm,pvals,0); % db3 as last argument for normal plot 
[F1,A1,C1,I1] = plot_FNC(FNC1, CLIM, LABEL, RSN_I,'',T,MOD,cm,Fdatt1,0); 
%% SZ
SZ_i = counts_SZ;
datt = zeros(47);    
%datt(idx) = log(SZ_i+eps);  % tval
datt(idx) = (SZ_i);  % tval
datt2 = rot90(fliplr(datt));
Fdatt = zeros(47);
%Fdatt(idx) = log(SZ_i+eps);
Fdatt(idx) = (SZ_i);
FNC2 = datt2+Fdatt; % Adding two matrices 47X47. One have lower traingle all incorrected values and one have upper traingle with FDR corrected.  
%FNC = datt2+datt;
%pvals = tvalMat{db3,1}; % pval
%CLIM = [-max(abs(HC_i(:))) max(abs(HC_i(:)))];
CLIM = [0 max((Fdatt(:)))];
T = 'SZ pairwise activity';
%fName = ['t-val',num2str(db3) '.fig'];
%[F,A,C,I] = plot_FNC(rot90(fliplr(datt)), CLIM, LABEL, RSN_I,'',T,MOD,cm,pvals,0); % db3 as last argument for normal plot 
[F2,A2,C2,I2] = plot_FNC(FNC2, CLIM, LABEL, RSN_I,'',T,MOD,cm,log(SZ_i+eps),0);

%% differeces 
FF=zeros(47);
FFF =zeros(47);
FF(idx) = HC_i;
pal1 = rot90(fliplr(FF));
FF = FF+pal1;
FFF(idx) = SZ_i;
pal2 = rot90(fliplr(FFF));
FFF = FFF +pal2;
F1_F2   = FF-FFF; 
clear FNCC
FNCC = F1_F2;
%FNCC = log(F1_F2+eps);
CLIM = [-max(abs(FNCC(:))) max(abs(FNCC(:)))];
T = 'HC-SZ activity difference';
%fName = ['t-val',num2str(db3) '.fig'];
%[F,A,C,I] = plot_FNC(rot90(fliplr(datt)), CLIM, LABEL, RSN_I,'',T,MOD,cm,pvals,0); % db3 as last argument for normal plot 
[F,A,C,I] = plot_FNC(FNCC, CLIM, LABEL, RSN_I,'',T,MOD,cm,FNCC,0);

%% ranking for kendall tau
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/SZPaircount.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD/HCPaircount.mat

% sort SZ
[~,sortedSZ_idx] = sort(SZ_i,'descend');
%sortedSZs = SZ_i(sortedSZ_idx);        % SZ sorted hight to low counts  
%sortedHCs = HC_i(sortedSZ_idx);        % HC sorted accrodingly
%sortedSZ_HC = [sortedSZs sortedHCs];

% sort HC
[~,sortedHC_idx] = sort(HC_i,'descend');
ranking_sSZ_HC = zeros(1081,2);
ranking_sHC_SZ = zeros(1081,2);

for i = 1:1081
    ranking_sSZ_HC(i,1) = i;
    ranking_sSZ_HC(i,2) = find(sortedHC_idx==sortedSZ_idx(i)); % According to SZ
    
    ranking_sHC_SZ(i,1) = i;
    ranking_sHC_SZ(i,2) = find(sortedSZ_idx==sortedHC_idx(i)); % According to HC 
end

%sortedHCs1 = HC_i(sortedHC_idx);        % HC sorted hight to low counts  
%sortedSZs1 = SZ_i(sortedHC_idx);        % SZ sorted accrodingly 
%sortedHC_SZ = [sortedSZs1 sortedHCs1]; 
%% Z score calculation using tau
T = 0.011094; % updated EMD 
%T = 0.045373;   % HC-sorted or SZ sorted 
%T = 0.012992 ; % Overall
z = (3*T*sqrt(1081*1080))/sqrt(2*(2*1081+5));

%% dig into the most frequent shapes 
clear
clc
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBEMD.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/shapestosubjectclassMapping.mat % shapes represnts which group SZ or HC
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allshapestoPairmapping.mat % shapes coming from which pair 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/SZPaircounts.mat % counts SZ: how many times occured
load /Users/mrahaman1/Documents/Statelet_V2/fResults/HCPaircount.mat  % counts HC: how many times occured
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/pDensity_allDOMshapes.mat % probability 
SZshapes_Pairs = shape_to_pair(shape_to_subcl==0);
HCshapes_Pairs = shape_to_pair(shape_to_subcl==1);
SZshapesIdx = find(shape_to_subcl==0);
HCshapesIdx = find(shape_to_subcl==1);

%[~,sortedSZ_idx] = sort(SZ_i,'descend');
%domshapes_occ_SZ = find(SZshapes_Pairs==sortedSZ_idx(1)); % Pair 519 occur[441,1661,1673,2180,2464,2569,2570,2571,2572,2729]
%realshapesmapping_SZ = SZshapesIdx(domshapes_occ_SZ); 
%rsm_SZ = allSUBshapes(realshapesmapping_SZ,:);
%% HC sorting w.r.t. most dominant 
[~,sortedHC_idx] = sort(HC_i,'descend');
domshapes_occ_HC = find(HCshapes_Pairs==sortedHC_idx(1)); % pair 38 occur[167,1392,1513,1514,1923,2372,3066,3298,4330,4389]
realshapesmapping_HC = HCshapesIdx(domshapes_occ_HC); 
%rsm_HC = allSUBshapes(realshapesmapping_HC,:);
%meanrsm_HC = mean(rsm_HC);
[~,mdx] = max(p(realshapesmapping_HC)); % 4330th shape
sims_for_dHC = simscore(realshapesmapping_HC(mdx),HCshapesIdx);
otherPair = setdiff(1:1081,shape_to_pair(realshapesmapping_HC(mdx)));
[~,sorted_shapes_HC] = sort(sims_for_dHC);
flags = zeros(size(sorted_shapes_HC,2),1);

dS_occ_HC = find(HCshapes_Pairs==HCshapes_Pairs(sorted_shapes_HC(1))); % initialize the flags for the dominant pair [pair to shape id mapping]  
flags(dS_occ_HC) = 1;
count =1;

for i = 2:size(sorted_shapes_HC,2)
    if(flags(sorted_shapes_HC(i))==0)
    distances(count) = sims_for_dHC(i);
    ipair = HCshapes_Pairs(sorted_shapes_HC(i)); % shape id to pair 
    dS_occ_HC = find(HCshapes_Pairs==ipair); 
    flags(dS_occ_HC) = 1;
    dSmapping_HC = HCshapesIdx(dS_occ_HC);
    [~,pds] = max(p(dSmapping_HC));
    pairshapes(count) = dSmapping_HC(pds); % store the shape
    mappedpairs(count) = ipair;            % store the pair 
    
    count = count+1;
    end
end
HC_sortedMatch = allSUBshapes(pairshapes,:);

%% Identical frequency testing

clear
clc

%{
load /Users/mrahaman1/Documents/Statelet_V2/fResults/shapestosubjectclassMapping.mat % shapes represnts which group SZ or HC
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allshapestoPairmapping.mat % shapes coming from which pair 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/SZPaircounts.mat % counts SZ: how many times occured
load /Users/mrahaman1/Documents/Statelet_V2/fResults/HCPaircount.mat  % counts HC: how many times occured
load /Users/mrahaman1/Documents/Statelet_V2/fResults/pDensity_allDOMshapes.mat % probability 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes_real.mat
%scatter(SZ_i,HC_i)
%}
% freq[SZ,HC]= [6,6],[5,5],[4,4] [5,4] [4,5],[3,3]
% freq[SZ,HC]= [15,15],[14,14],[12,12][10,10] 
cd('/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD')
load allSUBshapes.mat
load allSUBshapes_PD.mat
load allSUBsimscore.mat
load SZPaircount.mat
load HCPaircount
for i=1:length(allSUBshapes)
shape_to_subcl(i) = allSUBshapes{i,4};
shape_to_pair(i)  = allSUBshapes{i,3};
allSUBshapes_real{i} = allSUBshapes{i,1};
end

equalappearance = [9,8,7];
for i=1:length(equalappearance)
    SZ_f_i  = find(SZ_i == equalappearance(i));
    HC_f_i  = find(HC_i == equalappearance(i));
    f_SZ_HC = intersect(SZ_f_i,HC_f_i);  
    figure()
    for j=1:length(f_SZ_HC)
        subplot(5,5,j)    
    %subplot(length(f_SZ_HC),length(f_SZ_HC),j)    
    shapes_SZ_HC = find(shape_to_pair==f_SZ_HC(j));
    groups   =  shape_to_subcl(shapes_SZ_HC);
    
    SZidx   = find(groups==0);
    shapes_SZ = shapes_SZ_HC(SZidx);
    [~,hfreqshapes_SZ_i] = max(p(shapes_SZ));
    hfs_SZ(j) = shapes_SZ(hfreqshapes_SZ_i);
    
    HCidx   = find(groups==1);
    shapes_HC = shapes_SZ_HC(HCidx);
    [~,hfreqshapes_HC_i] = max(p(shapes_HC));
    hfs_HC(j) = shapes_HC(hfreqshapes_HC_i);
    
    tempSZ = allSUBshapes_real{hfs_SZ(j)};
    %figure()
    plot(tempSZ,'r');
    hold on
    tempHC = allSUBshapes_real{hfs_HC(j)};
    plot(tempHC,'b');
    maxl = max(length(tempSZ),length(tempHC));
    plot(zeros(maxl,1),'k')
    box off
    axis tight
    axis off
    end
    
 
end

% 10-10
SZ_with_f10 = find(SZ_i==13);
HC_with_f10 = find(HC_i==13);
f10_10       = intersect(SZ_with_f10,HC_with_f10);    


for i=1:length(f10_10)
    shapes10_10 = find(shape_to_pair==f10_10(i));
    groups   =  shape_to_subcl(shapes10_10);
    
    SZidx   = find(groups==0);
    shapes10_10_SZ = shapes10_10(SZidx);
    [~,hfreqshapes10_10_SZ_i] = max(p(shapes10_10_SZ));
    hfs10_10_SZ(i) = shapes10_10_SZ(hfreqshapes10_10_SZ_i);
    
    HCidx   = find(groups==1);
    shapes10_10_HC = shapes10_10(HCidx);
    [~,hfreqshapes10_10_HC_i] = max(p(shapes10_10_HC));
    hfs10_10_HC(i) = shapes10_10_HC(hfreqshapes10_10_HC_i);
end


% 5-5
SZ_with_f5 = find(SZ_i==5);
HC_with_f5 = find(HC_i==5);
f5_5       = intersect(SZ_with_f5,HC_with_f5);  
for i=1:length(f5_5)
    shapes5_5 = find(shape_to_pair==f5_5(i));
    groups   =  shape_to_subcl(shapes5_5);
    
    SZidx   = find(groups==0);
    shapes5_5_SZ = shapes5_5(SZidx);
    [~,hfreqshapes5_5_SZ_i] = max(p(shapes5_5_SZ));
    hfs5_5_SZ(i) = shapes5_5_SZ(hfreqshapes5_5_SZ_i);
    
    HCidx   = find(groups==1);
    shapes5_5_HC = shapes5_5(HCidx);
    [~,hfreqshapes5_5_HC_i] = max(p(shapes5_5_HC));
    hfs5_5_HC(i) = shapes5_5_HC(hfreqshapes5_5_HC_i);
end


% 4-4
SZ_with_f4 = find(SZ_i==4);
HC_with_f4 = find(HC_i==4);
f4_4       = intersect(SZ_with_f4,HC_with_f4);  
for i=1:length(f4_4)
    shapes4_4 = find(shape_to_pair==f4_4(i));
    groups   =  shape_to_subcl(shapes4_4);
    
    SZidx   = find(groups==0);
    shapes4_4_SZ = shapes4_4(SZidx);
    [~,hfreqshapes4_4_SZ_i] = max(p(shapes4_4_SZ));
    hfs4_4_SZ(i) = shapes4_4_SZ(hfreqshapes4_4_SZ_i);
    
    HCidx   = find(groups==1);
    shapes4_4_HC = shapes4_4(HCidx);
    [~,hfreqshapes4_4_HC_i] = max(p(shapes4_4_HC));
    hfs4_4_HC(i) = shapes4_4_HC(hfreqshapes4_4_HC_i);
end

% 4-5
SZ_with_f4 = find(SZ_i==4);
HC_with_f5 = find(HC_i==5);
f4_5       = intersect(SZ_with_f4,HC_with_f5);  
for i=1:length(f4_5)
    shapes4_5 = find(shape_to_pair==f4_5(i));
    groups   =  shape_to_subcl(shapes4_5);
    
    SZidx   = find(groups==0);
    shapes4_5_SZ = shapes4_5(SZidx);
    [~,hfreqshapes4_5_SZ_i] = max(p(shapes4_5_SZ));
    hfs4_5_SZ(i) = shapes4_5_SZ(hfreqshapes4_5_SZ_i);
    
    HCidx   = find(groups==1);
    shapes4_5_HC = shapes4_5(HCidx);
    [~,hfreqshapes4_5_HC_i] = max(p(shapes4_5_HC));
    hfs4_5_HC(i) = shapes4_5_HC(hfreqshapes4_5_HC_i);
end

% 5-4
SZ_with_f5 = find(SZ_i==5);
HC_with_f4 = find(HC_i==4);
f5_4       = intersect(SZ_with_f5,HC_with_f4);  
for i=1:length(f5_4)
    shapes5_4 = find(shape_to_pair==f5_4(i));
    groups   =  shape_to_subcl(shapes5_4);
    
    SZidx   = find(groups==0);
    shapes5_4_SZ = shapes5_4(SZidx);
    [~,hfreqshapes5_4_SZ_i] = max(p(shapes5_4_SZ));
    hfs5_4_SZ(i) = shapes5_4_SZ(hfreqshapes5_4_SZ_i);
    
    HCidx   = find(groups==1);
    shapes5_4_HC = shapes5_4(HCidx);
    [~,hfreqshapes5_4_HC_i] = max(p(shapes5_4_HC));
    hfs5_4_HC(i) = shapes5_4_HC(hfreqshapes5_4_HC_i);
end

% 3-3
SZ_with_f3 = find(SZ_i==3);
HC_with_f3 = find(HC_i==3);
f3_3       = intersect(SZ_with_f3,HC_with_f3);  
for i=1:length(f3_3)
    shapes3_3 = find(shape_to_pair==f3_3(i));
    groups   =  shape_to_subcl(shapes3_3);
    
    SZidx   = find(groups==0);
    shapes3_3_SZ = shapes3_3(SZidx);
    [~,hfreqshapes3_3_SZ_i] = max(p(shapes3_3_SZ));
    hfs3_3_SZ(i) = shapes3_3_SZ(hfreqshapes3_3_SZ_i);
    
    HCidx   = find(groups==1);
    shapes3_3_HC = shapes3_3(HCidx);
    [~,hfreqshapes3_3_HC_i] = max(p(shapes3_3_HC));
    hfs3_3_HC(i) = shapes3_3_HC(hfreqshapes3_3_HC_i);
end
hfs_SZ_all = [hfs10_10_SZ hfs5_5_SZ hfs4_4_SZ hfs4_5_SZ hfs5_4_SZ];
hfs_HC_all = [hfs10_10_HC hfs5_5_HC hfs4_4_HC hfs4_5_HC hfs5_4_HC];
%$$$$$$$$$$$$$$$$$$$$$$$$$ plot this
figure()
pl = [1,7,8,13,14,15,16,17,19,20,21,22,25,26,27,28,29,30]; % subplot ID
 for i =1:length(hfs_HC_all)
 subplot(6,6,pl(i))
 tempSZ = allSUBshapes_real{hfs_SZ_all(i)};
 %tempSZ = allSUBshapes(hfs_SZ_all(i),:);
 plot(tempSZ,'r');
 hold on
 tempHC = allSUBshapes_real{hfs_HC_all(i)};
 %tempHC = allSUBshapes(hfs_HC_all(i),:);
 plot(tempHC,'b');
 %plot(allSUBshapes_real{hfs_HC_all(i)},'b');
 %title(['Pair: ' num2str(shape_to_pair(hfs6_6_HC(i)))]);
 %legend('SZ','HC')
 %set(gca,'FontSize',20)
 maxl = max(length(tempSZ),length(tempHC));
 plot(zeros(maxl,1),'k')
 box off
 axis tight
 axis off
%saveas(F,fullfile(figDir,fName));
%saveas('/Users/mrahaman1/Documents/Statelet_V2/fResults/domshapes/hfs6',gcf);
%close()
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()

for i =1:length(hfs_HC_all)
 %tempSZ = allSUBshapes_real{hfs_SZ_all(i)};
 tempHC = allSUBshapes_real{hfs_HC_all(i)};
 %plot(tempSZ,'r');
 plot(tempHC,'b');
 hold on
end
%{
figure(2)
for i =1:length(hfs4_4_HC)  
subplot(3,2,i)    
plot(allSUBshapes(hfs4_4_SZ(i),:),'r');
hold on
plot(allSUBshapes(hfs4_4_HC(i),:),'b');
title(['Pair: ' num2str(shape_to_pair(hfs4_4_HC(i)))]);
legend('SZ','HC')
set(gca,'FontSize',20)
%saveas(F,fullfile(figDir,fName));
%saveas('/Users/mrahaman1/Documents/Statelet_V2/fResults/domshapes/hfs6',gcf);
%close()
end
%}

%% Try
Ltvals = zeros(1081,1);
Ltvals(288) = 1;
Ltvals(400) = 2;
Ltvals(534) = 3;
Ltvals(656) = 4;
Ltvals(767) = 5;

datt = zeros(47);    
datt(idx) = Ltvals;  % tval
datt2 = rot90(fliplr(datt));
Fdatt = zeros(47);
Fdatt(idx) = Ltvals; 

FNC = datt2+Fdatt; 
CLIM = [-max(abs(Ltvals(:))) max(abs(Ltvals(:)))];
T = 't-values';
[F,A,C,I] = plot_FNC(FNC, CLIM, LABEL, RSN_I,'',T,MOD,cm,Ltvals,0);
%% opposite for HC and SZ yet the same pairs
clear
clc
%{
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBEMD.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/shapestosubjectclassMapping.mat % shapes represents which group SZ or HC
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allshapestoPairmapping.mat % shapes coming from which pair 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes.mat
load /Users/mrahaman1/Documents/Statelet_V2/fResults/pDensity_allDOMshapes.mat % probability 
%}
cd('/Users/mrahaman1/Documents/Statelet_V2/fResults/updatedEMD')
load allSUBshapes.mat
load allSUBshapes_PD.mat
load allSUBsimscore.mat
load SZPaircount.mat
load HCPaircount
for i=1:length(allSUBshapes)
shape_to_subcl(i) = allSUBshapes{i,4};
shape_to_pair(i)  = allSUBshapes{i,3};
allSUBshapes_real{i} = allSUBshapes{i,1};
end

SZshapes_Pairs = shape_to_pair(shape_to_subcl==0);
HCshapes_Pairs = shape_to_pair(shape_to_subcl==1);
SZshapesIdx = find(shape_to_subcl==0);
HCshapesIdx = find(shape_to_subcl==1);
common_Pair = intersect(unique(SZshapes_Pairs),unique(HCshapes_Pairs));

for i =1:length(common_Pair)
    szids = SZshapesIdx((SZshapes_Pairs==common_Pair(i)));
    hcids = HCshapesIdx((HCshapes_Pairs==common_Pair(i)));
    %for ii =1:length(szids)
     %   if(simscore())
    %end
    emdscore = (simscore(szids,hcids));
    emdscore1 = (simscore(hcids,szids));
    maxemd = max(emdscore(:));
    maxemd1 = max(emdscore1(:));
    [x,y]=find(emdscore==maxemd);
    [x1,y1]=find(emdscore1==maxemd1);
    if(maxemd>maxemd1)
      shapespair = [szids(x(1)) hcids(y(1))];
      RevPairshapes{i,1} = shapespair;     % shapes pair show highest emd 
      RevPairemd(i,1) = maxemd;  
    else
      shapespair = [szids(y1(1)) hcids(x1(1))]; 
      RevPairshapes{i,1} = shapespair;     % shapes pair show highest emd 
      RevPairemd(i,1) = maxemd1;  
    end
    %RevPairshapes{i,1} = shapespair;     % shapes pair show highest emd 
    %RevPairemd(i,1) = maxemd;            % highest emd
    RevPair(i,1) = common_Pair(i);       % pair
end

[~,tals] = sort(RevPairemd);
kk = RevPair(tals(1:10));

for t =1:10
ppp(t,:) = RevPairshapes{tals(t),1};
end

figure()
for l =1:length(ppp)
    subplot(2,5,l)
    plot(allSUBshapes_real{ppp(l,1)},'r');
    hold on
    plot(allSUBshapes_real{ppp(l,2)},'b');
    maxl = max(length(allSUBshapes_real{ppp(l,1)}),length(allSUBshapes_real{ppp(l,2)}));
    plot(zeros(maxl,1),'k')
    title(['Pair: ' num2str(shape_to_pair(ppp(l,1)))])
   % box off
    axis tight
   axis off
end

%% subcortical pairs [1,2,3,4,47,48,49,92,93,136]
load /Users/mrahaman1/Documents/Statelet_V2/fResults/pDensity_allDOMshapes.mat % probability 
sc = [1,2,3,4,47,48,49,92,93,136];
HCsc = HC_i(sc);
hsc = sc(HCsc~=0);
kl = [1 2 47];
for c =1:length(kl)
    tilo = find(HCshapes_Pairs == kl(c));
    tsh = HCshapesIdx(tilo);
    [~,tt] = max(p(tsh));
    %hcsc_i(c)
    pal(c) = tsh(tt);
end

SZsc = SZ_i(sc);
ssc = sc(SZsc~=0);

for c=1:length(kl)
    tlo = find(SZshapes_Pairs == (kl(c)));
    ts = SZshapesIdx(tlo);
    [~,ttt] = max(p(ts));
    pal1(c) = ts(ttt);
    %szsc_i(c) = 
end


find(HCshapes_Pairs ==1)

figure()
for l =1:length(pal1)
    subplot(1,3,l)
    plot(allSUBshapes(pal(l),:),'b');
    hold on
    plot(allSUBshapes(pal1(l),:),'r');
end

%%

