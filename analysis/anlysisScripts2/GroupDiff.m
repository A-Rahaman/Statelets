%% Compute the pairwise group diffferces of variables max conn, mean conn, pair counts: # of apperances in the grouo dynamics, probability density 
% The probability density by a pair is the sum of PD by all the shapes coming from that pair
% This has been done for group wise dynamics HC/SZ seperate (HC: 9285 and SZ: 8706 shapes) 
clear
clc
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/idx1081.mat  % Lower 1081 idx (idx)
%load /Users/mrahaman1/Documents/mrahaman/d_NBiC/Data/idxUpper.mat % Upper 1081 idx (IND1)
%% Time decays 

load HC_avggraph_edges.mat
load SZ_avggraph_edges.mat
HC_i = HC_edge_weights;
SZ_i = SZ_edge_weights;

%% PD  
load HC_pairwise_probability.mat
load SZ_pairwise_probability.mat


%% Counts
load Peakcounts_HC.mat
load Peakcounts_SZ.mat
% load sz_paircounts_conv.mat
% load hc_paircounts_conv.mat
% SZ_i = sz_paircount;
% HC_i = hc_paircount;

%% Max conn
%load HCPairwisemaxConn.mat
%load SZPairwisemaxConn.mat
load pairwiseleaders_HC.mat
load pairwiseleaders_SZ.mat

HC_i = zeros(1081,1);
SZ_i = zeros(1081,1);
for i =1:size(pairwiseleader_HC,2)
    HC_i(i) = max(abs(pairwiseleader_HC{1,i}));
    SZ_i(i) = max(abs(pairwiseleader_SZ{1,i}));
end
%% Mean Conn
%load HCPairwisemeanConn.mat
%load SZPairwisemeanConn.mat
load pairwiseleaders_HC.mat
load pairwiseleaders_SZ.mat
HC_i = zeros(1081,1);
SZ_i = zeros(1081,1);
for i =1:size(pairwiseleader_HC,2)
    HC_i(i) = mean(pairwiseleader_HC{1,i});
    SZ_i(i) = mean(pairwiseleader_SZ{1,i});
end
%% Compute HC Matrix
datt1 = zeros(47);    
%datt1(idx) = log(HC_i+eps); 
datt1(idx) = (HC_i);  
datt21 = rot90(fliplr(datt1));
Fdatt1 = zeros(47);
%Fdatt1(idx) = log(HC_i+eps);
Fdatt1(idx) = (HC_i);
FNC1 = datt21+Fdatt1; 


% SZ Matrix
datt = zeros(47);    
%datt(idx) = log(SZ_i+eps);  
datt(idx) = (SZ_i);  % tval
datt2 = rot90(fliplr(datt));
Fdatt = zeros(47);
%Fdatt(idx) = log(SZ_i+eps);
Fdatt(idx) = (SZ_i);
FNC2 = datt2+Fdatt; 

% Differece Matrix 
FF=zeros(47);
FF(idx) = HC_i - SZ_i;
FF = rot90(fliplr(FF)); 
pal1 = zeros(47);
pal1(idx) = HC_i - SZ_i;
FNCC = FF+pal1;


%% Plot all
load /Users/mrahaman1/Documents/StateLets/Plotting_Eswar/reqData_for_plot_FNC.mat
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/data/bluewhitered'))
addpath(genpath('/Users/mrahaman1/Documents/StateLets/Plotting_Eswar'));
MOD = fmod_RSN;
LABEL = L;
cm = '';
M1 = max(max(abs(FNC1(:))),max(abs(FNC2(:))));
M2 = max(M1,max(abs(FNCC(:))));
CLIM = [0 M2]; 
%CLIM = [-M2 M2]; 
%CLIM = [0 0.7];
%CLIM = [];
T1 = 'HC';
T2 = 'SZ';
T3 = 'HC-SZ';

%cmap = cold(129);
%cmap = cmap(end:-1,:);
cmap = jet(20);
cmap = flipud(cmap(1:10,:));
cmap(1,:) = [1,1,1];
%colormap(cmap);

[F,A,C,I] = plot_FNC(FNCC, CLIM, LABEL, RSN_I,'',T3,MOD,cm,FNCC,0); 
%colormap(bluewhitered(256)), colorbar
%colormap(cmap), colorbar;
[F1,A1,C1,I1] = plot_FNC(FNC1, CLIM, LABEL, RSN_I,'',T1,MOD,cm,FNC1,0); 
%colormap(bluewhitered(256)), colorbar
%colormap(cmap), colorbar;
[F2,A2,C2,I2] = plot_FNC(FNC2, CLIM, LABEL, RSN_I,'',T2,MOD,cm,FNC2,0);
%colormap(bluewhitered(256)), colorbar
%colormap(cmap), colorbar
%% Groupwise most dominant 
clear
clc
load allPairEMD_HC.mat
load allPairshapes_HC.mat
thHC_SS = (simscore>0 & simscore<0.8); 
sumthHC_SS = sum(thHC_SS);
%HC_shapes_replication_index = sumthHC_SS/9285;
[hcmaxrep,hcidx1] = sort(sumthHC_SS,'descend');
hcmaxrep(1)
hcidx1(1)
%% SZ
load allPairEMD_SZ.mat
load allPairshapes_SZ.mat
thSZ_SS = (simscore>0 & simscore<0.8); 
sumthSZ_SS = sum(thSZ_SS);
%HC_shapes_replication_index = sumthHC_SS/9285;
[szmaxrep,szidx1] =sort(sumthSZ_SS,'descend');
szmaxrep(1)
szidx1(1)
figure()
for k=1:50
subplot(5,10,k)
plot(HC_allPairshapes{hcidx1(k),1},'b')
hold on
plot(SZ_allPairshapes{szidx1(k),1},'r')
box off
axis tight
axis off 
end
%%





