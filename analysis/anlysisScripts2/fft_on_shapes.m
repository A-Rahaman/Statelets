clear
clc

%load /Users/mrahaman1/Documents/StateLets/metaData/Wavelets_314X1081.mat
load /Users/mrahaman1/Documents/StateLets/metaData/SM_SUB_wise_leaderShapelets_uLen.mat
Wavelets = SUB_wise_shapelets_uLen;
for i=1:size(Wavelets,1)
    i
    for j=1:size(Wavelets,2)
        w = Wavelets{i,j} - mean(Wavelets{i,j}); % removing the mean 
        %Freqs(i,j,:) = fft(Wavelets{i,j},10);  % fft on 
        Freqs(i,j,:) = fft(w,9);  % fft on 
    end
end


%%
load /Users/mrahaman1/Documents/Statelet_V2/data/HCsubjectIdx.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/SZsubjectIdx.mat

SZ_info = squeeze(abs(Freqs(SZIdx,:,:)));
HC_info = squeeze(abs(Freqs(CtIdx,:,:)));

for i=1:9
    TMP_HC = HC_info(:,:,i);
    TMP_SZ = SZ_info(:,:,i);
    [~,pvalss(i)] = ttest2(TMP_HC(:),TMP_SZ(:));
end