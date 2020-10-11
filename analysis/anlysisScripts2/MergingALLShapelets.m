clear
clc
addpath(genpath('/export/mialab/users/mrahaman/StateLet/Scripts_Updated'));
addpath(genpath('/export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts'));
load '/export/mialab/users/mrahaman/StateLet/dFNCs.mat'
load /export/mialab/users/mrahaman/StateLet/Data/Candidates_10-50.mat
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/Info_HCs.mat
load /export/mialab/users/mrahaman/StateLet/Data/CtsubjectIdx.mat
%Wavelets = zeros(314,1081,50);
Wavelets = cell(314,1081);
for i =1:size(Info_HCs,1)
    shapelets = Info_HCs{i}.repWavelet_Details;
    data  = squeeze(rFdyn(CtIdx(i),:,:));  % Subject's data   
    
    for j =1:size(shapelets,1)
        data_j            = data(:,j);
        cands_j           = candidates{shapelets{j,1}};
        pp=data_j(cands_j{shapelets{j,4}});
        %Wavelets(CtIdx(i),j,1:length(pp)) = pp';
        Wavelets{CtIdx(i),j}= pp;
        
    end  
end

% SZ's
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/Info_SZs.mat
load /export/mialab/users/mrahaman/StateLet/Data/SZsubjectIdx.mat
for i =1:size(Info_SZs,1)
    shapelets = Info_SZs{i}.repWavelet_Details;
    data  = squeeze(rFdyn(SZIdx(i),:,:));  % Subject's data   
    
    for j =1:size(shapelets,1)
        data_j            = data(:,j);
        cands_j           = candidates{shapelets{j,1}};
        pp=data_j(cands_j{shapelets{j,4}});
        %Wavelets(SZIdx(i),j,1:length(pp)) = pp';
        Wavelets{SZIdx(i),j}= pp;
        
    end  
end

