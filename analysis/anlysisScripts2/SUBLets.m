clear
clc
addpath(genpath('/export/mialab/users/mrahaman/StateLet/Scripts_Updated'));
load '/export/mialab/users/mrahaman/StateLet/dFNCs.mat';
load /export/mialab/users/mrahaman/StateLet/Data/Candidates_10-50.mat
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/ViZResults_HC_SZ/Info_HCs.mat
load /export/mialab/users/mrahaman/StateLet/Data/CtsubjectIdx.mat
addpath(genpath('/export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts'));
Sublets = zeros(314,1);
oldSigLen = 0;
for i = 1:size(Info_HCs,1)
    tic
    shapelets = Info_HCs{i}.repWavelet_Details;
    data  = squeeze(rFdyn(CtIdx(i),:,:));  % Subject's data
    [mergedSig,tIDS,remLengths] = Chop_and_Merge(data,shapelets,candidates);
    SigLen = 1:length(mergedSig);
    if(length(mergedSig)~= oldSigLen)     % If equal then new candidates aren't necessary.   
    newCandidates = genCandidates(SigLen,50,10);
    end
    uremLengths = unique(remLengths);
    Sublets(CtIdx(i),1) = getShapelets(mergedSig,newCandidates,uremLengths); 
    toc
    oldSigLen = length(mergedSig);
end
