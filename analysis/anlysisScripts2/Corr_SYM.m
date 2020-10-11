clear
clc
load /Users/mrahaman1/Documents/StateLets/metaData/Statelets_SM_k5.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/SZsubjectIdx.mat
load /Users/mrahaman1/Documents/mrahaman/d_NBiC/Data/PanssScore_allSubjects.mat
load /Users/mrahaman1/Documents/StateLets/metaData/SM_SUB_wise_leaderShapelets_uLen.mat


for i = 1:size(Statelets,1)
SUB =  unique(Statelets{i,1});
SZSUB = intersect(SZIdx,SUB);  
PPos = PanssScore(:,1);
PNeg = PanssScore(:,2);
PGen = PanssScore(:,3);
PosScore = PPos(SZSUB); % for all subjects
%NegScore = PNeg(SZSUB); % for all subjects
%GenScore = PGen(SZSUB); % for all subjects
idx = find(isnan(PosScore)==1);
if(~isempty(idx))
    SZSUB(idx(:)) = [];
end

PosScore = PPos(SZSUB); % for all subjects
NegScore = PNeg(SZSUB); % for all subjects
GenScore = PGen(SZSUB); % for all subjects
mean_sub_shapelet = zeros(length(SZSUB),1);
for sz =1:length(SZSUB)
    for sh = 1:size(SUB_wise_shapelets_uLen,2)
        mean_shapelet(sh) = mean(SUB_wise_shapelets_uLen{SZSUB(sz),sh});
    end
    mean_sub_shapelet(sz,1) = mean(mean_shapelet);
end
fprintf("%d  %d",length(mean_sub_shapelet),length(PosScore));
 
corr_StateLes(i,1) = corr(mean_sub_shapelet,PosScore,'Type','Spearman');
corr_StateLes(i,2) = corr(mean_sub_shapelet,NegScore,'Type','Spearman');
corr_StateLes(i,3) = corr(mean_sub_shapelet,GenScore,'Type','Spearman');
end

figure(20)
for k = 1:size(corr_StateLes,2)
plot(corr_StateLes(:,k))
hold on
end


    