clear
clc
load /Users/mrahaman1/Documents/StateLets/metaData/CORR/Statelets_CoRR_K7.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/SZsubjectIdx.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/Data/HCsubjectIdx.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/Wavelets_314X1081X50.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/Wavelets_314X1081.mat
%load /Users/mrahaman1/Documents/StateLets/metaData/SM_SUB_wise_leaderShapelets_uLen.mat
load /Users/mrahaman1/Documents/StateLets/metaData/SM_SUB_wise_leaderShapelets.mat

%SZs = intersect(unique(Statelets{6,1}),SZIdx);
%HCs = intersect(unique(Statelets{6,1}),CtIdx);
%% dwell time

for st = 1:size(Statelets,1)
    %SZs = intersect(unique(Statelets{st,1}),SZIdx);
    %HCs = intersect(unique(Statelets{st,1}),CtIdx);
    subs = Statelets{st,1};
    shapes = Statelets{st,2};    
    c1 = 1;
    c2 = 1;
    dT_SZ = [];
    dT_HC = [];
    for i =1:length(subs)
        if(intersect(subs(i),CtIdx))
        %dT_HC(c1) = length(SUB_wise_shapelets_uLen{subs(i),shapes(i)});
        
        c1 = c1 +1;
        else
            dT_SZ(c2) = length(SUB_wise_shapelets_uLen{subs(i),shapes(i)});
            c2=c2+1;
        end
    
    end
    meandT_SZ(st,1) = mean(dT_SZ);
    meandT_HC(st,1) = mean(dT_HC);
    
end 
%{
for st = 1:size(Statelets,1)
SZs = intersect(unique(Statelets{6,1}),SZIdx);
HCs = intersect(unique(Statelets{6,1}),CtIdx);
figure()
for i =1:length(SZs)
    plot(squeeze(Wavelets(SZs(i),833,:)),'r');
    hold on
end
hold on
for i =1:length(HCs)
    plot(squeeze(Wavelets(HCs(i),833,:)),'b');
    hold on
end
end
%}
%% Average shapelets 
sublets = zeros(314,28);
for i = 1:size(SUB_wise_shapelets,1)
    sublets(i,:) = mean(squeeze(SUB_wise_shapelets(i,:,:)));
end
SZsublets = sublets(SZIdx,:);
HCsublets = sublets(CtIdx,:);
avgSZsublets = mean(SZsublets);
avgHCsublets = mean(HCsublets);

figure();
plot(avgHCsublets,'b');
hold on
plot(avgSZsublets,'r');

%%

