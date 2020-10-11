simSCOREs = zeros(314,1081,1080);
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts/simSCORE_1-184.mat
for i = 1:181
    for j =1:1081
        x = squeeze(simSCORE(i,j,:));
        k = find(isnan(x)==0);
        simSCOREs(i,j,:) = x(k);
    end
end
% Load the files one by one and perform the addition.  
simSCOREs = simSCOREs+simSCORE;
save('simSCOREs','simSCOREs','-v7.3');