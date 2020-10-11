%% Getting top 10 shapes with max probability density  
clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load allPairshapes_HC.mat;
for kl =1:size(allPairshapes,1)
shapes_HC{kl} = allPairshapes{kl,1};
end
load allPairshapes_SZ.mat
for kl =1:size(allPairshapes,1)
shapes_HC{kl} = allPairshapes{kl,1};
end

load ('allPairshapesPD_HC.mat');
p_HC = p;
load ('allPairshapesPD_SZ.mat');
p_SZ = p;

[~,hidx] = sort(p_HC,'descend');
t10_p_HC = hidx(1:10);

[~,sidx] = sort(p_SZ,'descend');
t10_p_SZ = sidx(1:10);

figure()
for k=1:length(t10_p_SZ)
    subplot(2,5,k)
    tempSZ = shapes_HC{t10_p_SZ(k)};
    plot(tempSZ,'r');
    hold on
    tempHC = shapes_HC{t10_p_HC(k)};
    plot(tempHC,'b');
    maxl = max(length(tempSZ),length(tempHC));
    plot(zeros(maxl,1),'k')
    box off
    axis tight
    axis off 
end
%% get top 10 high frequency pair from HC and just get the leader for those pair from SZ

clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load allPairshapes_HC.mat;
for kl =1:size(allPairshapes,1)
shapes_HC{kl} = allPairshapes{kl,1};
pair_HC(kl) = allPairshapes{kl,4};
end
load ('allPairshapesPD_HC.mat');
p_HC = p;
[~,hidx] = sort(p_HC,'descend');
t10_s_HC = hidx(1:10);
t10_p_HC = pair_HC(t10_s_HC);
load pairwiseleaders_SZ.mat

figure()
for k=1:length(t10_p_HC)
    subplot(2,5,k)
    tempSZ = pairwiseleader_SZ{t10_p_HC(k)}; 
    plot(tempSZ,'r');
    hold on
    tempHC = shapes_HC{t10_s_HC(k)};
    plot(tempHC,'b');
    maxl = max(length(tempSZ),length(tempHC));
    plot(zeros(maxl,1),'k')
    title(['Pair: ' num2str(t10_p_HC(k))]);
    box off
    axis tight
    axis off 
end


%% SZ sorted 

clear
clc
cd ('/Users/mrahaman1/Documents/Statelet_V2/fResults/PairWiseAnalysis')
load allPairshapes_SZ.mat;
for kl =1:size(allPairshapes,1)
shapes_SZ{kl} = allPairshapes{kl,1};
pair_SZ(kl) = allPairshapes{kl,4};
end
load ('allPairshapesPD_SZ.mat');
p_SZ = p;
[~,sidx] = sort(p_SZ,'descend');
t10_s_SZ = sidx(1:10);
t10_p_SZ = pair_SZ(t10_s_SZ);
load pairwiseleaders_HC.mat

figure()
for k=1:length(t10_p_SZ)
    subplot(2,5,k)
    tempHC = pairwiseleader_HC{t10_p_SZ(k)}; 
    plot(tempHC,'b');
    hold on
    tempSZ = shapes_SZ{t10_s_SZ(k)};
    plot(tempSZ,'r');
    maxl = max(length(tempSZ),length(tempHC));
    plot(zeros(maxl,1),'k')
    title(['Pair: ' num2str(t10_p_SZ(k))]);
    box off
    axis tight
    axis off 
end
%%

