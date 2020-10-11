
%% HC merging
clear
clc
Info_HCs = {};
%HCDir = 'T:\mialab\users\mrahaman\StateLet\WeightedDrift_implementation\FullRun\Results_HC';
HCDir = '/export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/Results_HC';
HCDirs = dir(HCDir);
HCFreq = zeros(1081,length(HCDirs)-2);
HCavglen = zeros(1081,length(HCDirs)-2);
i=1;
for t = 3:length(HCDirs)
HC_i = load(fullfile(HCDirs(t).folder,HCDirs(t).name));    
for j = 1:size(HC_i.repWavelet_Details,1)
HCFreq(j,i) = HC_i.repWavelet_Details{j,5}; % Frquency in 5th field 
HCavglen(j,i) = HC_i.repWavelet_Details{j,6}; % Average length of shapelet is in 6th field 
end
Info_HCs {i,1} = HC_i.repWavelet_Details;
i=i+1;
end
cd 'T:\mialab\users\mrahaman\StateLet\WeightedDrift_implementation\FullRun\ViZResults_HC_SZ';
save('ShapeletFreq_HC','HCFreq');
save('ShapeletMeanLength_HC','HCavglen');
save('Info_HCs','Info_HCs');

%%  SZ merging 
Info_SZs = {};
SZDir = 'T:\mialab\users\mrahaman\StateLet\WeightedDrift_implementation\FullRun\Results_SZ';
SZDirs = dir(SZDir);
SZFreq = zeros(1081,length(SZDirs));
SZavglen = zeros(1081,length(SZDirs));
i=1;
for t = 3:length(SZDirs)
SZ_i = load(fullfile(SZDirs(t).folder,SZDirs(t).name));    
for j = 1:size(SZ_i.repWavelet_Details,1)
SZFreq(j,i) = SZ_i.repWavelet_Details{j,5};
SZavglen(j,i) = SZ_i.repWavelet_Details{j,6};
end
Info_SZs {i,1} = SZ_i.repWavelet_Details;
i=i+1;
end
cd 'T:\mialab\users\mrahaman\StateLet\WeightedDrift_implementation\FullRun\ViZResults_HC_SZ';
save('ShapeletFreq_SZ','SZFreq');
save('ShapeletMeanLength_SZ','SZavglen');
save('Info_SZs','Info_SZs');
%%