clear
clc
%load C:\Users\arahaman\Dropbox\StateLet\tempScripts\HCNAME.mat
%load C:\Users\arahaman\Dropbox\StateLet\tempScripts\SZNAME.mat
Info_HCs = {};
Info_SZs = {};
%SZDir = 'T:\mialab\users\mrahaman\StateLet\WeightedDrift_implementation\FullRun\Results_SZ';
HCDir = 'T:\mialab\users\mrahaman\StateLet\WeightedDrift_implementation\FullRun\Results_HC';
%SZDirs = dir(SZDir);
HCDirs = dir(HCDir);

HCFreq = zeros(1081,30);
SZFreq = zeros(1081,30);
HCavglen = zeros(1081,30);
SZavglen = zeros(1081,30);

for i =1:30
HC_ref = load(fullfile('T:\mialab\users\mrahaman\StateLet\WeightedDrift_implementation\FullRun\Results_HC',HCNAME{i}));    
SZ_ref = load(fullfile('T:\mialab\users\mrahaman\StateLet\WeightedDrift_implementation\FullRun\Results_SZ',SZNAME{i})); 
for j = 1:1081
HCFreq(j,i) = HC_ref.repWavelet_Details{j,5};   % Frquency in 5th field 
HCavglen(j,i) = HC_ref.repWavelet_Details{j,6}; % Average length of shapelet is in 6th field 
SZFreq(j,i) = SZ_ref.repWavelet_Details{j,5};
SZavglen(j,i) = SZ_ref.repWavelet_Details{j,6}; 
end
Info_HCs {i,1} = HC_ref;
Info_SZs {i,1} = SZ_ref;
end

