clear
clc
dataDir = '/Users/mrahaman1/Documents/Statelet_V2/fResults/subWISEEMDMats';          % distances 
dirs    = dir(dataDir);
ids = zeros(size(dirs,1)-2,1);
for i = 3:size(dirs,1)
    naam = dirs(i).name;
    id_of_naam = naam(12:end-4);
    ids(i-2) = str2num(id_of_naam);
end
remids = setdiff(1:314,ids);