clear
clc
addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/anlysisScripts'))
load /Users/mrahaman1/Documents/Statelet_V2/fResults/shapestosubjectclassMapping.mat % shapes represnts which group SZ or HC
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allshapestoPairmapping.mat % shapes coming from which pair 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/pDensity_allDOMshapes.mat % probability 
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes_real.mat % real length
%load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBEMD.mat          % EMD          
load /Users/mrahaman1/Documents/Statelet_V2/fResults/allSUBshapes.mat % extrapolated 
shape = initShape();
for i=1:size(p,2)
    shape(i).extrapolated = allSUBshapes(i,:);
    shape(i).real_length= allSUBshapes_real{i};
    shape(i).length = length(allSUBshapes_real{i});   
    shape(i).maxconnectivity= max(abs(allSUBshapes_real{i})); % Need to take abs becuase neagative and positive both are connectivity
    shape(i).meanconnectivity=mean(abs(allSUBshapes_real{i}));
    shape(i).pair = shape_to_pair(i);
    shape(i).subjectclass = shape_to_subcl(i);
    shape(i).probabilitydensity = p(i);                  
end