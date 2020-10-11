clear
clc
addpath(genpath('/export/mialab/users/mrahaman/StateLet/Scripts_Updated'));
load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts/Wavelets_314X1081.mat
shapes = zeros(314,1081,50);
for i =1:size(Wavelets,1)
    for j =1:size(Wavelets,2)
        x = Wavelets{i,j}; 
        if length(x) < 50
            x = interp1( x, linspace( 1, length(x), 50 ) );   % Resample x
        end
       shapes(i,j,:)= x;
    end
end
save('shapes_314x1081x50','shapes','-v7.3');

data1 = reshape(shapes,314*1081,50);
[P, cen, SUMDD, DD]= kmeans(data1, 5,'Distance','correlation','MaxIter',10000,'Replicates',100);
PMat = reshape(P,314,1081);
%[idx,C,sumd,D] = kmeans(X,20,'Distance','correlation','MaxIter',10000,'Replicates',100);
