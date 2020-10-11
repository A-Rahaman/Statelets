clear
clc
%load C:\Users\munna\Dropbox\StateLet\Data\Wavelets_314X1081.mat
%load /export/mialab/users/mrahaman/StateLet/WeightedDrift_implementation/FullRun/tempScripts/Wavelets_314X1081.mat
load /Users/mrahaman1/Documents/StateLets/metaData/Wavelets_314X1081.mat
simSCORE = zeros(314,1081,1080);

for i =1:size(Wavelets,1)
    i
    for j=1:size(Wavelets,2)
       x = Wavelets{i,j};
       c=1;
       for k =1:size(Wavelets,2)
           if j~=k
               y = Wavelets{i,k};
               % Check if y is longer
               if length(x) < length(y)
                 x = interp1( x, linspace( 1, length(x), length(y) ) );   % Resample x
                 x= x';
               else
                  y = interp1( y, linspace( 1, length(y), length(x) ) );   % Resample y
                  y = y';
               end
               simSCORE(i,j,c) = corr(x,y,'Type','Spearman');
               c=c+1;
           end
       end
    end
end
save('SimScores_Spearman','simSCORE','-v7.3');