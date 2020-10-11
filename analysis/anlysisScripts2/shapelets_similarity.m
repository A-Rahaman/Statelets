addpath(genpath('/Users/mrahaman1/Documents/Statelet_V2/scripts/shapeletAlgorithmScripts'));       
load /Users/mrahaman1/Documents/Statelet_V2/data/tempLeaders.mat
for j = 922%1081
    %%
    xx = pairWISELeadershapes{j,1};
    %if(length(x)<50)
     %   x = interp1( x, linspace( 1, length(x), 50 ) );   % Resample x
      %  x= x';
    %end
    c = 1;
    strechedshapes = cell(1080,1);
    for k =1:size(pairWISELeadershapes,1)
        x = pairWISELeadershapes{j,1};
         if j~=k
             %y = pairWISELeadershapes{k,1};
             y= pairWISELeadershapes{k,1};
               % Check if y is longer
             if length(x) < length(y)
                 x = interp1( x, linspace( 1, length(x), length(y) ) );   % Resample x
                 x = x';
             else
                  y = interp1( y, linspace( 1, length(y), length(x) ) );   % Resample y
                  y = y';
             end
             x = x-min(x);
             x = x/sum(x);
             y = y-min(y);
             y = y/sum(y);i
             strechedshapes{c,1} = y;
             simscore(c) = dEMD(x,y);
             c=c+1;   
         end
    end
    
    [sortedSim,sidx] = sort(simscore);
    %indx = setdiff(1:1081,j); 
    %sortedrealindx = indx(sidx);
    %%
    %{
    figure();
    for f =1:8    
    subplot(8,1,f)
    plot(strechedshapes{sidx(f)});
    end
    %}
end