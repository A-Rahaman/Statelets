clear
clc
addpath(genpath('/Users/mrahaman1/Documents/StateLets/scripts/Scripts_Updated'));
load /Users/mrahaman1/Documents/StateLets/metaData/SimScores_Spearman.mat
load /Users/mrahaman1/Documents/StateLets/metaData/Wavelets_314X1081.mat
genSIMSCORE = zeros(314,1081);
shSim = zeros(314,1081,2);
replicated_in_pairs = cell(314,1081);
%% Get the best shapelets per pairs
%{
% First method for choosing best shapelet

for i =1:size(simSCORE,1)
    for j =1:size(simSCORE,2)
        x = squeeze(simSCORE(i,j,:));
        idx = setdiff(1:1081,j);
        ids = find(x >= 0.5);
        replicated_in_pairs{i,j} = idx(ids); % Ok  
        shSim(i,j,1)     = mean(x(ids));     % 
        shSim(i,j,2)     = length(ids);
        genSIMSCORE(i,j) = (mean(x(ids))+length(ids)/1081)/2; % Effect of similarity values and the number of pairs it covers. 
    end
end
save('genSIMSCORE_SM','genSIMSCORE','-v7.3');
save('simScore-popularity_SM','shSim','-v7.3');
save('replications_SM','replicated_in_pairs','-v7.3');
%}

% Second method based on replicability of the pair  
for i =1:size(simSCORE,1)
    for j =1:size(simSCORE,2)
        x = squeeze(simSCORE(i,j,:));
        idx = setdiff(1:1081,j);
        ids = find(x >= 0.7); % See the similarity/Correlation 
        replicated_in_pairs{i,j} = idx(ids);   
        shSim(i,j,1)     = mean(x(ids));     
        shSim(i,j,2)     = length(ids);
        %genSIMSCORE(i,j) = (mean(x(ids))+length(ids)/1081)/2; % Effect of similarity values and the number of pairs it covers. 
        genSIMSCORE(i,j) = length(ids);
    end
end

save('genSIMSCORE_REP','genSIMSCORE','-v7.3');
save('simScore-popularity_REP','shSim','-v7.3');
save('replications_REP','replicated_in_pairs','-v7.3');
%% Pick 5 most diversified shapelets 
clear
clc


%{
Get Five diversified 
load /Users/mrahaman1/Documents/StateLets/metaData/genSIMSCORE_SM
load /Users/mrahaman1/Documents/StateLets/metaData/Wavelets_314X1081.mat
load /Users/mrahaman1/Documents/StateLets/metaData/replications_SM.mat

SUB_shLets_idx = zeros(314,5);
SUB_wise_shapelets = zeros(314,5,28);
SUB_wise_shapelets_uLen = cell(314,5);
reps_COMM = cell(314,5);

for i = 1:size(genSIMSCORE,1)
    i
    %fprintf("Subject %d \n",i);    
    sims_i = genSIMSCORE(i,:);
    %reps_i = shSim(i,:,2);
    %[~,idss] = maxk(sims_i,5); % Do something here keep diversity 
    %%% diversity %%%%
    [vals,inds] = sort(sims_i,'descend');
    K_simscores(1) = vals(1); 
    idss(1)        = inds(1);
    for j = 2:length(vals)
        flag = 1;
        for l = 1:length(idss)
            %%% generalized length
               x = Wavelets{i,inds(j)}; % The one being tested  
               y = Wavelets{i,idss(l)}; % lth enlisted 
               % Check if y is longer
               if length(x) < length(y)
                 x = interp1( x, linspace( 1, length(x), length(y) ) );   % Resample x
                 x = x';
               else
                  y = interp1( y, linspace( 1, length(y), length(x) ) );   % Resample y
                  y = y';
               end
           %%%
            if(corr(x,y,'Type','Spearman')> 0.5)
            %fprintf("Correlation high at %d \n",j);    
            flag=0;
            break;
            end
        end
        if(flag~=0)
            K_simscores(end+1) = vals(j); 
            idss(end+1)        = inds(j);
        end
        if(length(idss)==5)
            break;
        end
    end
    %%%%
    SUB_shLets_idx(i,:) = idss; % Index of most popular shapelets of subject i 
    for k = 1:length(idss)
        % The median of the all the lengths is 28. So, we would interpolate/extrapolate all the signal to that lengths  
         x = Wavelets{i,idss(k)}; 
         if length(x) < 28
            x = interp1( x, linspace( 1, length(x), 28 ) );   % Resample x
         elseif length(x) > 28
            x = imresize(x',[1,28]);                          % Downsampling x
         end
       SUB_wise_shapelets_uLen{i,k} = Wavelets{i,idss(k)}; 
       SUB_wise_shapelets(i,k,:)= x;
       reps_COMM{i,k} = union(replicated_in_pairs{i,idss(k)},idss(k)); % Put the pair itself in replications list
    end
end
save('SM_replicability_of_Leaders','reps_COMM','-v7.3');
save('SM_SUB_wise_leaderShapelets_uLen','SUB_wise_shapelets_uLen','-v7.3');
save('SM_Subject_wise_shIndex','SUB_shLets_idx','-v7.3');
save('SM_SUB_wise_leaderShapelets','SUB_wise_shapelets','-v7.3');
%}

% Get Five with max diversity in term of replicability 
load /Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/genSIMSCORE_REP
load /Users/mrahaman1/Documents/StateLets/metaData/Wavelets_314X1081.mat
load /Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/replications_REP.mat

SUB_shLets_idx = zeros(314,5);
SUB_wise_shapelets = zeros(314,5,28);
SUB_wise_shapelets_uLen = cell(314,5);
reps_COMM = cell(314,5);
for i = 1:size(genSIMSCORE,1)
    i
    %fprintf("Subject %d \n",i);    
    sims_i = genSIMSCORE(i,:);
    %reps_i = shSim(i,:,2);
    [~,idss] = maxk(sims_i,5); % Do something here to keep diversity 
    SUB_shLets_idx(i,:) = idss; % Index of most popular shapelets of subject i 
    for k = 1:length(idss)
        % The median of the all the lengths is 28. So, we would interpolate/extrapolate all the signal to that lengths  
         x = Wavelets{i,idss(k)}; 
         if length(x) < 28
            x = interp1( x, linspace( 1, length(x), 28 ) );   % Resample x
         elseif length(x) > 28
            x = imresize(x',[1,28]);                          % Downsampling x
         end
       SUB_wise_shapelets_uLen{i,k} = Wavelets{i,idss(k)}; 
       SUB_wise_shapelets(i,k,:)= x;
       reps_COMM{i,k} = union(replicated_in_pairs{i,idss(k)},idss(k)); % Put the pair itself in replications list
    end
end
save('REP_replicability_of_Leaders','reps_COMM','-v7.3');
save('REP_SUB_wise_leaderShapelets_uLen','SUB_wise_shapelets_uLen','-v7.3');
save('REP_Subject_wise_shIndex','SUB_shLets_idx','-v7.3');
save('REP_SUB_wise_leaderShapelets','SUB_wise_shapelets','-v7.3');
%% testing the optimized k

% Optimize for best K and it is 5 as well. 
clear
clc
%load /Users/mrahaman1/Documents/StateLets/metaData/SM_SUB_wise_leaderShapelets.mat % shapelets
%load /Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/REP_SUB_wise_leaderShapelets.mat
load /Users/mrahaman1/Documents/Statelet_V2/data/EMDresults/allSUBshapes.mat
%data = reshape(SUB_wise_shapelets,1570,28);
data = allSUBshapes;
Ks = 2:20;
for k = 1:length(Ks)   

fprintf("Running for k = %d\n",Ks(k));
% Method 1 : Elbow method for determining optimal k
clear C IDX SUMD D      
[IDX, C, SUMD, D]= kmeans(data, Ks(k),'Distance','correlation','MaxIter',10000,'Replicates',5);
idx{k} = IDX;
distance = zeros(size(D,2),1);
% SSE Calculation for all k clusters
for s = 1:size(D,2)
windows_for_s_cluster = find(IDX==s); % all windows of that type
alldist = D(:,s); % all distances for that cluster
dist_for_sclst_windows = alldist(windows_for_s_cluster); % Distances for that type of windows 
distance(s) = sum(dist_for_sclst_windows.^2); % Squared sum of distances for each cluster
end
DIST{k} = distance; % Storing all distances for each k value
SSE(k) = sum(distance); % For each value of k
end
% Method 2: Silhouette score calculation for determining optimal k
%{
clear idx5 silh5
% City block is the manhattan distance 
idx5 = kmeans(X,Ks(k),'Distance','correlation');
idx{k} = idx5;
figure
[silh5,h] = silhouette(X,idx5,'correlation');
silhouette_values{k} = silh5;
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value'
ylabel 'Cluster'
saveas(h,['K',num2str(Ks(k))]);
close
%}
%}
%% Running the final k-means on subjectwise shapelets  

clear
clc
%load /Users/mrahaman1/Documents/StateLets/metaData/SM_SUB_wise_leaderShapelets.mat % shapelets
%load /Users/mrahaman1/Documents/StateLets/metaData/SM_replicability_of_Leaders.mat

load /Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/REP_SUB_wise_leaderShapelets.mat
load /Users/mrahaman1/Documents/StateLets/metaData/Replicability_based_sim/REP_replicability_of_Leaders.mat
data = reshape(SUB_wise_shapelets,1570,28);
k = 7;
clear C IDX SUMD D  
%[IDX, C, SUMD, D]= kmeans(data, k,'Distance','sqeuclidean','MaxIter',10000,'Replicates',100);
[IDX, C, SUMD, D]= kmeans(data, k,'Distance','correlation','MaxIter',10000,'Replicates',100);
IDXX = reshape(IDX,314,5); % Return back to the original data dimension 

Statelets = cell(k,3);

for i = 1:k
[s,sh]= find(IDXX==i);    
pairs = reps_COMM{s(1),sh(1)};
for p =1:length(s)
    pairs = union(pairs,reps_COMM{s(p),sh(p)});
    %pairs = intersect(pairs,reps_COMM{s(p),sh(p)});
end
Statelets{i,1} = s;     % Subejcts
Statelets{i,2} = sh;    % Shapelets 
Statelets{i,3} = pairs; % Pairs 
end

%%
figure(1)
subplot(2,4,1)
plot(C(1,:))
title('1')

subplot(2,4,2)
plot(C(2,:))
title('2')

subplot(2,4,3)
plot(C(3,:))
title('3')

subplot(2,4,4)
plot(C(4,:))
title('4')

subplot(2,4,8)
plot(C(5,:))
title('5')

subplot(2,4,7)
plot(C(6,:))
title('6')

subplot(2,4,6)
plot(C(7,:))
title('7')
% Subjects in Statelets for CORR_K7
Statelet - 1 ---- HC : 63, SZ : 67
           2 ---- HC : 74, SZ : 74
           3 ---- HC : 130 SZ : 115
           4 ---- HC : 127   132
           5 ----  39 29
           6 ----  70 68
           7 ----  65 57
           
           
