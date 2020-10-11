function fr = getsubWISEshapes_and_EMDscore(s) 
addpath(genpath('/home/users/mrahaman1/Statelet_V2/scripts/shapeletAlgorithmScripts'));
load /home/users/mrahaman1/Statelet_V2/data/candidates_22to50.mat
load /home/users/mrahaman1/Statelet_V2/data/dFNCs.mat;
dataDir = '/home/users/mrahaman1/Statelet_V2/results/EMDresults/shapes_EMD'; % for EMD
subDir = fullfile(dataDir,['info_SUB_' num2str(s,'%03.f')]);
load(subDir);
subData  = squeeze(rFdyn(s,:,:));
%% get the shapes and store, for all the pairs of subejct
c = 1; % shapes counter
pairWISELeadershapes = cell(1081,1);
for p =1:size(subData,2)
cands_p = candidates{1,repWavelet_Details{p,1}};
pSig = subData(:,p);                                                % signal of that pair, p 
occur = repWavelet_Details{p,3};                                    % how many times that occures 
pairWISELeadershapes{p} = pSig(cands_p{repWavelet_Details{p,4}});   % leader
    for j =1:length(occur)
        pairWISEshapes{p,j} = pSig(cands_p{occur(j)});              % Pairwise distributed
        x = pairWISEshapes{p,j};
       if(length(x)<50)
         x = interp1( x, linspace( 1, length(x), 50 ) );            % Resample x
         x = x';
       end
        allshapes_s{c,1} = x;                                        % All the shapes 
        c = c+1;
    end
end
save(['/home/users/mrahaman1/Statelet_V2/results/EMDresults/shapelets_arrangedsubWISE_EMD/sub_' num2str(s,'%03.f')],'pairWISEshapes');
save(['/home/users/mrahaman1/Statelet_V2/results/EMDresults/allshapesof_a_sub_EMD/sub_' num2str(s,'%03.f')],'allshapes_s');
save(['/home/users/mrahaman1/Statelet_V2/results/EMDresults/pairWISELeadershapes_EMD/sub_' num2str(s,'%03.f')],'pairWISELeadershapes');
%% get the EMD scores for each of the shape of a subject
for j = 1:size(allshapes_s,1)
  x = allshapes_s{j,1};  
  c = 1;
  for k = 1:size(allshapes_s,1)
    if j~=k
    y = allshapes_s{k,1};
    x = x-min(x);
    x = x/sum(x);
    y = y-min(y);
    y = y/sum(y);
    simscore(j,c) = dEMD(x,y);
    c = c+1;   
    end 
  end
end
save(['/home/users/mrahaman1/Statelet_V2/results/EMDresults/EMDscore_subWISE/EMD_of_sub_' num2str(s,'%03.f')],'simscore');
end
