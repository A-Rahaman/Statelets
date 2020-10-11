function fr = subjectWISEshapes(s) 
load /Users/mrahaman1/Documents/Statelet_V2/data/candidates_22to50.mat
load /Users/mrahaman1/Documents/mrahaman/StateLet/dFNCs.mat;
dataDir = '/home/users/mrahaman1/Statelet_V2/data/shapelets';
subDir = fullfile(dataDir,['info_SUB_' num2str(s,'%03.f')]);
load(subDir);
subData  = squeeze(rFdyn(s,:,:));
%% For all the pairs of subejct
c = 1; % shapes counter
pairWISELeadershapes = cell(1081,1);
for p =1:size(subData,2)
cands_p = candidates{1,repWavelet_Details{p,1}};
pSig = subData(:,p);                                                % signal of that pair, p 
occur = repWavelet_Details{p,3};                                    % how many times that occures 
pairWISELeadershapes{p} = pSig(cands_p{repWavelet_Details{p,4}});   % leader
    for j =1:length(occur)
        pairWISEshapes{p,j} = pSig(cands_p{occur(j)});              % Pairwise distributed 
        allshapes_s{c,1} = pSig(cands_p{occur(j)});                 % All the shapes 
        c = c+1;
    end
end
save(['/home/users/mrahaman1/Statelet_V2/results/shapelets_arrangedsubWISE/sub_' num2str(s,'%03.f')],'pairWISEshapes');
save(['/home/users/mrahaman1/Statelet_V2/results/allshapesof_a_sub/sub_' num2str(s,'%03.f')],'allshapes_s');
save(['/home/users/mrahaman1/Statelet_V2/results/leaders/sub_' num2str(s,'%03.f')],'allshapes_s');
end
