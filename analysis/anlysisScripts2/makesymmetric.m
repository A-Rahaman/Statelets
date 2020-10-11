dataDir = '/Users/mrahaman1/Documents/Statelet_V2/data/EMDresults/EMDscore_subWISE'; % for EMD

for i =1:314
    subDir = fullfile(dataDir,['EMD_of_sub_' num2str(i,'%03.f')]);
    load(subDir);
    simscores = zeros(size(simscore,1));
    for j=1:size(simscore,1)
    idx = setdiff(1:size(simscore,1),j);    
    simscores(j,idx) = simscore(j,:);
    end
    save(['/Users/mrahaman1/Documents/Statelet_V2/data/EMDresults/EMDscore_subWISE_symm/EMD_of_sub_' num2str(i,'%03.f')],'simscores');
    
end