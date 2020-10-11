function pr = co_occurences_forserver(i)
%% co-occurences of the shapelets. lets hammer the temporal information of extracted shapelets 
% For the cooccurrences, for a given subject, a given pair, get each occurrence of shapelet and go down across all other pairs 
% using that time information. Count how many other matched - I mean occurred as well probably with a different shape.  
% Fr = #time shapelet of pair j (1 to 1080) occurred/ #time shapelet of pair i occurred  
% Get 1081 X 1080 matrix for each subject  
% Analyze it for HC and SZ separately. Get one mean occurrences across all the HCs and SZs  
% For each group you get a mean occurrence matrix of 1081 X 1080 
% Print that as you print tvals.  
% If you get higher value there that means it is happening across the subjects consistently  
% Hopefully, get some modularity in it.      
% t-test on occurrences for group differences     
% by //munna//
% Dated 11 July 2019
load /home/users/mrahaman1/Statelet_V2/data/candidates_22to50.mat
co_occurences = zeros(1081,1080);
lags = zeros(1081,1080);
dirr = '/home/users/mrahaman1/Statelet_V2/fResults/co-occurences';
dirr1 = '/home/users/mrahaman1/Statelet_V2/fResults/laggings';
load(fullfile('/home/users/mrahaman1/Statelet_V2/fResults/subinfo',['info_SUB_',num2str(i,'%03.f')]));
     Info_i  = repWavelet_Details;
     for j = 1:size(Info_i,1)
        cands_ij = candidates{1,Info_i{j,1}};                              % Candidates for pair j and subject i
        occur_ij = Info_i{j,3};                                            % Occurences for pair j of subject i 
        pairs = 1:1081;
        k_pairs  = pairs(pairs~=j); 
        for k = 1:size(k_pairs,2)
            cands_ik = candidates{1,Info_i{k_pairs(k),1}};                 % Candidates for pair k and subject i
            occur_ik = Info_i{k_pairs(k),3};                               % Occurences for pair k of subject i 
            cooc_count = 0;                                                % Co-occurence count init
            lag =0;
            for oc_j = 1:length(occur_ij)
                for oc_k =1:length(occur_ik)
                    oPercent = (length(cands_ij{occur_ij(oc_j)})*25)/100;
                    if(length(intersect(cands_ij{occur_ij(oc_j)}, cands_ik{occur_ik(oc_k)}))>=oPercent)
                        cooc_count = cooc_count+1;                         % Co-occurence count incremented
                        sh1 = cands_ij{occur_ij(oc_j)};
                        sh2 = cands_ik{occur_ik(oc_k)};
                        if(sh1(1)-sh2(1)>0)
                            lag=lag+1;
                        end
                        break;
                    end
                end
            end
            co_occurences(j,k) = cooc_count/length(occur_ij);         % Co-occurences betn pair j and k of sub i
            lags(j,k)           = lag/length(occur_ij); 
         end
     end
      save(fullfile(dirr,['co-occur_',num2str(i),'.mat']),'co_occurences');
      save(fullfile(dirr1,['lag_',num2str(i),'.mat']),'lags');
end