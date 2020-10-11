function lminimas = findLocalMinimas(slices,cand_j,cand_jj,sigmoid)
% Unlike having an average across all the distance, We focuse on the best matches in different part of the signal  
% Our similarity measurement function returns differences instead
% Go to each slice of the 'slices'. Pick one minima per slice 
% each slice corresponds to a different part of the time series but a specific length 

itr = 0;
lminimas = zeros(length(slices),1);
for sl = 1:length(slices)
    tempSlice = slices{sl};
    lSigmas = [];
    lc = 1;
    for cn = (itr+1):length(cand_jj)
        if(length(intersect(tempSlice,cand_j{cand_jj(cn)}))>=length(tempSlice)/2)
            lSigmas(lc) = sigmoid(cn);
            itr = cn;
            lc = lc+1;
        end
    end
    %fprintf("Length of local sigmas %d\n",length(lSigmas));
    if(~isempty(lSigmas))
    lminimas(sl) = min(lSigmas);
    end
end
end