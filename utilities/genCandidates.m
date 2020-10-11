function candidates = genCandidates(data,max_len,min_len)
% generates all possible consecutive subsequnces between min_len to max_len of data
% sliding over the time points 1 by one 
% data = is a vector [1:136].The values we want subsequences of   
range = max_len-min_len+1;
cand = {};
l = max_len;
while l>=min_len
    for k = 1:(length(data)-l)
        %fprintf("in\n")
        cand{end+1,1} = data(k:k+l-1); % index start from 1 so, K+L creates a length of L+1 where we're looking at size L.
        %So, we require k+l-1
        %fprintf("%s",cand{1,1});
        %count = count+1;
        %cand{end+1,2} = labels(k:k+l);
    end
    candidates{:,range} = cand;
    cand={};
    l=l-1;
    range = range-1;
end
%candidates = cand;
end