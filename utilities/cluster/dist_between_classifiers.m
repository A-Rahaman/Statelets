function d = dist_between_classifiers(ID1, ID2)

s = 0;
iter = 0;
% loop through all i,j possible pairs (too large to vectorize)
for ii = 1:length(ID1)
    for jj = ii+1:length(ID2)
        iter = iter+1;
        V1 = ID1(ii) == ID1(jj);
        V2 = ID2(ii) == ID2(jj);
        s = s + [(V1+V2) == 1 ];
    end
end
d = s/iter; 

%k = length(unique([ID1(:); ID2(:)]));
% %% Determine the normalizing probability 
% p1 = zeros(1,k);
% p2 = zeros(1,k);
% for ii = 1:k
%    p1(ii) = sum(ID1 == ii)/length(ID1);
%    p2(ii) = sum(ID2 == ii)/length(ID2);
% end
% 
% %prob of V1/V2 being 1
% p1_1 = sum(p1.*p1);
% p2_1 = sum(p2.*p2);
% 
% %prob of their sum being 1
% pchance = p1_1*(1-p2_1) + p2_1*(1-p1_1);
% 
% s_n = pchance;
% s_raw = s;


% pairs = nchoosek(1:length(ID1),2);
% INDICATOR1 = ID1(pairs(:,1)) == ID1(pairs(:,2));
% INDICATOR2 = ID2(pairs(:,1)) == ID2(pairs(:,2));
% 
% s = sum([INDICATOR1 + INDICATOR2] == 1); 