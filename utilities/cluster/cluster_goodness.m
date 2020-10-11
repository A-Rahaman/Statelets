function [SSE, R] = cluster_goodness(D, IDX)


[nobs, k] = size(D);

for ii = 1:k
    Din(ii) = sum(D(IDX == ii,ii).^2); % dispersion in cluster 
    Dout(ii) = sum(D(IDX ~= ii,ii).^2); % sum of squared distances
end
SSE = sum(Din);
R = mean(Din./Dout);


