function [cBIC, RSS] = clusterBIC(D, IDX)

[nobs, k] = size(D);

for ii = 1:k
    Din(ii) = sum(D(IDX == ii,ii).^2); % sum of squared distances
end
RSS = sum(Din);

cBIC = nobs*log(RSS/nobs) + k*log(nobs);
