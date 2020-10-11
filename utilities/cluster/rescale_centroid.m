function  C = rescale_centroid(X,IDX)

k =length(unique(IDX));
C = zeros(k,size(X,2));
for ii = 1:k
    C(ii,:) = mean(X(IDX == ii,:));
end