function[y,idx] = maxk(A,k)
[AA,ids]=sort((A),'descend');
idx = ids(1:k);
y=AA(1:k);
end