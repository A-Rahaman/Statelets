function[y,idx] = mink(A,k)
[AA,ids]=sort(A);
idx = ids(1:k);
y=AA(1:k);
end