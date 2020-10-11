function stvec = compute_statprob(TM)

[V,D] = eig(TM');
stIND = find(abs(diag(D) -1) < 1000*eps); % find the eigenvector with eigenvalue of 1
if isempty(stIND)
    stvec = [];
else
    stvec = V(:,stIND)/sum(V(:,stIND));
end
