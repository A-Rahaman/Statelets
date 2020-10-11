plot(allshapes_s{1})
sims = simscore(1,:);
[sortedSIM,sdx] = sort(sims);
sx = setdiff(1:2978,1);
SIMidx = sx(sdx);
figure();
%subplot(1,10,1);

for i=1:10
    subplot(10,1,i)
    plot(allshapes_s{SIMidx(i)});
    hold on
end
c = 1;
for i=2:length(allshapes_s)
sigmoid(c) = weightedDrift1(allshapes_s{1},allshapes_s{i}); 
c = c+1;
end
[sSigma,sigIdx] = sort(sigmoid,'descend');
sigx = setdiff(1:2978,1);
SIGidx = sigx(sigIdx);

figure();
for i=1:10
    subplot(10,1,i)
    plot(allshapes_s{SIGidx(i)});
    hold on
end