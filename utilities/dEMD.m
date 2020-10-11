function escore = dEMD(x,y)
prevemd = 0;
emd = zeros(length(x),1);
for i = 1:length(x)
    emd(i) = (x(i)+prevemd)-y(i);
    prevemd =  emd(i);     
end
escore = sum(abs(emd(:))); % [ https://wikimedia.org/api/rest_v1/media/math/render/svg/e15df1f0ad325482b1a23321c52a9f9ca06c09b1 ]
end

