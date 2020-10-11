function escore = dEMD(x,y)
prevemd = 0;
escore = 0;
for i = 1:length(x)
    emd = (x(i)+prevemd)-y(i);
    escore = escore+abs(emd); % will tlak about the abs
    prevemd =  abs(emd);     
end
%escore = sum(abs(emd(:)));
end
%% Preprocessing required before calculating EMD
% x = x-min(x);
% x = x/sum(x);
% y = y-min(y);
% y = y/sum(y);
%%
