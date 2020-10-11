function c = initShape()
%% To initialize structure 'Shape'

measBuff = struct('extrapolated',     [], ...
                  'real_length',      [],...
                  'length',           0,...
                  'maxconnectivity',  0,...
                  'meanconnectivity', 0,...
                  'pair',             0,...
                  'subjectclass',     0,...
                  'probabilitydensity',0); 
c = measBuff;
end