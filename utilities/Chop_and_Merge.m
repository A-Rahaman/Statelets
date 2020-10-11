function [mergedSig,ids,remLengths] = Chop_and_Merge(data,shapelets,candidates)
mergedSig = [];
MAXSIMS = zeros(size(shapelets,1),1);
lengths = zeros(size(shapelets,1),1);
for sh = 1:size(shapelets,1)
tid = shapelets{sh,1}; 
lengths(sh) = tid;
MAXSIMS(sh) = shapelets{sh,2}(tid); % Similarity Values
end

ids = find(MAXSIMS>=mean(MAXSIMS));
remLengths = lengths(ids);

for i = 1:size(ids,1)
%shapelets_i_occur = shapelets{i,3};
data_i            = data(:,ids(i));
cands_i           = candidates{shapelets{ids(i),1}};
        pp = [];
        %pp = data_i(cands_i{shapelets_i_occur(j)});
        pp = data_i(cands_i{shapelets{ids(i),4}});
        mergedSig = horzcat(mergedSig,pp');
end
end