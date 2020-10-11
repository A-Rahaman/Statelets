SUB_shapelet_Len = zeros(314,5);
for i =1:size(SUB_wise_shapelets_uLen)
    for j=1:size(SUB_wise_shapelets_uLen,2)
    SUB_shapelet_Len(i,j) = length(SUB_wise_shapelets_uLen{i,j});
    end
end
median(SUB_shapelet_Len(:))
hist(SUB_shapelet_Len(:))
mean(SUB_shapelet_Len(:))