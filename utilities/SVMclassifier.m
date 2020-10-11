function r = SVMclassifier(X,labels)
% Intended to classify based on the length of the statelet  
% X is the input data #subject X #predictors 
% for length probably 1 avg length per subject or 1081 (per pair)values per subjects
% Return the classification loss 
% Accuracy is 1 - loss

M = fitcsvm(X, labels);    % trainng the model
CM = crossval(M);          % Cross validation/ testing
r = abs(1-kfoldLoss(CM));  % k-fold loss 
end