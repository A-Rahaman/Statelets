function fr = combinationUtil(arr, n, r,index, data, i)
%{
Create all possible subsequences (with size r) of a given time series (of length n) without any repitation of items   
arr   = Array of time points/window indices. For example,(1:136)
data  = Temporary array to store current subsequence 
index = Current index in data. Count number of items in a subsequence 
r     = Size of a subsequence to be printed 
Input: arr = [10 20 30 40 50], r = 3, n = 5, index = 1, i = 1
Output: 
10 20 30 
10 20 40 
10 20 50 
10 30 40 
10 30 50 
10 40 50 
20 30 40 
20 30 50 
20 40 50 
30 40 50 
Dated: 6th Feb, 2019
by munna
%}
% Current subsequence is ready! Store it to a global repository.
% We need a global one since the function is recursive 
global subseq;
if (index-1) == r % -1 is importnat since it is required to get index r+1 time incremented before printing r items
    for j = 1:r
        fprintf("%d ", data(j));
    end
    % Store the subsequences in a global cell array according to their size  
    subseq{end+1,length(data)} = data;
    fprintf("\n");
    return;
end

% When no more elements to put in data
if i > n % starts indexing from 1. So no ">="
    return;
end

% current is included, put next at next
% location
data(index) = arr(i);
combinationUtil(arr, n, r, index + 1,data, i + 1);

% current is excluded, replace it with
% next (Note that i+1 is passed, but
%index is not changed)
combinationUtil(arr, n, r, index, data, i + 1);
end

