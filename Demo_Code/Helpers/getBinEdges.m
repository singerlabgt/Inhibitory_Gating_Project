function out = getBinEdges(bins)
%This function takes input of 1xN (or Nx1) array that is the bin edges used
%for histcounts, and outputs (N-1)x2 array that describes the edge bounds
%of each bin. Here the N-1 matches the number of bins 
%NJ created 06.12.2020
out = nan(length(bins)-1,2);

for ii = 1:length(bins)-1
    out(ii,:) = [bins(ii), bins(ii+1)];
end
    