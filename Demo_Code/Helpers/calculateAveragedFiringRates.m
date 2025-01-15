function [averagedFiringRates, nEntries] = calculateAveragedFiringRates(ratemap, zonetrial)
    % Can calculate the mean of any matrix A of size mxn based on index vector B mx1, and produce the
    % resulting matrix of size max(B)xn.
    
    % Input:
    % - ratemap: Matrix of firing rates (mTrial x nBin)
    % - zonetrial: Vector indicating the times the animal enters a specific zone (mTrial x 1)
    % Output:
    % - averagedFiringRates: Matrix of averaged firing rates (max(zonetrial) x nBin)

    % Get unique zone entry times
    uniqueZoneEntryTimes = unique(zonetrial);

    % Initialize a matrix to store the averaged firing rates
    averagedFiringRates = zeros(length(uniqueZoneEntryTimes), size(ratemap, 2));
    nEntries = zeros(length(uniqueZoneEntryTimes), 1);
    % Loop through each unique zone entry time
    for i = 1:length(uniqueZoneEntryTimes)
        % Extract firing rates for the current zone entry time
        currentZoneTrials = zonetrial == uniqueZoneEntryTimes(i);
        nEntries(i) = sum(currentZoneTrials);
        selectedTrials = ratemap(currentZoneTrials, :);
        
        % Calculate the mean firing rate for the current zone entry time
        averagedFiringRates(i, :) = nanmean(selectedTrials, 1);
    end
end
