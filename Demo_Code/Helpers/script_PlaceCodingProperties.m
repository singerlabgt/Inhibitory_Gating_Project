%% script_PlaceCodingProperties
% This script takes the raw binned firing rate maps and shuffled rate maps
% to calculate spatial information and identify "spatially modulated"
% units. To determine spatially modualted cells, used 95th percentile or
% above the spatial information distribution from shuffled data. Takes a
% long time to run all 1000 shuffled data. 

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Nuri\Code\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Wide Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','WS Interneuron','Pyramidal Cell'}; %names used in our manuscript
for ct = 1:length(celltypes)
    cellT{ct} = find(...
        strcmp(cell_metrics.putativeCellType, celltypes{ct})...
        & ~ismember(1:length(cell_metrics.cellID), cell_metrics.tags.Bad));
end

%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);
sessions = unique(allindex(:,1:2),'rows');

%% create place coding output structure 
numShuffles = 1000;
environments = params.environments;
trialTypes = {'distance2RZ','laps'}; 

for iSess = 1:size(sessions,1)
% for iSess = 1:length(sessions)
    animal = sessions(iSess,1);
    recDay = sessions(iSess,2);
    index = allindex(allindex(:,1) == animal & allindex(:,2) == recDay, :);
    novelDay = unique(index(:,7));
    sessStr = [params.iden num2str(animal) '_' num2str(recDay)];

    for iType = 1:length(trialTypes) %lap-by-lap or distance-to-reward trials
        tt = trialTypes{iType};
        
        %load binned firing rate maps with behavioral info
        fdir = fullfile(dirs.data2load, ['singlesess_ratemap_' tt]);
        if iType > 1
            fname = getlatestfile_with_string(fdir,...
                [sessStr '_createTrialByTrialRateMaps_Day ' num2str(novelDay)]);
        else
            fname = getlatestfile_with_string(fdir,...
                [sessStr '_createTrialByTrialDist2RewardRateMaps_Day ' num2str(novelDay)]);
        end
        load(fullfile(fdir, fname));

        binOptions = outmap.binOptions;
        for iE = 1:length(environments)
            env = environments{iE};

            for edg = 1:length(binOptions)
                binSize = binOptions(edg);

                %dynamic binsize name
                dynBinsize = ['bin' num2str(binSize)];

                %way to deal with no-novel sessions
                if ~isempty(outmap.(dynBinsize).(env).labels)

                    %pick out only RZ-centered trals for creating distance2RZ ratemaps
                    if strcmp(tt, 'Distance2RZ')
                        trials_rz_centered = find(outmap.(dynBinsize).(env).labels(:,2) == 1);
                    else %otherwise include all trials for consideration
                        trials_rz_centered = (1:size(outmap.(dynBinsize).(env).rawCount,1))';
                    end

                    %pick out trials during running
                    smoothSpeed = outmap.(dynBinsize).(env).smoothSpeed;
                    trials_below_speed = smoothSpeed < params.speedTh; %during movement

                    if ~isempty(outmap.(dynBinsize).(env).ratemap) && ~isempty(trials_rz_centered)
                        for iUnit = 1:length(outmap.unitID)
                            rawCount = outmap.(dynBinsize).(env).rawCount(:,:,iUnit);
                            rawCount(trials_below_speed) = 0; %remove spikes in low speed bins 
                            rawCount = rawCount(trials_rz_centered,:);
                            %smoothe AFTER removing spikes that occurred
                            %durng low speeds
                            smoothCount = smoothGaussianMultiple(rawCount, 0.5);
                            if size(outmap.(dynBinsize).(env).rawOccup,3) > 1
                                rawOccup = outmap.(dynBinsize).(env).rawOccup(trials_rz_centered,:,iUnit);
                                smoothOccup = outmap.(dynBinsize).(env).smoothOccup(trials_rz_centered,:,iUnit);
                            else
                                rawOccup = outmap.(dynBinsize).(env).rawOccup(trials_rz_centered,:);
                                smoothOccup = outmap.(dynBinsize).(env).smoothOccup(trials_rz_centered,:);
                            end

                            %in rare cases where only one trial passes speed threshold, do not sum across one trial; NJ 08/20/22
                            if length(trials_rz_centered) > 1
                                smoothCountSum = nansum(smoothCount);
                                smoothOccupSum = nansum(smoothOccup);
                            else
                                smoothCountSum = smoothCount;
                                smoothOccupSum = smoothOccup;
                            end
                            ratemap_original = smoothCountSum ./ smoothOccupSum;

                            % %quantify spaital info from original map
                            % binInfo = zeros(1,length(ratemap_original));
                            % %loop through each bin
                            % for binIdx = 1:length(ratemap_original)
                            %     binFR = ratemap_original(binIdx) ./nanmean(ratemap_original); %FRi/meanFR
                            %     binOcc = smoothOccupSum(binIdx) ./ nansum(smoothOccupSum); %occ/totaltime
                            %     binInfo(binIdx) = binOcc*binFR*log2(binFR);
                            % end
                            % spatialinfo_original = nansum(binInfo);
                            %% added xz 08212024
                            normOcc = smoothOccupSum./nansum(smoothOccupSum);
                            binInfo2 = zeros(1,length(ratemap_original));
                            for binIdx = 1:length(ratemap_original)
                                binFR = smoothCountSum(binIdx)/(nansum(smoothCountSum.*normOcc));
                                binOcc = smoothOccupSum(binIdx)/nansum(smoothOccupSum);
                                binInfo2(binIdx) = binOcc*binFR*log2(binFR);
                            end
                            spatialinfo_original = nansum(binInfo2);
                            %% end xz 08212024

                            %quantify spaital info from shuffled map                            
                            ratemap_shuffled = zeros(numShuffles, length(rawCount));
                            spatialinfo_shuffled = zeros(numShuffles, 1);
                            for iShuff = 1:numShuffles
                                shuffledSpikes = zeros(size(rawCount));
                                shuffledOccup = zeros(size(rawOccup));
                                for iSess = 1:size(shuffledSpikes,1) %shuffle both spike counts and occupancy map
                                    tempCount = rawCount(iSess, :);
                                    shuffledSpikes(iSess,:) = tempCount(randperm(length(rawCount(iSess,:))));
                                    tempOccup = rawOccup(iSess, :);
                                    shuffledOccup(iSess,:) = tempOccup(randperm(length(rawOccup(iSess,:))));
                                end
                                %% xz 8/21/24   
                                shuffledCountSum = nansum(smoothGaussianMultiple(shuffledSpikes,2));
                                shuffledOccupSum = nansum(smoothGaussianMultiple(shuffledOccup,2));
                                ratemap_shuffled = shuffledCountSum ./ shuffledOccupSum;
                                normOcc = shuffledOccupSum./nansum(shuffledOccupSum);
                                binInfo2 = zeros(1,length(ratemap_shuffled));
                                for binIdx = 1:length(ratemap_shuffled)
                                    binFR = shuffledCountSum(binIdx)/(nansum(shuffledCountSum.*normOcc));
                                    binOcc = shuffledOccupSum(binIdx)/nansum(shuffledOccupSum);
                                    binInfo2(binIdx) = binOcc*binFR*log2(binFR);
                                end
                                spatialinfo_shuffled(iShuff) = nansum(binInfo2);
                                %% end xz 8/21/24
                                % shuffledCount = nansum(smoothGaussianMultiple(shuffledSpikes,2));
                                % shuffledOccup = nansum(smoothGaussianMultiple(shuffledOccup,2));
                                % shuffledOccupSum = nansum(shuffledOccup);
                                % ratemap_shuffled = shuffledCount ./ shuffledOccup;
                                % 
                                % binInfo = zeros(1,length(ratemap_shuffled));
                                % %loop through each bin
                                % for binIdx = 1:length(ratemap_shuffled)
                                %     binFR = ratemap_shuffled(binIdx) ./nanmean(ratemap_shuffled); %FRi/meanFR
                                %     binOcc = smoothOccupSum(binIdx) ./ nansum(smoothOccupSum); %occ/totaltime
                                %     binInfo(binIdx) = binOcc*binFR*log2(binFR);
                                % end
                                % spatialinfo_shuffled(iShuff) = nansum(binInfo);
                            end
                            %export spatial info for place cell criteria
                            outmap.(dynBinsize).(env).spatialinfo_shuffled(:,iUnit) = spatialinfo_shuffled;
                            outmap.(dynBinsize).(env).spatialinfo_original(iUnit) = spatialinfo_original;
                            outmap.(dynBinsize).(env).ratemap_original(iUnit,:) = ratemap_original;
                            %export peak firing rate for place cell criteria
                            outmap.(dynBinsize).(env).peak_firingrate(iUnit) = max(ratemap_original);
                        end
                    end
                end
            end
        end
        savefiledir = fullfile(dirs.saveoutputstruct, 'placecoding', tt);
        if ~isdir(savefiledir)
            mkdir(savefiledir)
        end
        filename = fullfile(savefiledir,[params.iden num2str(animal) '_' num2str(recDay) '_' datestr(now, 'yymmdd')]);
        save(filename, 'outmap')
        clear outmap
    end
end

%collect ratemaps from all sessions
for iSess = 1:length(sessions)
    ratemapdir = fullfile(dirs.saveoutputstruct, 'placecoding','Laps');
    fname = getlatestfile_with_string(ratemapdir,[params.iden num2str(sessions(iSess,1)) '_' num2str(sessions(iSess,2))]);
    s_ratemap(iSess) = load(fullfile(ratemapdir, fname));
end

%% calculate map correlations using place-modulated units 
prev = 1;
for iSess = 1:length(sessions)
    numUnits = length(s_ratemap(iSess).outmap.unitID);
    for iUnit = 1:numUnits
        famLbl = s_ratemap(iSess).outmap.bin5.fam.labels; 
        ratemap = s_ratemap(iSess).outmap.bin5.fam.ratemap(:,:,iUnit);        
        r1 = nanmean(ratemap(famLbl(:,1)==1,:),1); %1st half of session 
        r2 = nanmean(ratemap(famLbl(:,1)==2,:),1); %2nd half of session 
        %correlate between first and second halves of sessions in the same environment
        MapCorrelation.fam_within(prev) = diag(corrcoef(r1,r2),1); 
        %correlate session-averaged ratemap to each trial 
        MapCorrelation.fam_trials(prev) = nanmean(arrayfun( @(x) diag(corrcoef(ratemap(x,:), nanmean(ratemap,1)),1), 1:size(ratemap,1)));
        
        if ~isempty(s_ratemap(iSess).outmap.bin5.nov.ratemap)
            novLbl = s_ratemap(iSess).outmap.bin5.nov.labels; 
            ratemap = s_ratemap(iSess).outmap.bin5.nov.ratemap(:,:,iUnit);
            r3 = nanmean(ratemap(novLbl(:,1)==1,:),1); %1st half of session
            r4 = nanmean(ratemap(novLbl(:,1)==2,:),1); %2nd half of session
            MapCorrelation.nov_within(prev) = diag(corrcoef(r3,r4),1);
            MapCorrelation.between(prev) = diag(corrcoef(r1,r3),1); %correlate first halves of two environment
            MapCorrelation.nov_trials(prev) = nanmean(arrayfun( @(x) diag(corrcoef(ratemap(x,:), nanmean(ratemap,1)),1), 1:size(ratemap,1)));

        else
            MapCorrelation.nov_within(prev) = nan;
            MapCorrelation.between(prev) = nan;
            MapCorrelation.nov_trials(prev) = nan;
        end
        prev = prev + 1;
    end
end

%% save output structure 
save(fullfile(dirs.data2load, ['placecodingout_' datestr(now,'yymmdd'), '.mat']),...
    'sessions','s_ratemap','thetaRZ_fam','thetaRZ_nov','MapCorrelation','-v7.3');
