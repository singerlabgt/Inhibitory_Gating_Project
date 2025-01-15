clear; close all;
%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

% load Neuron information (bad, celltype, sessioninfo, etc. for each day)
load(fullfile(maindir, "Demo_Data/", "NeuroInfo.mat"));
% load neurons with significant increase or decrease in AZRZ (for PV
% selection)
load(fullfile(maindir, "Demo_Data/", "allsess_shuffledSig_distance2RZ.mat"));
% load neuron firing rates per day (for PV selection)
load(fullfile(maindir, "Demo_Data/", "allsess_raw_vs_residuals_distance2RZ.mat"));
% load ratemap for each day 
outmapdir = fullfile(maindir, "Demo_Data/", "singlesess_ratemap_distance2RZ/");
%%
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.animals);
allindex = allindex(animal_idx,:); %filter based on animals to include 
allindex = allindex(animal_idx,:); %filter based on animals to include 

[sessions, sessID] = unique(allindex(:,[1:2, 7]),'rows'); %session = [animalID, recording date, novelty day]

%%
load_data = 1; % spatial information calculation can take some time. Load existing data by specifying load_data = 1;
if load_data == 0
    correct_only = 1;
    decrease_only = 0; % select neurons hvagin significant decrease at goals
    increase_only = 0;
    separatezone = 0;
    scale = 1;
    spatialinfo_analysis = 1; % in addition to firing rate traces
    ratemap_analysis = 1; % in addition to firing rate traces
    
    allday_ratemap = zeros(0, 65);
    allday_trialidx = zeros(0, 1);
    allday_novenv = zeros(0, 1);
    allday_novelexposureday = zeros(0, 1);
    allday_iZone = zeros(0, 1);
    allday_zonetrial = zeros(0, 1); % trial num on each zone
    allday_spatialinfo = zeros(0, 6); % env, trial block, novelexposure day, cell, SI, SI shuffled
    allday_ratemapcorr = zeros(0, 5); % env, trial block, novelexposure day, cell, ratemapCorr
    trialblocksize = [25 40];
    
    for iDay = 1:length(NeuronInfo)
        disp(['processing day' num2str(iDay)])
        fn = getlatestfile_with_string(outmapdir, NeuronInfo(iDay).SessionName);
        if length(fn) < 1
            continue
        end
        load(fullfile(outmapdir, fn))
        % identify the bad cells from NeuronInfo
        notbadidx = find(ismember(NeuronInfo(iDay).TSIDtotal, NeuronInfo(iDay).TSID_notbad));
        notbadCluID = NeuronInfo(iDay).CluIDtotal(notbadidx); % check if match allsess
        % identify session related cells in allsess
        sessionname = NeuronInfo(iDay).SessionName;
        sessionname_split = strsplit(sessionname, '_');
        animalnum = str2num(sessionname_split{1}(2:end));
        day = str2num(sessionname_split{2});
        index = allindex(allindex(:, 1) == animalnum & allindex(:, 2) == day, :);
        env = unique(index(:, 6));
        
        % only include WT animals
        if ~ismember(animalnum, params.WTmice)
            continue
        end
    
        currsess_idx = find(allsess_sessinfo_sh(:, 1)==animalnum & allsess_sessinfo_sh(:, 2)==day);
        if isempty(currsess_idx)
            continue
        end
        if min(allsess_unitID_sh(currsess_idx)' == notbadCluID) == 0 
            disp(['NeuronInfo notbadcluID and allsess_unitID_sh does not match on iDay ', num2str(iDay)])
        end
        % select neurons:
        % not bad, nov decreasing -15-20 deg, low firing, pyramidal cells
        PYR_units = find(ismember(allsess_unitType_sh(currsess_idx), 'Pyramidal Cell'));
        
        if decrease_only == 1
            AZRZ_decrease_units = find(sum(allsess_sigDecrease_fam(currsess_idx, 11:16), 2) > 0);
            % unit2incl =  intersect(notbadidx, intersect(AZRZ_decrease_units, PYR_units));
            unit2incl =  intersect(AZRZ_decrease_units, PYR_units);
        elseif increase_only == 1
            AZRZ_increase_units = find(sum(allsess_sigIncrease_fam(currsess_idx, 11:16), 2) > 0);
            unit2incl =  intersect(AZRZ_increase_units, PYR_units);
        else
            unit2incl = PYR_units;
        end
        if isempty(unit2incl)
            continue
        end
        % select trials:
        if correct_only == 1
            select_trials = find(outmap.bin2.nov.labels(:, 2) == 1 & outmap.bin2.nov.labels(:, 3) == 1); % RZ-centered and receive reward
        else
            select_trials = find(outmap.bin2.nov.labels(:, 2) == 1 & outmap.bin2.nov.labels(:, 3) == 0); % only select RZ_centered
        end
        if isempty(select_trials) || length(select_trials) < 5
            continue
        end
        % for each unit (including good or bad neurons of all types), regress out speed and lick rate
        ratemap_nov_residual = zeros(size(outmap.bin2.nov.ratemap));
        x1 = reshape([outmap.bin2.fam.smoothSpeed; outmap.bin2.nov.smoothSpeed], [], 1);
        x2 = reshape([outmap.bin2.fam.smoothLickrate; outmap.bin2.nov.smoothLickrate], [], 1);
        X = [ones(size(x1)) x1 x2 x1.*x2];
        for iCell = 1:size(outmap.bin2.fam.ratemap, 3)
            f_rate = outmap.bin2.fam.ratemap(:, :, iCell); 
            n_rate = outmap.bin2.nov.ratemap(:, :, iCell);
            combined_rate = [f_rate; n_rate];
            y = reshape(combined_rate, [], 1);
            [b, bint, r, rint, stats] = regress(y, X);
            tempRes = reshape(r, size(combined_rate));
            famRes = tempRes(1:size(f_rate, 1), :);
            novRes = tempRes(size(f_rate, 1)+1:end, :);
            ratemap_nov_residual(:, :, iCell) = novRes(:, :);
        end
        ratemap_nov_residual = ratemap_nov_residual(select_trials, :, notbadidx(unit2incl));
        %% calculate spatial info
        if spatialinfo_analysis == 1 && ~isempty(select_trials)
            rawCounts = outmap.bin2.nov.rawCount(select_trials, :, notbadidx(unit2incl));
            rawOccups = outmap.bin2.nov.rawOccup(select_trials, :); 
            smoothCounts = outmap.bin2.nov.smoothCount(select_trials, :, notbadidx(unit2incl));
            smoothOccups = outmap.bin2.nov.smoothOccup(select_trials, :);
            novelexposureday = unique(allsess_sessinfo_sh(currsess_idx, 3));
        
            if ismember(2, env)
                trialblock = 25; track = 2;
            elseif ismember(3, env)
                trialblock = 40; track = 3;
            end
    
            for iBlock = 1:floor(length(select_trials)/trialblock)
                
    
                select_trials_tmp = trialblock*(iBlock-1)+1 : trialblock*(iBlock);
                
            
                for iCell = 1:length(notbadidx(unit2incl))
                    temp_spatialinfo = zeros(1, 6);
                    rawCount = rawCounts(select_trials_tmp, :, iCell);
                    rawOccup = rawOccups(select_trials_tmp, :);
                    smoothCount = smoothCounts(select_trials_tmp, :, iCell);
                    smoothOccup = smoothOccups(select_trials_tmp, :);
                    testtest = reshape(smoothOccup, [], 1); 
                    if sum(isnan(testtest)) > 0
                        disp(['some negative occupancy'])
                    end
                    if size(smoothCount, 1) > 1
                        smoothCountSum = nansum(smoothCount); % average across trials, output 1xnBin
                        smoothOccupSum = nansum(smoothOccup);
                    else
                        smoothCountSum = smoothCount;
                        smoothOccupSum = smoothOccup;
                    end
                    ratemap_original = smoothCountSum ./ smoothOccupSum;
                    
                    %% abby 8/20/24
                    normOcc = smoothOccupSum./nansum(smoothOccupSum);
                    binInfo2 = zeros(1,length(ratemap_original));
                    for binIdx = 1:length(ratemap_original)
                        binFR = smoothCountSum(binIdx)/(nansum(smoothCountSum.*normOcc));
                        binOcc = smoothOccupSum(binIdx)/nansum(smoothOccupSum);
                        binInfo2(binIdx) = binOcc*binFR*log2(binFR);
                    end
                    spatialinfo_original = nansum(binInfo2);
                    % method 3. % NOTE BY XZ 08212024: same calculation as
                    % previous one but in different syntax
                    posPDF = smoothOccupSum;
                    map = smoothCountSum;
                    posPDF = posPDF./nansum(nansum(posPDF));
                    meanrate = nansum(nansum(map.*posPDF));
                    normMap = map./meanrate;
                    temp = normMap.*log2(normMap);
                    spatialinfo_method3 = nansum(nansum(posPDF.*temp));

                    %% end abby 8/20/24


                    % %quantify spaital info from original map
                    % binInfo = zeros(1,length(ratemap_original));
                    % %loop through each bin
                    % for binIdx = 1:length(ratemap_original)
                    %     binFR = ratemap_original(binIdx) ./nanmean(ratemap_original); %FRi/meanFR
                    %     binOcc = smoothOccupSum(binIdx) ./ nansum(smoothOccupSum); %occ/totaltime
                    %     binInfo(binIdx) = binOcc*binFR*log2(binFR);
                    % end
                    % spatialinfo_original = nansum(binInfo);
                    % if spatialinfo_original < 0
                    %     A  = 1;
                    % end

                    %quantify spaital info from shuffled map    
                    numShuffles = 1000;
                    ratemap_shuffled = zeros(numShuffles, length(rawCount));
                    spatialinfo_shuffled = zeros(numShuffles, 1);
                    for iShuff = 1:numShuffles
                        shuffledSpikes = zeros(size(rawCount));
                        shuffledOccup = zeros(size(rawOccup));
                        
                        for ii = 1:size(shuffledSpikes,1) %shuffle both spike counts and occupancy map
                            tempCount = rawCount(ii, :);
                            shuffledSpikes(ii,:) = tempCount(randperm(length(rawCount(ii,:))));
                            tempOccup = rawOccup(ii, :);
                            shuffledOccup(ii,:) = tempOccup(randperm(length(rawOccup(ii,:))));
                        end
                        shuffledCountSum = nansum(smoothGaussianMultiple(shuffledSpikes,2));
                        shuffledOccupSum = nansum(smoothGaussianMultiple(shuffledOccup,2));
                        ratemap_shuffled = shuffledCountSum ./ shuffledOccupSum;
            

                        %% xz 8/21/24
                        normOcc = shuffledOccupSum./nansum(shuffledOccupSum);
                        binInfo2 = zeros(1,length(ratemap_shuffled));
                        for binIdx = 1:length(ratemap_shuffled)
                            binFR = shuffledCountSum(binIdx)/(nansum(shuffledCountSum.*normOcc));
                            binOcc = shuffledOccupSum(binIdx)/nansum(shuffledOccupSum);
                            binInfo2(binIdx) = binOcc*binFR*log2(binFR);
                        end
                        spatialinfo_shuffled(iShuff) = nansum(binInfo2);
                        
    
                        %% end xz 8/21/24
                        



                        % %loop through each bin
                        % for binIdx = 1:length(ratemap_shuffled)
                        %     binFR = ratemap_shuffled(binIdx) ./nanmean(ratemap_shuffled); %FRi/meanFR
                        %     binOcc = smoothOccupSum(binIdx) ./ nansum(smoothOccupSum); %occ/totaltime
                        %     binInfo(binIdx) = binOcc*binFR*log2(binFR);
                        % end
                        % spatialinfo_shuffled(iShuff) = nansum(binInfo);
                    end
                    temp_spatialinfo(:, 5) = spatialinfo_original; % scalar for each cell 
                    temp_spatialinfo(:, 6) = prctile(spatialinfo_shuffled, 95); % threshold_95th
                    temp_spatialinfo(:, 7) = animalnum; % added by xz 092524 to LMM
                    temp_spatialinfo(:, 4) = iCell;
                    temp_spatialinfo(:, 3) = novelexposureday;
                    temp_spatialinfo(:, 2) = iBlock;
                    temp_spatialinfo(:, 1) = track;
                    allday_spatialinfo = [allday_spatialinfo; temp_spatialinfo];
                end
            end
        end
        %% calculate rate map correlation
        if ratemap_analysis == 1 && ~isempty(select_trials)
            ratemap_all = outmap.bin2.nov.ratemap(select_trials, :, notbadidx(unit2incl));
            novelexposureday = unique(allsess_sessinfo_sh(currsess_idx, 3));
            if ismember(2, env)
                trialblock = 25; track = 2;
            elseif ismember(3, env)
                trialblock = 40; track = 3;
            end
            for iBlock = 1:floor(length(select_trials)/trialblock)
                select_trials_tmp = trialblock*(iBlock-1)+1 : trialblock*(iBlock);
                for iCell = 1:length(notbadidx(unit2incl))
                    temp_ratemapcorr = zeros(1, 5);
                    ratemap = ratemap_all(select_trials_tmp, :, iCell);
                    temp_ratemapcorr(:, 5) = nanmean(arrayfun( @(x) diag(corrcoef(ratemap(x,:), nanmean(ratemap,1)),1), 1:size(ratemap,1)));
                    temp_ratemapcorr(:, 6) = animalnum; % added by xz 092524 to LMM
                    temp_ratemapcorr(:, 4) = iCell;
                    temp_ratemapcorr(:, 3) = novelexposureday;
                    temp_ratemapcorr(:, 2) = iBlock;
                    temp_ratemapcorr(:, 1) = track;
                    allday_ratemapcorr = [allday_ratemapcorr; temp_ratemapcorr];
                end
            end
        end
        %% append data to all day table
        
        % averaging all cells
        ratemap_nov_residual = nanmean(ratemap_nov_residual, 3);
    
        allday_ratemap = [allday_ratemap; ratemap_nov_residual];
        allday_trialidx = [allday_trialidx; (1:size(ratemap_nov_residual, 1))']; 
        allday_novelexposureday = [allday_novelexposureday; repelem(unique(allsess_sessinfo_sh(currsess_idx, 3)), size(ratemap_nov_residual, 1), 1)];
        % get nov env
        index = allindex(allindex(:, 1) == animalnum & allindex(:, 2) == day, :);
        env = unique(index(:, 6));
        if ismember(2, env)
            allday_novenv = [allday_novenv; repelem(2, size(ratemap_nov_residual, 1), 1)];
        elseif ismember(3, env)
            allday_novenv = [allday_novenv; repelem(3, size(ratemap_nov_residual, 1), 1)];
        end
    
        zones = outmap.bin2.nov.labels(select_trials, 4);
        allday_iZone = [allday_iZone; zones];
        zonetrial = zeros(size(zones));
        for iZone = 1:3
            iZone_trialidx = find(zones == iZone);
            zonetrial(iZone_trialidx) = 1:length(iZone_trialidx);
        end
        allday_zonetrial = [allday_zonetrial; zonetrial];
    
    end
    save(fullfile(maindir, 'Demo_Data', 'allday_spatialinfo'), "allday_spatialinfo")
    save(fullfile(maindir, 'Demo_Data', 'allday_ratemapcorr'), "allday_ratemapcorr")
else
    load(fullfile(maindir, 'Demo_Data', 'allday_spatialinfo_092524.mat'))
    load(fullfile(maindir, 'Demo_Data', 'allday_ratemapcorr_092524.mat'))
end
%% plot trialblock spatial info  
% plot track B
fig = figure;
fig.Position = [225 432 519 192]; box off; t=tiledlayout(1, 2);
SI_edges = min(allday_spatialinfo(:, 5)):0.1:max(allday_spatialinfo(:, 5));
numColors = 4;
reds = [linspace(1, 0.5, numColors)', zeros(numColors, 1), zeros(numColors, 1)];
PConly = 0;
disp('Spatial info')
for ii = 1
    ax = nexttile; hold on
    for iBlock = 1:4 % code should be optimized later
        if PConly == 1
            temp = allday_spatialinfo(allday_spatialinfo(:, 1) == 2 & allday_spatialinfo(:, 2) == iBlock ...
                & allday_spatialinfo(:, 3) == ii & allday_spatialinfo(:, 5) > allday_spatialinfo(:, 6), 5);
        else
            temp = allday_spatialinfo(allday_spatialinfo(:, 1) == 2 & allday_spatialinfo(:, 2) == iBlock ...
                & allday_spatialinfo(:, 3) == ii, 5);
        end
        % disp(['number of selected cells on day ' num2str(ii) ', iBlock ' num2str(iBlock) ' is: ' num2str(length(temp(~isnan(temp))))])
        disp(['number of selected cells on day ' num2str(ii) ', iBlock ' num2str(iBlock) ' is: ' num2str(length(temp(~isnan(temp)))) ', percentile = ' ...
            num2str(prctile(temp, [0, 25, 50, 75, 100]))])
        temp = histcounts(temp(~isnan(temp)),SI_edges,'Normalization','probability');
        plot(ax, SI_edges(2:end),cumsum(temp),'Color', reds(iBlock,:),'LineWidth',2);
        alpha(ax, 1/iBlock)
        ylim(ax, [0 1])
    end
    title(['env=2, ' 'day=' num2str(ii) ])
    ylabel('Cumulative prop. of pyr cells')
    xlabel('Spatial info. (bits/spike)')
    set(ax, 'TickDir', 'out');
end



%% plot trialblock ratemap correlation
% plot track B
box off; 
disp('Ratemap correlation')
RCorr_edges = -0.2 : 0.05 : 1;
numColors = 4;
reds = [linspace(1, 0.5, numColors)', zeros(numColors, 1), zeros(numColors, 1)];
for ii = 1
    ax = nexttile; hold on
    for iBlock = 1:4 % code should be optimized later
        if PConly == 1
            temp = allday_ratemapcorr(allday_ratemapcorr(:, 1) == 2 & allday_ratemapcorr(:, 2) == iBlock ...
                & allday_ratemapcorr(:, 3) == ii & allday_spatialinfo(:, 5) > allday_spatialinfo(:, 6), 5);
        else
            temp = allday_ratemapcorr(allday_ratemapcorr(:, 1) == 2 & allday_ratemapcorr(:, 2) == iBlock ...
            & allday_ratemapcorr(:, 3) == ii, 5);
        end
        disp(['number of selected cells on day ' num2str(ii) ', iBlock ' num2str(iBlock) ' is: ' num2str(length(temp(~isnan(temp)))) ', percentile = ' ...
            num2str(prctile(temp, [0, 25, 50, 75, 100]))])
        temp = histcounts(temp(~isnan(temp)),RCorr_edges,'Normalization','probability');
        plot(ax, RCorr_edges(2:end),cumsum(temp),'Color', reds(iBlock,:),'LineWidth',2);
        alpha(ax, 1/iBlock)
        ylim(ax, [0 1])
    end
    title(['env=2, ' 'day=' num2str(ii) ])
    ylabel('Cumulative prop. of pyr cells')
    xlabel('Ratemap Corr')
    set(ax, 'TickDir', 'out');
end

makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure07_D';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')