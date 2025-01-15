%% expect scatter plot: Y axis is the decrease index, X axis is trial block (25 for track B, day 1; 40 for track C, day 2). Each dot is a neuron averaging trials within a trial block

clear; close all;
%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

% load Neuron information (bad, celltype, sessioninfo, etc. for each day)
load(fullfile(dirs.data2load, "NeuroInfo.mat"))
% load neurons with significant increase or decrease in AZRZ (for PV
% selection)
load(fullfile(dirs.data2load, "allsess_shuffledSig_distance2RZ.mat"))
% load neuron firing rates per day (for PV selection)
load(fullfile(dirs.data2load, "allsess_raw_vs_residuals_distance2RZ.mat"))
% load ratemap for each day 
outmapdir = fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/");
%%
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.animals);
allindex = allindex(animal_idx,:); %filter based on animals to include 

[sessions, sessID] = unique(allindex(:,[1:2, 7]),'rows'); %session = [animalID, recording date, novelty day]

%%
correct_only = 1;
scale = 1;

trackb_preRZ_DecreaseIndex = cell([], 1); trackb_postRZ_DecreaseIndex = cell([], 1); 
trackc_preRZ_DecreaseIndex = cell([], 1); trackc_postRZ_DecreaseIndex = cell([], 1); 
% for calculating DecreaseIndex:
inc_minus_dec = 1;
normalize = 1;

RZrangestart = -30;
RZrangeend = -1*RZrangestart;
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
    % not bad, nov decreasing -15-20 deg, high firing (max > 30), NS interneurons
    NS_units = find(ismember(allsess_unitType_sh(currsess_idx), 'Narrow Interneuron'));
    unit2incl = intersect(notbadidx, NS_units);
   
    highfiring_units = find(mean(allsess_mean_fam_raw_all(currsess_idx,:), 2) > 20.1099);
    unit2incl = intersect(unit2incl, highfiring_units);
    
    if isempty(unit2incl)
        continue
    end
    % select trials:
    if correct_only == 1
        select_trials = find(outmap.bin2.nov.labels(:, 2) == 1 & outmap.bin2.nov.labels(:, 3) == 1); % RZ-centered and receive reward
    else
        select_trials = find(outmap.bin2.nov.labels(:, 2) == 1); % only select RZ_centered
    end
    if isempty(select_trials) || length(select_trials) < 5
        continue
    end
    % for each unit (including good or bad neurons of all types), regress out speed and lick rate
    % refer to Nuri's code:
    % Y:\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023\Code\Helpers\script_MultipleLinearRegression.m
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
    
    %% append data to all day table
    ratemap_nov_residual = ratemap_nov_residual(select_trials, :, unit2incl);
    % average across zones
    zones = outmap.bin2.nov.labels(select_trials, 4);
    zonetrial = zeros(size(zones));
    for iZone = 1:3
        iZone_trialidx = find(zones == iZone);
        zonetrial(iZone_trialidx) = 1:length(iZone_trialidx);
    end
    [averagedFiringRates, nEntries] = calculateAveragedFiringRates_zoneXspaceXneuron(ratemap_nov_residual, zonetrial);
    % average within zone trial blocks 
    novelexposureday = unique(allsess_sessinfo_sh(currsess_idx, 3));
    index = allindex(allindex(:, 1) == animalnum & allindex(:, 2) == day, :);
    env = unique(index(:, 6));
    if novelexposureday == 1 && ismember(2, env)
        trialblock = 25; 
    elseif novelexposureday == 2 && ismember(3, env)
        trialblock = 40;
    else
        continue
    end
    %% Use trial block t-test p-value to calculate DecreaseIndex: number of spatial bins
    % each cell of DecreaseIndex is a trial block, each element is a unit, the value
    % is the number of significantly decreased spatial bin
    ratemap_nov_residual_trialblock = arrayfun(@(iBlock) mean(ratemap_nov_residual(trialblock*(iBlock-1)+1 : trialblock*(iBlock), :, :), 1), 1:floor(size(ratemap_nov_residual, 1)/trialblock), 'UniformOutput', false);
    ratemap_nov_residual_trialblock = cell2mat(ratemap_nov_residual_trialblock');
    for iBlock = 1:1:size(ratemap_nov_residual_trialblock, 1)
        for iCell = 1:size(ratemap_nov_residual_trialblock, 3)
            RZrange = [RZrangestart RZrangeend]; RZstart = find(outmap.bin2.binEdges(2:end) == RZrange(1)); RZend = find(outmap.bin2.binEdges(2:end) == RZrange(2));
            RZloc_relative = find(outmap.bin2.binEdges(2:end) == 0) - RZstart;
            temp = ratemap_nov_residual(trialblock*(iBlock-1)+1 : trialblock*(iBlock), RZstart : RZend, iCell);
            vMax = max(temp,[],2);
            vMin = min(temp,[],2);
            temp = (temp - vMin) ./ (vMax - vMin);
            if scale == 1
                temp = temp - nanmean(temp(:, 1:2), 2);
            end
            preRZ_DecreaseIndex = 0; postRZ_DecreaseIndex = 0;
            for iBin = 1:RZloc_relative
                [h,p,~,stats] = ttest(temp(:,iBin));
                if p * length(1:RZloc_relative) < 0.05
                    preRZ_DecreaseIndex = preRZ_DecreaseIndex - 1;
                end
            end
            for iBin = RZloc_relative + 1 : size(temp, 2)
                [h,p,~,stats] = ttest(temp(:,iBin));
                if p * length(1:RZloc_relative) < 0.05
                    postRZ_DecreaseIndex = postRZ_DecreaseIndex - 1;
                end
            end
            if ismember(2, env)
                if length(trackb_preRZ_DecreaseIndex) < iBlock || isempty(trackb_preRZ_DecreaseIndex{iBlock})
                    trackb_preRZ_DecreaseIndex{iBlock} = [];
                end
                if length(trackb_postRZ_DecreaseIndex) < iBlock || isempty(trackb_postRZ_DecreaseIndex{iBlock})
                    trackb_postRZ_DecreaseIndex{iBlock} = [];
                end
                trackb_preRZ_DecreaseIndex{iBlock} = [trackb_preRZ_DecreaseIndex{iBlock}, preRZ_DecreaseIndex];
                trackb_postRZ_DecreaseIndex{iBlock} = [trackb_postRZ_DecreaseIndex{iBlock}, postRZ_DecreaseIndex];
            elseif ismember(3, env)
                if length(trackc_preRZ_DecreaseIndex) < iBlock || isempty(trackc_preRZ_DecreaseIndex{iBlock})
                    trackc_preRZ_DecreaseIndex{iBlock} = [];
                end
                if length(trackc_postRZ_DecreaseIndex) < iBlock || isempty(trackc_postRZ_DecreaseIndex{iBlock})
                    trackc_postRZ_DecreaseIndex{iBlock} = [];
                end
                trackc_preRZ_DecreaseIndex{iBlock} = [trackc_preRZ_DecreaseIndex{iBlock}, preRZ_DecreaseIndex];
                trackc_postRZ_DecreaseIndex{iBlock} = [trackc_postRZ_DecreaseIndex{iBlock}, postRZ_DecreaseIndex];
            end
        end
    end
end

%% ED Fig. 7C
figure('Unit', 'inches', 'Position', [1.1354 3.7188 3.3646 1.8958]); 
t = tiledlayout(1, 1);
nexttile; 
data = [];
errors = [];
scatter_x = cell(length(trackb_preRZ_DecreaseIndex), 2);
scatter_y = cell(length(trackb_preRZ_DecreaseIndex), 2);
for i = 1:length(trackb_preRZ_DecreaseIndex)
    data = [data; -[mean(trackb_preRZ_DecreaseIndex{i}), mean(trackb_postRZ_DecreaseIndex{i})]];
    stderr_preRZ = std(trackb_preRZ_DecreaseIndex{i}) / sqrt(length(trackb_preRZ_DecreaseIndex{i}));
    stderr_postRZ = std(trackb_postRZ_DecreaseIndex{i}) / sqrt(length(trackb_postRZ_DecreaseIndex{i}));
    errors = [errors; [stderr_preRZ, stderr_postRZ]];
    % for scatter plot
    jitter1 = rand(1, length(trackb_preRZ_DecreaseIndex{i}));
    jitter2 = rand(1, length(trackb_postRZ_DecreaseIndex{i}));
    scatter_x{i, 1} = repelem(i - 0.15, length(trackb_preRZ_DecreaseIndex{i})) + (jitter1 - mean(jitter1)) / 6;
    scatter_x{i, 2} = repelem(i + 0.15, length(trackb_postRZ_DecreaseIndex{i})) + (jitter2 - mean(jitter2)) / 6;
    scatter_y{i, 1} = -trackb_preRZ_DecreaseIndex{i};
    scatter_y{i, 2} = -trackb_postRZ_DecreaseIndex{i};
end
hBar = bar(data); hold on; 
x = nan(size(data));
for k = 1:size(data,2)
    x(:,k) = hBar(k).XEndPoints;
end
errorbar(x, data, zeros(size(errors)), errors, 'k', 'linestyle', 'none');
for i = 1:length(trackb_preRZ_DecreaseIndex)
    % scatter(scatter_x{i, 1}, scatter_y{i, 1}, 50,[0,    0.4470,    0.7410])
    % scatter(scatter_x{i, 2}, scatter_y{i, 2}, 50,[0.8500,    0.3250,    0.0980])
    s1=scatter(scatter_x{i, 1}, scatter_y{i, 1}, 10,   'MarkerEdgeColor','k', 'LineWidth',0.5);
    s2=scatter(scatter_x{i, 2}, scatter_y{i, 2}, 10,  'MarkerEdgeColor','k', 'LineWidth',0.5);
    % alpha(s1, 0.5)
    % alpha(s2, 0.5)
end
legend({'preRZ','postRZ'}, "Location", 'northeastoutside')
ylim([0 12])
title('Track B'); xlabel('trial block'); ylabel('# spatial bins with decrease')

makefigurepretty(gcf)
figname = 'ExtendedDataFigure07_C';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')