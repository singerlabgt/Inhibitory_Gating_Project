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
allindex = allindex(animal_idx,:); %filter based on animals to include 

[sessions, sessID] = unique(allindex(:,[1:2, 7]),'rows'); %session = [animalID, recording date, novelty day]

%%
correct_only = 1;
scale = 1;

allday_ratemap = zeros(0, 65);
allday_trialidx = zeros(0, 1);
allday_novenv = zeros(0, 1);
allday_novelexposureday = zeros(0, 1);
allday_iZone = zeros(0, 1);
allday_zonetrial = zeros(0, 1); % trial num on each zone
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
    unit2incl = NS_units;
    highfiring_units = find(mean(allsess_mean_fam_raw_all(currsess_idx,:), 2) > 20.1099);
    unit2incl = intersect(unit2incl, highfiring_units); 
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
    ratemap_nov_residual = ratemap_nov_residual(select_trials, :, notbadidx(unit2incl));
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

%% Fig. 3G
fig = figure('units','inch','position',[0 0 6 2]);
t = tiledlayout(1,3,'TileSpacing','compact');

RZrange = [-30,30]; 
RZstart = find(outmap.bin2.binEdges(2:end) == RZrange(1)); RZend = find(outmap.bin2.binEdges(2:end) == RZrange(2));
RZloc_relative = find(outmap.bin2.binEdges(2:end) == 0) - RZstart;
x_axis_original = outmap.bin2.binEdges(2:end);
tilelocs = [1 2 3; 4 5 6];
for env = 2
    nexttile(tilelocs(env-1,1)); hold on; box off; 
    if env == 2
        select_exposureday = 1; trialblock = 25;
    elseif env == 3
        select_exposureday = 2; trialblock = 40;
    end
    day1_select_trials = allday_novelexposureday == select_exposureday & allday_novenv == env;
    novelday1_ratemap = allday_ratemap(day1_select_trials, :); novelday1_zonetrial = allday_trialidx(day1_select_trials);
    [novelday1_ratemap, novelday1_nEntries] = calculateAveragedFiringRates(novelday1_ratemap, novelday1_zonetrial);
    %normlize
    novelday1_ratemap = novelday1_ratemap(:, RZstart : RZend);
    vMax = max(novelday1_ratemap,[],2);
    vMin = min(novelday1_ratemap,[],2);
    novelday1_ratemap = (novelday1_ratemap - vMin) ./ (vMax - vMin);
    novelday1_ratemap_noscale = novelday1_ratemap;
    if scale == 1
        novelday1_ratemap = novelday1_ratemap - nanmean(novelday1_ratemap(:, 1:2), 2);
    end
    x_axis_plot = x_axis_original(RZstart : RZend);

    imagesc(x_axis_plot,1:trialblock, novelday1_ratemap_noscale(1:trialblock, :))
    xlim([-30 30])
    xticks([-30 0 30])
    ylim([1 trialblock])
    xlabel('Distance to novel RZ (deg)')
    ylabel('Trial');set(gca, 'YDir', 'reverse')
    title('first block')

    % plot the first trial of day1
    % fig = figure();
    nexttile(tilelocs(env-1,3)); 
    ops.ax = gca;
    x_axis_original = outmap.bin2.binEdges(2:end);
    % first block
    iBlock = 1;
    ops.x_axis = mean([x_axis_original(RZstart : RZend); x_axis_original(RZstart : RZend)]); %mean(getBinEdges(time_binEdges),2);
    ops.color_area = params.colors_narrowInt(1, :);
    ops.color_line = params.colors_narrowInt(1, :);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(novelday1_ratemap(1:trialblock, :), ops); hold on;

    %t-test with Bonferroni correction for multiple comparisons
    x_axis_plot = x_axis_original(RZstart : RZend);
    temp = novelday1_ratemap(1:trialblock, :);
    h = nan(size(temp,2),1);
    p = nan(size(temp,2),1);
    stats = cell(size(temp,2),1);
    for iT = 1:size(temp,2)
        % [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT), 0, 'Tail', 'left');
        [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT), 0, 'Tail', 'left');
    end
    sig = find(p < 0.05 / numel(p));
    arrayfun( @(ii) scatter(gca, x_axis_plot(sig(ii)), 0.4,...
        'Marker','|','MarkerEdgeColor',params.colors_narrowInt(1, :)),1:length(sig));

    ops.color_area = params.colors_narrowInt(3, :);
    ops.color_line = params.colors_narrowInt(3, :);
    % last block
    if env == 2
        iBlock = floor(size(novelday1_ratemap, 1)/trialblock);
    else
        iBlock = floor(size(novelday1_ratemap, 1)/trialblock) - 1; % track C has too few trials at last block, so take the last second block
    end
    % last_trials = trialblock*(floor(size(novelday1_ratemap, 1)/(trialblock)-1)-1)+1 : trialblock*(floor(size(novelday1_ratemap, 1)/trialblock)-1);
    last_trials = trialblock*(iBlock-1)+1 : trialblock*(iBlock);
    plot_areaerrorbar(novelday1_ratemap(last_trials, :), ops); hold on

    %t-test with Bonferroni correction for multiple comparisons
    x_axis_plot = x_axis_original(RZstart : RZend);
    temp = novelday1_ratemap(last_trials, :);
    h = nan(size(temp,2),1);
    p = nan(size(temp,2),1);
    stats = cell(size(temp,2),1);
    for iT = 1:size(temp,2)
        % [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT), 0, 'Tail', 'left');
        [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT), 0, 'Tail', 'left');
    end
    sig = find(p < 0.05 / numel(p));
    arrayfun( @(ii) scatter(gca, x_axis_plot(sig(ii)), 0.5,...
        'Marker','|','MarkerEdgeColor',params.colors_narrowInt(3, :)),1:length(sig));

    xline(0, '--'); yline(0, '--'); 
    xlim([-30 30])
    xticks([-30 0 30])
    ylim([-0.65 0.65])
    yticks([-0.6 -0.4 -0.2 0 0.2 0.4 0.6])
    xlabel('Distance to novel RZ (deg)')
    ylabel('Normalized residual firing rate'); 
    title(['env' num2str(env) ', trialblock = ' num2str(trialblock)])    

    nexttile(tilelocs(env-1,2)); hold on;
    imagesc(x_axis_plot,1:trialblock, novelday1_ratemap_noscale(last_trials, :))
    xlim([-30 30])
    xticks([-30 0 30])
    ylim([1 length(last_trials)])
    xlabel('Distance to novel RZ (deg)')
    ylabel('Trial'); set(gca, 'YDir', 'reverse')
    title('last block')
    colormap(gray)
end
%% save figure
makefigurepretty(gcf)
figname = 'Figure04_FG';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
