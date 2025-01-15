%% load data
% Jeong et al. 2023 MANUSCRIPT - FIGURE 04
% NJeong 03/23/2023

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%
%load default parameters
[dirs, params] = getDefaultParameters(maindir); 
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);
sessions = unique(allindex(:,1:2),'rows');

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Wide Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','WS Interneuron','Pyramidal Cell'}; %names used in our manuscript
iBad = cell_metrics.tags.Bad; 
iPV = setdiff(cell_metrics.groundTruthClassification.lightsensitive, iBad); 
iNarrow = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Narrow')), iBad);
iWide = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Wide')), iBad);
iPyr = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Pyr')), iBad);

%load the most updated distance-based residual firing rate map across
%distance relative to reward zone (nUnits x nBins)
filename = getlatestfile_with_string(dirs.data2load, 'distance2RZ.mat');
load(fullfile(dirs.data2load, filename));

%identify cell types for all units 
cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
    1:length(celltypes),'UniformOutput',false);

%load place coding data structure 
load(fullfile(dirs.data2load, 'placecodingout_09212024.mat'));


%% create structure 
binEdges = getBinEdges(0:5:360);
nDeg = 10; %in degrees; distance before and after to count as 'goal'

%initialize 
concat_sessInfo = [];       %Nx3, concatenated list of [animal ID, date, novel day]
concat_rzBins.fam = [];     %Nx12, concatenated bin numbers for reward zone in familiar environment
concat_rzBins.nov = [];     %Nx12,concatenated bin numbers for reward zone in novel environment
spatialinfo_fam = [];       %Nx1, concatenated spatial information in familiar environment
spatialinfo_nov = [];       %Nx1, concatenated spatial information in novel environment
isPlaceMod_fam = [];        %Nx1 (binary), 1 if determined as spatially modulated in familiar environment, 0 if not
isPlaceMod_nov = [];        %Nx1 (binary), 1 if determined as spatially modulated in novel environment, 0 if not
famMap = [];                %Nx72, rate map in familiar environment across linearized track
novMap = [];                %Nx72, rate map in novel environment across linearized track
famMap_pertrial = {};
novMap_pertrial = {};
concat_rzBinsforplot.fam = [];     %Nx12, concatenated bin numbers for reward zone in familiar environment
concat_rzBinsforplot.nov = [];     %Nx12,concatenated bin numbers for reward zone in novel environment
concat_rzBinsextendAZRZ.fam = [];     %Nx12, concatenated bin numbers for reward zone in familiar environment
concat_rzBinsextendAZRZ.nov = [];     %Nx12,concatenated bin numbers for reward zone in novel environment
for ii = 1:length(sessions)
    currSess = sessions(ii,:);
    iSess = find(sessions(:,1) == currSess(1) & sessions(:,2) == currSess(2));
    spatialinfo_fam = vertcat(spatialinfo_fam, s_ratemap(iSess).outmap.bin5.fam.spatialinfo_original');
    
    unitID = s_ratemap(iSess).outmap.unitID;
    for iUnit = 1:length(unitID)
        tempInfo = s_ratemap(iSess).outmap.bin5.fam.spatialinfo_original(iUnit);
        threshold_95th = prctile(s_ratemap(iSess).outmap.bin5.fam.spatialinfo_shuffled(:,iUnit), 95); %95th percentile shuffled ratemap
        isPlaceMod_fam = [isPlaceMod_fam; tempInfo > threshold_95th];
        
        if isempty(s_ratemap(iSess).outmap.bin5.nov.ratemap)
            isPlaceMod_nov = [isPlaceMod_nov; nan];
        else
            tempInfo = s_ratemap(iSess).outmap.bin5.nov.spatialinfo_original(iUnit);
            threshold_95th = prctile(s_ratemap(iSess).outmap.bin5.nov.spatialinfo_shuffled(:,iUnit), 95); %95th percentile shuffled ratemap
            isPlaceMod_nov = [isPlaceMod_nov; tempInfo > threshold_95th];
        end
    end
    
    concat_sessInfo = vertcat(concat_sessInfo, repmat(currSess, length(unitID), 1));
    
    famMap = vertcat(famMap, s_ratemap(iSess).outmap.bin5.fam.ratemap_original);
    if isempty(s_ratemap(iSess).outmap.bin5.nov.ratemap)
        s_ratemap(iSess).outmap.bin5.nov.ratemap_original =...
            nan(size(s_ratemap(iSess).outmap.bin5.fam.ratemap_original));
        spatialinfo_nov = vertcat(spatialinfo_nov, nan(size( s_ratemap(iSess).outmap.bin5.fam.spatialinfo_original')));
    else
        spatialinfo_nov = vertcat(spatialinfo_nov, s_ratemap(iSess).outmap.bin5.nov.spatialinfo_original');
    end
    novMap = vertcat(novMap, s_ratemap(iSess).outmap.bin5.nov.ratemap_original); 

    % for plotting example cell (11/12/2023)
    clearvars m n p i
    tmp1 = s_ratemap(iSess).outmap.bin5.fam.ratemap;
    [m,n,p] = size(tmp1); tmpc1 = cell(1, p); for i = 1:p tmpc1{i} = tmp1(:, :, i); end
    famMap_pertrial = [famMap_pertrial, tmpc1];
    clearvars m n p i
    tmp2 = s_ratemap(iSess).outmap.bin5.nov.ratemap;
    [m,n,p] = size(tmp2); tmpc2 = cell(1, p); for i = 1:p tmpc2{i} = tmp2(:, :, i); end
    if ~isempty(tmp2)
        novMap_pertrial = [novMap_pertrial, tmpc2];
    else
        novMap_pertrial = [novMap_pertrial, cell([1, size(tmp1, 3)])];
    end
    
    numUnits = length(s_ratemap(ii).outmap.unitID);
    rz_fam = thetaRZ_fam(ii,:);
    % rz_fam = [wrapTo360(rz_fam - nDeg); wrapTo360(rz_fam + nDeg)]';
    rz_fam = [wrapTo360(rz_fam); wrapTo360(rz_fam + nDeg)]';
    temp = arrayfun(@ (x) find(binEdges(:,1)==rz_fam(x,1)),1:3);
    bins_fam = cell2mat( arrayfun(@ (x) temp(x):temp(x)+ (ceil(nDeg*2/5)-1), 1:3, 'UniformOutput', false));
    concat_rzBins.fam = vertcat(concat_rzBins.fam, repmat(bins_fam, numUnits, 1));
    bins_fam_forplot = cell2mat( arrayfun(@ (x) temp(x)-4:temp(x)+ 8, 1:3, 'UniformOutput', false)); % doesn't change the classification, just for visualization
    concat_rzBinsforplot.fam = vertcat(concat_rzBinsforplot.fam, repmat(bins_fam_forplot, numUnits, 1));
    bins_fam_extendAZRZ = cell2mat( arrayfun(@ (x) temp(x)-2:temp(x)+ (ceil(nDeg*2/5)-1)+2, 1:3, 'UniformOutput', false)); % for nongoal modulation - AS suggestion: 10 deg before AZ and 10 deg after RZ. 11/13/2023
    concat_rzBinsextendAZRZ.fam = vertcat(concat_rzBinsextendAZRZ.fam, repmat(bins_fam_extendAZRZ, numUnits, 1));
    
    rz_nov = thetaRZ_nov(ii,:);
    rz_nov = [wrapTo360(rz_nov - nDeg); wrapTo360(rz_nov + nDeg)]';
    % rz_nov = [wrapTo360(rz_nov); wrapTo360(rz_nov + nDeg)]'; % debug, generate RZ-only fig4EH
    if any(isnan(rz_nov)) %doing this mess bcuz of some fam-only animals
        bins_nov = nan(size(bins_fam));
        bins_nov_forplot = nan(size(bins_fam_forplot));
        bins_nov_extendAZRZ = nan(size(bins_fam_extendAZRZ));
    else
        temp = arrayfun(@ (x) find(binEdges(:,1)==rz_nov(x,1)),1:3);
        bins_nov = cell2mat( arrayfun(@ (x) temp(x):temp(x)+ (ceil(nDeg*2/5)-1), 1:3, 'UniformOutput', false));
        bins_nov_forplot = cell2mat( arrayfun(@ (x) temp(x)-4:temp(x)+ 8, 1:3, 'UniformOutput', false));
        bins_nov_extendAZRZ = cell2mat( arrayfun(@ (x) temp(x)-2:temp(x)+ (ceil(nDeg*2/5)-1)+2, 1:3, 'UniformOutput', false)); % for nongoal modulation - AS suggestion: 10 deg before AZ and 10 deg after RZ. 11/13/2023
        % bins_nov_extendAZRZ = cell2mat( arrayfun(@ (x) temp(x)-4:temp(x)+(ceil(nDeg*2/5)-1), 1:3, 'UniformOutput', false)); % debug, generate RZ-only fig4EH
    end
    concat_rzBins.nov = vertcat(concat_rzBins.nov, repmat(bins_nov, numUnits, 1));
    concat_rzBinsforplot.nov = vertcat(concat_rzBinsforplot.nov, repmat(bins_nov_forplot, numUnits, 1));
    concat_rzBinsextendAZRZ.nov = vertcat(concat_rzBinsextendAZRZ.nov, repmat(bins_nov_extendAZRZ, numUnits, 1));

end
%%
%identify goal-representing place fields
%count multiple fiels per cell, if available; find units with peak rate
%within goal areas, defined as 10 deg before & after reward zone start angle (i.e. AZ+RZ)
[~, I] = nanmax(famMap,[], 2);
iPeakInFamGoal = find(ismember(I, concat_rzBins.fam)); %indices of PFs with a peak firing in AZ/RZ in at least once
iNoPeakInFamGoal = find(~ismember(I, concat_rzBinsextendAZRZ.fam)); %indices of PFs with a peak firing outside AZ/RZ
[~, I] = nanmax(novMap,[], 2);
iPeakInNovGoal = find(ismember(I, concat_rzBins.nov));
iNoPeakInNovGoal = find(~ismember(I, concat_rzBinsextendAZRZ.nov));

%identify putative place cells (spatial info > shuffled info's 95th
%percentile)
famPCs = find(isPlaceMod_fam == 1);
novPCs = find(isPlaceMod_nov == 1);
%rate map correlation over trials
iSess_goalstim = unique(...
    allindex(ismember(allindex(:,1),params.goalshamMice) & allindex(:,6) ~=1 & allindex(:,9) ==1,[1:2,7]),'rows'); 
iSess_shamstim = unique(...
    allindex(ismember(allindex(:,1),params.goalshamMice) & allindex(:,6) ~=1 & allindex(:,9) ==2,[1:2,7]),'rows'); 
iCellgoalstim = find(ismember(concat_sessInfo(:,1:2), iSess_goalstim(:,1:2),'rows'));
iCellshamstim = find(ismember(concat_sessInfo(:,1:2), iSess_shamstim(:,1:2),'rows'));

dataCorr.famMap_pvPeakInFamGoalwGoalStim = MapCorrelation.fam_trials(intersect(iCellgoalstim, intersect(iPeakInFamGoal, iPyr)));
dataCorr.famMap_pvPeakInFamGoalwShamStim = MapCorrelation.fam_trials(intersect(iCellshamstim, intersect(iPeakInFamGoal, iPyr)));
dataCorr.famMap_pvPeakInNovGoalwGoalStim = MapCorrelation.fam_trials(intersect(iCellgoalstim, intersect(iPeakInNovGoal, iPyr)));
dataCorr.famMap_pvPeakInNovGoalwShamStim = MapCorrelation.fam_trials(intersect(iCellshamstim, intersect(iPeakInNovGoal, iPyr)));
dataCorr.famMap_pvNoPeakInFamGoalwGoalStim = MapCorrelation.fam_trials(intersect(iCellgoalstim, intersect(iNoPeakInFamGoal, iPyr)));
dataCorr.famMap_pvNoPeakInFamGoalwShamStim = MapCorrelation.fam_trials(intersect(iCellshamstim, intersect(iNoPeakInFamGoal, iPyr)));
dataCorr.famMap_pvNoPeakInNovGoalwGoalStim = MapCorrelation.fam_trials(intersect(iCellgoalstim, intersect(iNoPeakInNovGoal, iPyr)));
dataCorr.famMap_pvNoPeakInNovGoalwShamStim = MapCorrelation.fam_trials(intersect(iCellshamstim, intersect(iNoPeakInNovGoal, iPyr)));

dataCorr.novMap_pvPeakInFamGoalwGoalStim = MapCorrelation.nov_trials(intersect(iCellgoalstim, intersect(iPeakInFamGoal, iPyr)));
dataCorr.novMap_pvPeakInFamGoalwShamStim = MapCorrelation.nov_trials(intersect(iCellshamstim, intersect(iPeakInFamGoal, iPyr)));
dataCorr.novMap_pvPeakInNovGoalwGoalStim = MapCorrelation.nov_trials(intersect(iCellgoalstim, intersect(iPeakInNovGoal, iPyr)));
dataCorr.novMap_pvPeakInNovGoalwShamStim = MapCorrelation.nov_trials(intersect(iCellshamstim, intersect(iPeakInNovGoal, iPyr)));
dataCorr.novMap_pvNoPeakInFamGoalwGoalStim = MapCorrelation.nov_trials(intersect(iCellgoalstim, intersect(iNoPeakInFamGoal, iPyr)));
dataCorr.novMap_pvNoPeakInFamGoalwShamStim = MapCorrelation.nov_trials(intersect(iCellshamstim, intersect(iNoPeakInFamGoal, iPyr)));
dataCorr.novMap_pvNoPeakInNovGoalwGoalStim = MapCorrelation.nov_trials(intersect(iCellgoalstim, intersect(iNoPeakInNovGoal, iPyr)));
dataCorr.novMap_pvNoPeakInNovGoalwShamStim = MapCorrelation.nov_trials(intersect(iCellshamstim, intersect(iNoPeakInNovGoal, iPyr)));

%% Fig 4AB: plotting example goal/non-goal modulated cells.
% AS suggestion: collapse all RZs into one. Maybe change the scaling so
% it's easier to see the contrast
mice2incl = params.goalshamMice;
goalstimsess = unique(allindex(ismember(allindex(:,1),mice2incl) & allindex(:,6) ~= 1 & allindex(:,9) == 1, [1:2,7]),'rows');
shamstimsess = unique(allindex(ismember(allindex(:,1),mice2incl) & allindex(:,6) ~= 1 & allindex(:,9) == 2, [1:2,7]),'rows');
mice_pv = find(ismember(concat_sessInfo(:,1),params.goalshamMice));

% all track, but line plot
fig0 = figure('units','inch','position',[0, 0, 6.5, 4]);
t0 = tiledlayout('flow');
colors = params.colors_pyr;
%% goal stim
%% goal modulated cells

sessinfo = concat_sessInfo(intersect(find(ismember(concat_sessInfo(:,2), goalstimsess(:,2))), intersect(iPyr, mice_pv)), :); 
%total goal modulated Pyr cells in novel in PV mice
novGoalModPyr = concat_sessInfo( intersect(find(ismember(concat_sessInfo(:,2), goalstimsess(:,2))),intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv)), :);
novGoalModPyr_cells = intersect(find(ismember(concat_sessInfo(:,2), goalstimsess(:,2))),intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv));
length(novGoalModPyr_cells) % 117


for iCell = [7, 22]
% for iCell = 1:length(novGoalModPyr_cells)
    ax = nexttile;
    celli = novGoalModPyr_cells(iCell);
    temp = novMap_pertrial{celli};
    % temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;
    ops.ax     = ax;
    ops.x_axis = binEdges(:,1);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;
    rzs = concat_rzBins.nov(celli, [3, 7, 11]);
    for rz = rzs
        if rz > 72
            rz = rz - 72;
        end
        if rz <= 0
            rz = 72 - rz;
        end
    end
    title(['goal-stim, goal'])
    % set(gca, 'XTick',[1,180, 360],...
    %             'XTickLabel',[0, 180, 360]);

    if iCell == 7
        xlim([230-40 230+40]);
            xticks([230-40 230 230+40])

    else
        xlim([150-40 150+40])
            xticks([150-40 150 150+40])

    end
end
%% non-goal modulated cells
sessinfo = concat_sessInfo(intersect(find(ismember(concat_sessInfo(:,2), goalstimsess(:,2))), intersect(iPyr, mice_pv)), :); 
%total goal modulated Pyr cells in novel in PV mice
novNonGoalModPyr = concat_sessInfo( intersect(find(ismember(concat_sessInfo(:,2), goalstimsess(:,2))),intersect(intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv)), :);
novNonGoalModPyr_cells = intersect(find(ismember(concat_sessInfo(:,2), goalstimsess(:,2))),intersect(intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv));
length(novNonGoalModPyr_cells) % 137 -> 94 after extendAZRZ
% all track, but line plot
for iCell = [68, 75]%1:length(novNonGoalModPyr_cells)
    ax = nexttile;
    celli = novNonGoalModPyr_cells(iCell);
    temp = novMap_pertrial{celli};
    % temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;
    ops.ax     = ax;
    ops.x_axis = binEdges(:,1);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;
    rzs = concat_rzBins.nov(celli, [3, 7, 11]);
    for rz = rzs
        if rz > 72
            rz = rz - 72;
        end
        if rz <= 0
            rz = 72 - rz;
        end
    end
    if iCell == 68
        xlim([120-40 120+40])
    else
        xlim([120-40 120+40])
    end
        xticks([120-40 120 120+40])

    title(['goal, Non goal'])
    % set(gca, 'XTick',[1,180, 360],...
    %             'XTickLabel',[0, 180, 360]);
end
%% sham stim
%% goal modulated cells
sessinfo = concat_sessInfo(intersect(find(ismember(concat_sessInfo(:,2), shamstimsess(:,2))), intersect(iPyr, mice_pv)), :); 
%total goal modulated Pyr cells in novel in PV mice
novGoalModPyr = concat_sessInfo( intersect(find(ismember(concat_sessInfo(:,2), shamstimsess(:,2))),intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv)), :);
novGoalModPyr_cells = intersect(find(ismember(concat_sessInfo(:,2), shamstimsess(:,2))),intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv));
length(novGoalModPyr_cells) % 164

for iCell = [27, 38]%1:length(novGoalModPyr_cells)
    ax = nexttile;
    celli = novGoalModPyr_cells(iCell);
    temp = novMap_pertrial{celli};
    % temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;
    ops.ax     = ax;
    ops.x_axis = binEdges(:,1);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;
    rzs = concat_rzBins.nov(celli, [3, 7, 11]);
    for rz = rzs
        if rz > 72
            rz = rz - 72;
        end
        if rz <= 0
            rz = 72 - rz;
        end
    end
    if iCell == 27
        xlim([120-40 120+40])
    else
        xlim([120-40 120+40])
    end
        xticks([120-40 120 120+40])

    title(['sham, goal'])
    % set(gca, 'XTick',[1,180, 360],...
    %             'XTickLabel',[0, 180, 360]);
end
%% non-goal modulated cells
sessinfo = concat_sessInfo(intersect(find(ismember(concat_sessInfo(:,2), shamstimsess(:,2))), intersect(iPyr, mice_pv)), :); 
%total goal modulated Pyr cells in novel in PV mice
novNonGoalModPyr = concat_sessInfo( intersect(find(ismember(concat_sessInfo(:,2), shamstimsess(:,2))),intersect(intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv)), :);
novNonGoalModPyr_cells = intersect(find(ismember(concat_sessInfo(:,2), shamstimsess(:,2))),intersect(intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv));
length(novNonGoalModPyr_cells) %? -> 73 after extendAZRZ

% all track, but line plot

for iCell = [56, 68]%1:length(novNonGoalModPyr_cells)
    ax = nexttile;
    celli = novNonGoalModPyr_cells(iCell);
    temp = novMap_pertrial{celli};
    % temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;
    ops.ax     = ax;
    ops.x_axis = binEdges(:,1);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;
    
    rzs = concat_rzBins.nov(celli, [3, 7, 11]);
    for rz = rzs
        if rz > 72
            rz = rz - 72;
        end
        if rz <= 0
            rz = 72 - rz;
        end
    end
    if iCell == 56
        xlim([150-40 150+40])

    else
        xlim([150-40 150+40])
    end
    ylim([0 5]) 
    xticks([150-40 150 150+40])
    title(['sham, Non goal'])
end
%% save figure
makefigurepretty(gcf)
figname = 'ExtendedDataFigure08_BC';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
