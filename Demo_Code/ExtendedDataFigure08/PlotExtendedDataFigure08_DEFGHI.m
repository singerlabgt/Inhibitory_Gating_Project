%% Jeong et al. 2023 MANUSCRIPT - FIGURE 04
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
celltypes = {'Narrow Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','Pyramidal Cell'}; %names used in our manuscript
iBad = cell_metrics.tags.Bad; 
iPV = setdiff(cell_metrics.groundTruthClassification.lightsensitive, iBad); 
iNarrow = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Narrow')), iBad);
iPyr = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Pyr')), iBad);

%load the most updated distance-based residual firing rate map across
%distance relative to reward zone (nUnits x nBins)
filename = getlatestfile_with_string(dirs.data2load, 'allsess_raw_vs_residuals_distance2RZ.mat');
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
concat_rzBinsforplot.fam = [];     
concat_rzBinsforplot.nov = [];   
concat_rzBinsextendAZRZ.fam = [];     
concat_rzBinsextendAZRZ.nov = [];    
spatialinfo_fam = [];       %Nx1, concatenated spatial information in familiar environment
spatialinfo_nov = [];       %Nx1, concatenated spatial information in novel environment
isPlaceMod_fam = [];        %Nx1 (binary), 1 if determined as spatially modulated in familiar environment, 0 if not
isPlaceMod_nov = [];        %Nx1 (binary), 1 if determined as spatially modulated in novel environment, 0 if not
famMap = [];                %Nx72, rate map in familiar environment across linearized track
novMap = [];                %Nx72, rate map in novel environment across linearized track

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
    
    
    numUnits = length(s_ratemap(ii).outmap.unitID);
    rz_fam = thetaRZ_fam(ii,:);
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
    if any(isnan(rz_nov)) %doing this mess bcuz of some fam-only animals
        bins_nov = nan(size(bins_fam));
        bins_nov_forplot = nan(size(bins_fam_forplot));
        bins_nov_extendAZRZ = nan(size(bins_fam_extendAZRZ));
    else
        temp = arrayfun(@ (x) find(binEdges(:,1)==rz_nov(x,1)),1:3);
        bins_nov = cell2mat( arrayfun(@ (x) temp(x):temp(x)+ (ceil(nDeg*2/5)-1), 1:3, 'UniformOutput', false));
        bins_nov_forplot = cell2mat( arrayfun(@ (x) temp(x)-4:temp(x)+ 8, 1:3, 'UniformOutput', false));
        bins_nov_extendAZRZ = cell2mat( arrayfun(@ (x) temp(x)-2:temp(x)+ (ceil(nDeg*2/5)-1)+2, 1:3, 'UniformOutput', false)); % for nongoal modulation - AS suggestion: 10 deg before AZ and 10 deg after RZ. 11/13/2023
    end
    concat_rzBins.nov = vertcat(concat_rzBins.nov, repmat(bins_nov, numUnits, 1));
    concat_rzBinsforplot.nov = vertcat(concat_rzBinsforplot.nov, repmat(bins_nov_forplot, numUnits, 1));
    concat_rzBinsextendAZRZ.nov = vertcat(concat_rzBinsextendAZRZ.nov, repmat(bins_nov_extendAZRZ, numUnits, 1));
end

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


%% PLOTTING
fig = figure('units','inch','position',[0, 0, 6.5, 3.5]);
t = tiledlayout(2,3,'TileSpacing','compact');

%% FIGURE 4D: Correlation of rate maps by goal-representing place cells (goal stim vs sham stim) 
ax = nexttile; box off; hold on;
RCorr_edges = -0.1:0.05:1;
disp('stats for fig. 4D (ratemap correlaion), goal cells')
for ii = 1:3
    disp(['day' num2str(ii)])
    temp = histcounts(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellgoalstim, intersect(iPeakInNovGoal, intersect(iPyr, novPCs))))), ...
        RCorr_edges,'Normalization','probability');
    disp(['goal stim: cell num = ' num2str(length(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellgoalstim, intersect(iPeakInNovGoal, intersect(iPyr, novPCs))))))), ...
        ', prctile = ' num2str(prctile(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellgoalstim, intersect(iPeakInNovGoal, intersect(iPyr, novPCs))))), [0, 25, 50, 75, 100]))])
    plot(RCorr_edges(2:end),cumsum(temp),'Color', params.colors_goalstim(ii,:),'LineWidth',2);
    xlabel('Ratemap correlation'); ylabel('Cumulative proportion of goal cells'); 
    temp = histcounts(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellshamstim, intersect(iPeakInNovGoal, intersect(iPyr, novPCs))))), ...
        RCorr_edges,'Normalization','probability');
    disp(['sham stim: cell num = ' num2str(length(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellshamstim, intersect(iPeakInNovGoal, intersect(iPyr, novPCs))))))), ...
        ', prctile = ' num2str(prctile(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellshamstim, intersect(iPeakInNovGoal, intersect(iPyr, novPCs))))), [0, 25, 50, 75, 100]))])
    plot(RCorr_edges(2:end),cumsum(temp),'Color', params.colors_shamstim(ii,:),'LineWidth',2);
end
ylim([0 1])
yticks([0 0.5 1])
xticks([0 0.5 1])
ylim([0,1])
set(ax, 'TickDir', 'out');
xlabel('Ratemap correlation'); ylabel('Cumulative proportion of goal cells'); 


%% FIGURE 4E: spatial info distribution by goal-representing place cells (goal stim vs sham stim)
ax = nexttile; box off; hold on;
mice2incl = params.goalshamMice;
goalstimsess = unique(allindex(ismember(allindex(:,1),mice2incl) & allindex(:,6) ~= 1 & allindex(:,9) == 1, [1:2,7]),'rows');
shamstimsess = unique(allindex(ismember(allindex(:,1),mice2incl) & allindex(:,6) ~= 1 & allindex(:,9) == 2, [1:2,7]),'rows');
disp('stats for fig. 4E (spatial info), goal cells')
% compare novel days only 
SI_edges = 0:0.1:max(spatialinfo_nov(ismember(concat_sessInfo(:,1), mice2incl) & isPlaceMod_nov==1));
for ii = 1:3 %nov only - goalstim
    temp0 = intersect(iPeakInNovGoal, intersect(iPyr, novPCs));
    temp = histcounts(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        goalstimsess(:,2)) & concat_sessInfo(:,3) == ii))),SI_edges,'Normalization','probability');
    disp(['day' num2str(ii)])
    disp(['goal stim: cell num = ' num2str(length(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        goalstimsess(:,2)) & concat_sessInfo(:,3) == ii))))), ...
        ', prctile = ' num2str(prctile(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        goalstimsess(:,2)) & concat_sessInfo(:,3) == ii))), [0, 25, 50, 75, 100]))])
    plot(SI_edges(2:end),cumsum(temp),'Color', params.colors_goalstim(ii,:),'LineWidth',2);
end
for ii = 1:3 %nov only - shamstim
    temp0 = intersect(iPeakInNovGoal, intersect(iPyr, novPCs));
    temp = histcounts(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        shamstimsess(:,2)) & concat_sessInfo(:,3) == ii))),SI_edges,'Normalization','probability');
    disp(['day' num2str(ii)])
    disp(['sham stim: cell num = ' num2str(length(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        shamstimsess(:,2)) & concat_sessInfo(:,3) == ii))))), ...
        ', prctile = ' num2str(prctile(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        shamstimsess(:,2)) & concat_sessInfo(:,3) == ii))), [0, 25, 50, 75, 100]))])
    plot(SI_edges(2:end),cumsum(temp),'Color', params.colors_shamstim(ii,:),'LineWidth',2);
end
ylim([0 1])
yticks([0 0.5 1])
ylabel('Cumulative prop. of goal cells')
xlabel('Spatial info. (bits/spike)')
set(ax, 'TickDir', 'out');
legend({'Goal 1','Goal 2','Goal 3','Sham 1','Sham 2','Sham 3'},'Location','southeast','Box','off')


%% FIGURE 4F: Percent goal-representing place cells (goal stim vs sham stim)
mice_pv = find(ismember(concat_sessInfo(:,1),params.goalshamMice));
sessinfo = concat_sessInfo(intersect(iPyr, mice_pv), :); 
%total goal modulated Pyr cells in novel in PV mice
novGoalModPyr = concat_sessInfo( intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv), :);
for ii=1:3
    
    temp = sessinfo(sessinfo(:,3)==ii & ismember(sessinfo(:,2),goalstimsess(:,2)),:);
    uniqsess = unique(temp,'rows');
    %#pyr cells per sess
    numPyr_goalstim{ii} = arrayfun(@ (x) sum(temp(:,1) == uniqsess(x,1) & temp(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));    
        
    %#goal-representing pyr cells per sess
    numNovGoalPyr_goalstim{ii} = arrayfun( @(x) sum(novGoalModPyr(:,1) == uniqsess(x,1)...
        & novGoalModPyr(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    perc_novgoalPyr_goalstim_sessInfo{ii} = [uniqsess, ones(size(uniqsess,1),1) .*2, ones(size(uniqsess,1),1)]; %nov = 2, goalstim = 1
    
    temp = sessinfo(sessinfo(:,3)==ii & ismember(sessinfo(:,2),shamstimsess(:,2)),:);
    uniqsess = unique(temp,'rows');
    %#pyr cells per sess
    numPyr_shamstim{ii} = arrayfun(@ (x) sum(temp(:,1) == uniqsess(x,1) & temp(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    %#goal-representing pyr cells per sess
    numNovGoalPyr_shamstim{ii} = arrayfun( @(x) sum(novGoalModPyr(:,1) == uniqsess(x,1)...
        & novGoalModPyr(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    
    %percentage of goal-representing pyr cells per sess
    perc_novgoalPyr_goalstim{ii} = numNovGoalPyr_goalstim{ii} ./ numPyr_goalstim{ii} * 100;
    perc_novgoalPyr_shamstim{ii} = numNovGoalPyr_shamstim{ii} ./ numPyr_shamstim{ii} * 100;
end
ax = nexttile; box off; hold on;
meanp = arrayfun( @(x) mean(perc_novgoalPyr_goalstim{x}),1:3);
err = arrayfun( @(x) std(perc_novgoalPyr_goalstim{x}) ./ sqrt(length(perc_novgoalPyr_goalstim{x})),1:3);
errorbar(1:3, meanp, err,'LineWidth',2,'Color',params.colors_goalstim(3,:))
meanp = arrayfun( @(x) mean(perc_novgoalPyr_shamstim{x}),1:3);
err = arrayfun( @(x) std(perc_novgoalPyr_shamstim{x}) ./ sqrt(length(perc_novgoalPyr_shamstim{x})),1:3);
errorbar(1:3, meanp, err,'LineWidth',2,'Color',params.colors_shamstim(3,:))
legend({'Goal Stim','Sham Stim'},'Box','off','Location','northoutside')
xlabel('Day')
ylabel('Goal cells (%)')
set(ax, 'TickDir', 'out');
xlim([0.5 3.5]); xticks(1:3);


%% FIGURE 4G: Correlation of rate maps by non-goal-representing place cells (goal stim vs sham stim) 
ax = nexttile; box off; hold on;
RCorr_edges = -0.1:0.05:1;
disp('stats for fig. 4G (ratemap correlation), non-goal cells')
for ii = 1:3
    disp(['day' num2str(ii)])
    temp = histcounts(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellgoalstim, intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs))))), ...
        RCorr_edges,'Normalization','probability');
    disp(['goal stim: cell num = ' num2str(length(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellgoalstim, intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs))))))), ...
        ', prctile = ' num2str(prctile(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellgoalstim, intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs))))), [0, 25, 50, 75, 100]))])
    plot(RCorr_edges(2:end),cumsum(temp),'Color', params.colors_goalstim(ii,:),'LineWidth',2);
    xlabel('Ratemap correlation'); ylabel('Cumulative proportion of goal cells'); 
    temp = histcounts(MapCorrelation.nov_trials(intersect(iCellshamstim, intersect(find(concat_sessInfo(:,3) == ii), intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs))))), ...
        RCorr_edges,'Normalization','probability');
    disp(['sham stim: cell num = ' num2str(length(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellshamstim, intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs))))))), ...
        ', prctile = ' num2str(prctile(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iCellshamstim, intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs))))), [0, 25, 50, 75, 100]))])
    plot(RCorr_edges(2:end),cumsum(temp),'Color', params.colors_shamstim(ii,:),'LineWidth',2);

end
ylim([0 1])
yticks([0 0.5 1])
ylim([0, 1])
xticks([0 0.5 1])
xlabel('Ratemap correlation'); ylabel('Cumulative proportion of non-goal cells'); 
set(ax, 'TickDir', 'out');


%% FIGURE 4H: Distribution of spatial information by non-goal-representing place cells (goal stim vs sham stim)
disp('stats for fig. 4H (spatial info), non-goal cells')
ax = nexttile; box off; hold on;
SI_edges = 0:0.1:max(spatialinfo_nov(ismember(concat_sessInfo(:,1), mice2incl) & isPlaceMod_nov==1));
for ii = 1:3 %nov only - goalstim
    temp0 = intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs));
    temp = histcounts(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        goalstimsess(:,2)) & concat_sessInfo(:,3) == ii))),SI_edges,'Normalization','probability');
    disp(['day' num2str(ii)])
    disp(['goal stim: cell num = ' num2str(length(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        goalstimsess(:,2)) & concat_sessInfo(:,3) == ii))))), ...
        ', prctile = ' num2str(prctile(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        goalstimsess(:,2)) & concat_sessInfo(:,3) == ii))), [0, 25, 50, 75, 100]))])
    plot(SI_edges(2:end),cumsum(temp),'Color', params.colors_goalstim(ii,:),'LineWidth',2);
end
for ii = 1:3 %nov only - shamstim
    temp0 = intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs));
    temp = histcounts(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        shamstimsess(:,2)) & concat_sessInfo(:,3) == ii))),SI_edges,'Normalization','probability');
    disp(['day' num2str(ii)])
    disp(['sham stim: cell num = ' num2str(length(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        shamstimsess(:,2)) & concat_sessInfo(:,3) == ii))))), ...
        ', prctile = ' num2str(prctile(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), mice2incl) & ismember(concat_sessInfo(:,2),...
        shamstimsess(:,2)) & concat_sessInfo(:,3) == ii))), [0, 25, 50, 75, 100]))])
    plot(SI_edges(2:end),cumsum(temp),'Color', params.colors_shamstim(ii,:),'LineWidth',2);
end
ylabel('Cumulative prop. of non-goal cells')
xlabel('Spatial info. (bits/spike)')
legend({'Goal 1','Goal 2','Goal 3','Sham 1','Sham 2','Sham 3'},'Location','southeast','Box','off')
set(ax, 'TickDir', 'out');
ylim([0 1])
yticks([0 0.5 1])

%% FIGURE 4I: Percent non-goal-representing place cells (goal stim vs sham stim)
novNonGoalModPyr = concat_sessInfo( intersect(intersect(iNoPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv), :);
for ii=1:3
    %%%%% goal stim %%%%%
    temp = sessinfo(sessinfo(:,3)==ii & ismember(sessinfo(:,2),goalstimsess(:,2)),:);
    uniqsess = unique(temp,'rows');
    %#pyr cells per sess
    numPyr_goalstim{ii} = arrayfun(@ (x) sum(temp(:,1) == uniqsess(x,1) & temp(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));        
    %#goal-representing pyr cells per sess
    numNovNonGoalPyr_goalstim{ii} = arrayfun( @(x) sum(novNonGoalModPyr(:,1) == uniqsess(x,1)...
        & novNonGoalModPyr(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    
    %%%%% sham stim %%%%%
    temp = sessinfo(sessinfo(:,3)==ii & ismember(sessinfo(:,2),shamstimsess(:,2)),:);
    uniqsess = unique(temp,'rows');
    %#pyr cells per sess
    numPyr_shamstim{ii} = arrayfun(@ (x) sum(temp(:,1) == uniqsess(x,1) & temp(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    %#goal-representing pyr cells per sess
    numNovGoalPyr_shamstim{ii} = arrayfun( @(x) sum(novNonGoalModPyr(:,1) == uniqsess(x,1)...
        & novNonGoalModPyr(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
       
    %percentage of goal-representing pyr cells per sess
    perc_novNongoalPyr_goalstim{ii} = numNovNonGoalPyr_goalstim{ii} ./ numPyr_goalstim{ii} * 100;
    perc_novNongoalPyr_shamstim{ii} = numNovGoalPyr_shamstim{ii} ./ numPyr_shamstim{ii} * 100;    
end

ax = nexttile; box off; hold on;
meanp = arrayfun( @(x) mean(perc_novNongoalPyr_goalstim{x}),1:3);
err = arrayfun( @(x) std(perc_novNongoalPyr_goalstim{x}) ./ sqrt(length(perc_novNongoalPyr_goalstim{x})),1:3);
errorbar(1:3, meanp, err,'LineWidth',2,'Color',params.colors_goalstim(3,:))
meanp = arrayfun( @(x) mean(perc_novNongoalPyr_shamstim{x}),1:3);
err = arrayfun( @(x) std(perc_novNongoalPyr_shamstim{x}) ./ sqrt(length(perc_novNongoalPyr_shamstim{x})),1:3);
errorbar(1:3, meanp, err,'LineWidth',2,'Color',params.colors_shamstim(3,:))
xlabel('Day')
ylabel('Non-goal cells (%)')
xlim([0.5 3.5]); xticks(1:3);
set(ax, 'TickDir', 'out');


%% save figure
makefigurepretty(gcf)
figname = 'ExtendedDataFigure08_DEFGHI';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')