%% Jeong et al. 2023 MANUSCRIPT - FIGURE 03_C,D,E
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

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','Pyramidal Cell'}; %names used in our manuscript
iBad = cell_metrics.tags.Bad; 
iPV = setdiff(cell_metrics.groundTruthClassification.lightsensitive, iBad); 
iNarrow = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Narrow')), iBad);
iPyr = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Pyr')), iBad);

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
    rz_fam = [wrapTo360(rz_fam - nDeg); wrapTo360(rz_fam + nDeg)]';
    temp = arrayfun(@ (x) find(binEdges(:,1)==rz_fam(x,1)),1:3);
    bins_fam = cell2mat( arrayfun(@ (x) temp(x):temp(x)+ (ceil(nDeg*2/5)-1), 1:3, 'UniformOutput', false));
    concat_rzBins.fam = vertcat(concat_rzBins.fam, repmat(bins_fam, numUnits, 1));
    
    rz_nov = thetaRZ_nov(ii,:);
    rz_nov = [wrapTo360(rz_nov - nDeg); wrapTo360(rz_nov + nDeg)]';
    if any(isnan(rz_nov)) %doing this mess bcuz of some fam-only animals
        bins_nov = nan(size(bins_fam));
    else
        temp = arrayfun(@ (x) find(binEdges(:,1)==rz_nov(x,1)),1:3);
        bins_nov = cell2mat( arrayfun(@ (x) temp(x):temp(x)+ (ceil(nDeg*2/5)-1), 1:3, 'UniformOutput', false));
    end
    concat_rzBins.nov = vertcat(concat_rzBins.nov, repmat(bins_nov, numUnits, 1));
end


%% identify goal-representing place fields
%count multiple fields per cell, if available; find units with peak rate
%within goal areas, defined as 10 deg before & after reward zone start angle (i.e. AZ+RZ)
[~, I] = nanmax(famMap,[], 2);
iPeakInFamGoal = find(ismember(I, concat_rzBins.fam)); %indices of PFs with a peak firing in AZ/RZ in at least once
iNoPeakInFamGoal = find(~ismember(I, concat_rzBins.fam)); %indices of PFs with a peak firing outside AZ/RZ
[~, I] = nanmax(novMap,[], 2);
iPeakInNovGoal = find(ismember(I, concat_rzBins.nov));
iNoPeakInNovGoal = find(~ismember(I, concat_rzBins.nov));

%identify putative place cells (spatial info > shuffled info's 95th
%percentile)
famPCs = find(isPlaceMod_fam == 1);
novPCs = find(isPlaceMod_nov == 1);

%% rate map correlation over trials for wild-type mice
dataCorr.famMap_wtPeakInFamGoal = MapCorrelation.fam_trials(intersect(iPeakInFamGoal,intersect(iPyr, find(ismember(concat_sessInfo(:,1), params.WTmice))))); 
dataCorr.famMap_wtPeakInNovGoal = MapCorrelation.fam_trials(intersect(iPeakInNovGoal,intersect(iPyr, find(ismember(concat_sessInfo(:,1), params.WTmice))))); 
dataCorr.famMap_wtNoPeakInFamGoal = MapCorrelation.fam_trials(intersect(iNoPeakInFamGoal,intersect(iPyr, find(ismember(concat_sessInfo(:,1), params.WTmice))))); 
dataCorr.famMap_wtNoPeakInNovGoal = MapCorrelation.fam_trials(intersect(iNoPeakInNovGoal,intersect(iPyr, find(ismember(concat_sessInfo(:,1), params.WTmice))))); 


fig = figure('units','inch','position',[0, 0, 6.5, 2]);
t = tiledlayout(1,3,'TileSpacing','compact');

%% FIGURE 3C: rate map correlation over trials
nexttile; box off; hold on
RCorr_edges = -0.1:0.05:1;
for ii = 1:3
    temp = histcounts(MapCorrelation.fam_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iPeakInFamGoal, intersect(intersect(iPyr, famPCs), find(ismember(concat_sessInfo(:,1), params.WTmice)))))), ...
        RCorr_edges,'Normalization','probability');
    plot(RCorr_edges(2:end),cumsum(temp),'Color', params.colors_fam(ii,:),'LineWidth',2);
    disp(['ratemap corr, fam, day' num2str(ii) ', cellN=' num2str(size(MapCorrelation.fam_trials(intersect(find(concat_sessInfo(:,3) == ii), intersect(iPeakInFamGoal, intersect(intersect(iPyr, famPCs), find(ismember(concat_sessInfo(:,1), params.WTmice)))))), 2))])
    xlabel('Ratemap correlation'); ylabel('Cumulative proportion of goal cells'); 
end
for ii=1:3
    temp = histcounts(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii),intersect(iPeakInNovGoal, intersect(intersect(iPyr, novPCs), find(ismember(concat_sessInfo(:,1), params.WTmice)))))), ...
        RCorr_edges,'Normalization','probability');
    disp(['ratemap corr, nov, day' num2str(ii) ', cellN=' num2str(size(MapCorrelation.nov_trials(intersect(find(concat_sessInfo(:,3) == ii),intersect(iPeakInNovGoal, intersect(intersect(iPyr, novPCs), find(ismember(concat_sessInfo(:,1), params.WTmice)))))), 2))])

    plot(RCorr_edges(2:end),cumsum(temp),'Color', params.colors_nov(ii,:),'LineWidth',2);
end
ylim([0,1])
set(gca, 'TickDir', 'out')


%% FIGURE 3D: pyrmaidal cell spatial info distribution
nexttile; box off; hold on;
SI_edges = 0:0.1:max(spatialinfo_fam(ismember(concat_sessInfo(:,1), params.WTmice) & isPlaceMod_fam==1));
for ii = 1:3 %fam
    temp0 = intersect(iPeakInFamGoal, intersect(iPyr, famPCs));
    temp = histcounts(spatialinfo_fam( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), params.WTmice) ...
        & concat_sessInfo(:,3) == ii))),SI_edges,'Normalization','probability');
    disp(['SI, fam, day' num2str(ii) ', cellN=' num2str(size(spatialinfo_fam( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), params.WTmice) ...
        & concat_sessInfo(:,3) == ii))), 1))])

    plot(SI_edges(2:end),cumsum(temp),'Color', params.colors_fam(ii,:),'LineWidth',2);
end
for ii = 1:3 %nov
    temp0 = intersect(iPeakInNovGoal, intersect(iPyr, novPCs));
    temp = histcounts(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), params.WTmice) ...
        & concat_sessInfo(:,3) == ii))),SI_edges,'Normalization','probability');
    disp(['SI, nov, day' num2str(ii) ', cellN=' num2str(size(spatialinfo_nov( intersect(temp0,...
        find(ismember(concat_sessInfo(:,1), params.WTmice) ...
        & concat_sessInfo(:,3) == ii))), 1))])
    plot(SI_edges(2:end),cumsum(temp),'Color', params.colors_nov(ii,:),'LineWidth',2);
end
ylabel('Cumulative proportion of goal cells')
xlabel('Spatial information (bits/spike)')
ylim([0 1]);
legend({'Fam 1','Fam 2','Fam 3',...
    'Nov 1','Nov 2','Nov 3'},'Location','southeast','Box','off')


%% FIGURE 3E: goal-representing Pyr cells across days - averaged over sessions
%WT fam+nov
mice_wt = find(ismember(concat_sessInfo(:,1),params.WTmice));
sessinfo = concat_sessInfo(intersect(iPyr, mice_wt), :); 
for ii=1:3
    temp = sessinfo(sessinfo(:,3)==ii,:);
    uniqsess = unique(temp,'rows');
    %#pyr cells per sess
    numPyr_fam{ii} = arrayfun(@ (x) sum(temp(:,1) == uniqsess(x,1) & temp(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    numPyr_nov{ii} = arrayfun(@ (x) sum(temp(:,1) == uniqsess(x,1) & temp(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    
    famGoalModPyr = concat_sessInfo( intersect(intersect(iPeakInFamGoal, intersect(iPyr, famPCs)), mice_wt), :);
    novGoalModPyr = concat_sessInfo( intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_wt), :);
    %#goal-representing pyr cells per sess
    numFamGoalPyr{ii} = arrayfun( @(x) sum(famGoalModPyr(:,1) == uniqsess(x,1) & famGoalModPyr(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    numNovGoalPyr{ii} = arrayfun( @(x) sum(novGoalModPyr(:,1) == uniqsess(x,1) & novGoalModPyr(:,2) == uniqsess(x,2)), 1:size(uniqsess,1));
    
    %percentage of goal-representing pyr cells per sess
    perc_famgoalPyr{ii} = numFamGoalPyr{ii} ./ numPyr_fam{ii} * 100;
    perc_novgoalPyr{ii} = numNovGoalPyr{ii} ./ numPyr_nov{ii} * 100;
    perc_famgoalPyr_sessInfo{ii} = [uniqsess, ones(size(uniqsess,1),1)]; %fam = 1
    perc_novgoalPyr_sessInfo{ii} = [uniqsess, ones(size(uniqsess,1),1) .* 2]; %nov = 2
end

nexttile; hold on;
meanp = arrayfun( @(x) mean(perc_famgoalPyr{x}),1:3);
err = arrayfun( @(x) std(perc_famgoalPyr{x}) ./ sqrt(length(perc_famgoalPyr{x})),1:3);
errorbar(1:3, meanp, err,'LineWidth',2,'Color',params.colors_fam(end,:))
meanp = arrayfun( @(x) mean(perc_novgoalPyr{x}),1:3);
err = arrayfun( @(x) std(perc_novgoalPyr{x}) ./ sqrt(length(perc_novgoalPyr{x})),1:3);
errorbar(1:3, meanp, err,'LineWidth',2,'Color',params.colors_nov(2,:))
legend({'Fam','Nov'},'Box','off','Location','northeast')
xlabel('Day')
ylabel('Mean +/- SEM goal cells (%)')
xlim([0.5 3.5]);

%% save figure
makefigurepretty(gcf)
figname = 'Figure03_CDE';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')