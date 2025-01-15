clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%%
allRZs = 0; % for goal stim, if including all RZs or including only the stimulated zone (log/high, excluding 0mV)
%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), setdiff(params.animals, 4)),:); % exclude X4. because this animal only has sham stim. 
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


mice2incl = params.goalshamMice;
goalstimsess = unique(allindex(ismember(allindex(:,1),mice2incl) & allindex(:,6) ~= 1 & allindex(:,9) == 1, [1:2,7]),'rows');
shamstimsess = unique(allindex(ismember(allindex(:,1),mice2incl) & allindex(:,6) ~= 1 & allindex(:,9) == 2, [1:2,7]),'rows');
mice_pv = find(ismember(concat_sessInfo(:,1),params.goalshamMice));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% peak firing rate at RZ in goal-modulated cells: difference between goal-stim and sham-stim sessions
% NOV
novGoalModPyr_goalstim = intersect(find(ismember(concat_sessInfo(:,2), goalstimsess(:,2))),intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv));
novGoalModPyr_shamstim = intersect(find(ismember(concat_sessInfo(:,2), shamstimsess(:,2))),intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv));

length(novGoalModPyr_goalstim) % 117 place cells
length(novGoalModPyr_shamstim) % 163

% based on thetaRZ_nov and sessioninfo, get RZ-centered FR (-60~60 deg)
f = arrayfun(@(x) find(ismember(sessions, concat_sessInfo(x, :), 'rows')), 1:size(concat_sessInfo, 1));
novGoalModPyr_thetaRZ = thetaRZ_nov(f, :);
[nanr, nanc] = find(isnan(novGoalModPyr_thetaRZ));
novGoalModPyr_thetaRZ(nanr, nanc) = 10000;
novGoalModPyr_RZrange = [];
for iRZ = 1:3
    RZrange = vertcat(arrayfun(@(x) novGoalModPyr_thetaRZ(x, iRZ) - 60 : 5 : novGoalModPyr_thetaRZ(x, iRZ) + 65, 1:size(novGoalModPyr_thetaRZ, 1), 'UniformOutput', false))';
    novGoalModPyr_RZrange = [novGoalModPyr_RZrange, RZrange];
end
novGoalModPyr_RZrange = cell2mat(novGoalModPyr_RZrange);
novGoalModPyr_goalstim_RZrange = novGoalModPyr_RZrange(novGoalModPyr_goalstim, :);
novGoalModPyr_shamstim_RZrange = novGoalModPyr_RZrange(novGoalModPyr_shamstim, :);

novGoalModPyr_goalstim_RZrange(novGoalModPyr_goalstim_RZrange > 5000) = []; % goal cells are selected not on NAN thetaRZ days, so this step shouldn't affect the result
novGoalModPyr_shamstim_RZrange(novGoalModPyr_shamstim_RZrange > 5000) = [];

novGoalModPyr_goalstim_RZrange = wrapTo360(novGoalModPyr_goalstim_RZrange); novGoalModPyr_goalstim_RZrange(novGoalModPyr_goalstim_RZrange == 360) = 0;
novGoalModPyr_shamstim_RZrange = wrapTo360(novGoalModPyr_shamstim_RZrange); novGoalModPyr_shamstim_RZrange(novGoalModPyr_shamstim_RZrange == 360) = 0;

% convert to binEdge index
novGoalModPyr_goalstim_RZrange = novGoalModPyr_goalstim_RZrange ./ 5 + 1;
novGoalModPyr_shamstim_RZrange = novGoalModPyr_shamstim_RZrange ./ 5 + 1;
%% for goal stim, only include the reward zones with goal stim
if allRZs==0
    sess_novGoalModPyr_goalstim = concat_sessInfo(intersect(find(ismember(concat_sessInfo(:,2), goalstimsess(:,2))),intersect(intersect(iPeakInNovGoal, intersect(iPyr, novPCs)), mice_pv)), :);
    mask_goalstim = []; % only include the trials with goal stim
    for iSess = 1:size(sess_novGoalModPyr_goalstim, 1)
        if sess_novGoalModPyr_goalstim(iSess, 1) == 4
            iden = 'X';
        else
            iden = 'N';
        end
        fname=dir(fullfile(maindir, 'Demo_Data', 'singlesess_ratemap_distance2RZ', [iden num2str(sess_novGoalModPyr_goalstim(iSess, 1)) '_' num2str(sess_novGoalModPyr_goalstim(iSess, 2)) '_*.mat']));
        load(fullfile(maindir, 'Demo_Data', 'singlesess_ratemap_distance2RZ', fname.name));
        hasStim = arrayfun(@(iRZ) max(outmap.bin5.nov.stimVoltage{iRZ}) > 0, 1:3);
        trialspatialsize = size(outmap.bin5.nov.ratemap, 2);
        tmp_mask = [repelem(hasStim(1), trialspatialsize), repelem(hasStim(2), trialspatialsize), repelem(hasStim(3), trialspatialsize)];
        mask_goalstim = [mask_goalstim; tmp_mask];
    end
    tmp = nan([size(novGoalModPyr_goalstim_RZrange, 1), size(novGoalModPyr_goalstim_RZrange, 2)*2/3]);
    for iRow = 1:size(tmp, 1)
        tmp1 = novGoalModPyr_goalstim_RZrange(iRow, :) .* mask_goalstim(iRow, :);
        tmp(iRow, :) = tmp1(tmp1~=0);
    end
    novGoalModPyr_goalstim_RZrange = tmp;
end
novGoalModPyr_goalstim_fr = arrayfun(@(x) novMap(novGoalModPyr_goalstim(x), novGoalModPyr_goalstim_RZrange(x, :)), 1:size(novGoalModPyr_goalstim_RZrange, 1), 'UniformOutput', false);% firstly try avg trials, if bad then per trial.
novGoalModPyr_goalstim_fr = cell2mat(novGoalModPyr_goalstim_fr');
novGoalModPyr_shamstim_fr = arrayfun(@(x) novMap(novGoalModPyr_shamstim(x), novGoalModPyr_shamstim_RZrange(x, :)), 1:size(novGoalModPyr_shamstim_RZrange, 1), 'UniformOutput', false);% firstly try avg trials, if bad then per trial.
novGoalModPyr_shamstim_fr = cell2mat(novGoalModPyr_shamstim_fr');
% concat eachreward trial vertically
blocksize = size(novGoalModPyr_shamstim_fr, 2) / 3;
if allRZs == 0
    novGoalModPyr_goalstim_fr = [novGoalModPyr_goalstim_fr(:, 1:blocksize); novGoalModPyr_goalstim_fr(:, blocksize + 1 : 2 * blocksize)];
    novGoalModPyr_shamstim_fr = [novGoalModPyr_shamstim_fr(:, 1:blocksize); novGoalModPyr_shamstim_fr(:, blocksize + 1 : 2 * blocksize); novGoalModPyr_shamstim_fr(:, 2 * blocksize + 1 : 3 * blocksize)];
else
    novGoalModPyr_goalstim_fr = [novGoalModPyr_goalstim_fr(:, 1:blocksize); novGoalModPyr_goalstim_fr(:, blocksize + 1 : 2 * blocksize); novGoalModPyr_goalstim_fr(:, 2 * blocksize + 1 : 3 * blocksize)];
    novGoalModPyr_shamstim_fr = [novGoalModPyr_shamstim_fr(:, 1:blocksize); novGoalModPyr_shamstim_fr(:, blocksize + 1 : 2 * blocksize); novGoalModPyr_shamstim_fr(:, 2 * blocksize + 1 : 3 * blocksize)];
end
normalizeFR = 0; scaleFR = 0;
% normalize
if normalizeFR == 1
    novGoalModPyr_goalstim_fr = (novGoalModPyr_goalstim_fr - min(novGoalModPyr_goalstim_fr, [], 2)) ./ (max(novGoalModPyr_goalstim_fr, [], 2) - min(novGoalModPyr_goalstim_fr, [], 2));
    novGoalModPyr_shamstim_fr = (novGoalModPyr_shamstim_fr - min(novGoalModPyr_shamstim_fr, [], 2)) ./ (max(novGoalModPyr_shamstim_fr, [], 2) - min(novGoalModPyr_shamstim_fr, [], 2));
end
% scale
if scaleFR == 1
    novGoalModPyr_goalstim_fr = (novGoalModPyr_goalstim_fr - nanmean(novGoalModPyr_goalstim_fr(:, 1:2), 2)) .* 100;
    novGoalModPyr_shamstim_fr = (novGoalModPyr_shamstim_fr - nanmean(novGoalModPyr_shamstim_fr(:, 1:2), 2)) .* 100;
end
%%  plot
figure('Units', 'inches', 'Position', [5.4167 3.1771 5.55 1.65]);
 t = tiledlayout(1, 2);
ax = nexttile; hold on
for iCond = 1:2
    if iCond == 1
        temp = novGoalModPyr_goalstim_fr;
        colors = params.colors_goalstim;
    elseif iCond == 2
        temp = novGoalModPyr_shamstim_fr;
        colors = params.colors_shamstim;
    end

    ops.ax     = ax;
    ops.x_axis = position_binEdges(1:end-1);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;
    % ylim([-0.2, 0.6]);
end
%t-test with Bonferroni correction for multiple comparisons
h = nan(size(temp,2),1);
p = nan(size(temp,2),1);
stats = cell(size(temp,2),1);
for iT = 1:size(temp,2)
    [h(iT),p(iT),~,stats{iT}] = ttest2(novGoalModPyr_goalstim_fr(:,iT), novGoalModPyr_shamstim_fr(:, iT));
end
sig = find(p < 0.05 / numel(p));
arrayfun( @(ii) scatter(ax, position_binEdges(sig(ii)), 2.5,...
    'Marker','|','MarkerEdgeColor',params.colors_goalstim(end,:)),1:length(sig));
ylabel('Firing rate in Hz')
xlim([-65 65])
xticks([-60 0 60])


nexttile;
% compare peak firing rate using violins. peak firing rate is defined as
% max firing rte within goals
maxL = max(size(novGoalModPyr_goalstim_fr, 1), size(novGoalModPyr_shamstim_fr, 1));
temp = nan(maxL,2);
temp(1:size(novGoalModPyr_goalstim_fr, 1), 1) = max(novGoalModPyr_goalstim_fr(:,11:16), [], 2);
temp(1:size(novGoalModPyr_shamstim_fr, 1), 2) = max(novGoalModPyr_shamstim_fr(:,11:16), [], 2);
v1 = violinplot_half(temp);
colormat = [params.colors_goalstim(2,:); params.colors_shamstim(2,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 25;
end
ylabel('Peak firing rate (Hz)')
xlabel('Goal stim (blue) vs Sham stim (orange)')
xticks([])

disp(['there are ', num2str(size(novGoalModPyr_goalstim_fr,1)), ' goal place cells for peak fr'])
disp(['there are ', num2str(size(novGoalModPyr_shamstim_fr,1)), ' sham place cells for peak fr'])

makefigurepretty(gcf,1)
savefigALP([figdir '/'], 'SuppFigure07_AB', 'filetype', 'pdf')
