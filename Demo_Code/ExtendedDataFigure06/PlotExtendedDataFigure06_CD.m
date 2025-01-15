%% include both correct and incorrect trials, opto-tagged PV and high-firing NS units in WT NOV
clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 
baseBins = 1:2; %first 2 bins to average across; to be used as baseline firing rate

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end
%load data 
% load(fullfile('Y:\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\script\submission2\stimRatemapChangeBaseBin', 'AllMap'));
load(fullfile(dirs.data2load, "AllMap.mat"))
load(fullfile(dirs.data2load, "allsess_shuffledSig_distance2RZ.mat"));
load(fullfile(dirs.data2load, "allsess_raw_vs_residuals_distance2RZ.mat"), "allsess_mean_fam_raw_all");
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
iBad = cell_metrics.tags.Bad; %bad units to exclude
celltypes = {'PV Interneuron', 'Narrow Interneuron', 'Pyramidal Cell', 'Pyramidal Cell'};
for ct = 1:length(celltypes)
    iCT{ct} = setdiff(find(strcmp(cell_metrics.putativeCellType,celltypes{ct})), iBad)';
    if ct==1
        lightsensitive = cell_metrics.groundTruthClassification.lightsensitive;
        narrow = setdiff(find(startsWith(cell_metrics.putativeCellType,'Narrow')), iBad)';
        iCT{ct} = intersect(lightsensitive, narrow);
        
    end
    % only include high firing units for NS interneurons
    if ct == 2
        for iCell = 1:length(iCT{ct})
            cellid = iCT{ct}(iCell); unitID = cell_metrics.cluID(cellid);
            sessstr = cell_metrics.sessionName{cellid};
            sessstr_split = split(sessstr, '_');
            anim = sessstr_split{1}(2:end);
            recday = sessstr_split{2};
            tmpid = find(allsess_sessinfo_sh(:, 1) == str2num(anim) & allsess_sessinfo_sh(:, 2) == str2num(recday) & allsess_unitID_sh == unitID);
            if length(tmpid) ~= 1
                iCT{ct}(iCell) = -1;
                continue
            end
            if mean(allsess_mean_fam_raw_all(tmpid,:), 2) < 20.1099
                iCT{ct}(iCell) = -1;
            end
        end
        iCT{ct} = iCT{ct}(iCT{ct} > 0);        

    end
   
end
%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);
sessions = unique(allindex(:,[1:2,7]), 'rows'); %define session as one date
%determine animalID, recDay, novDay associated with each unit for indexing
numUnits = length(cell_metrics.sessionName);
cellMetSessInfo = nan(numUnits,3);
cellMetSessInfo_extracted = arrayfun( @(x) extract(cell_metrics.sessionName{x}, digitsPattern), 1:numUnits, 'UniformOutput', false);
cellMetSessInfo(:,1) = arrayfun( @(x) str2num(cellMetSessInfo_extracted{x}{1}), 1:numUnits); %#ok<ST2NM>
cellMetSessInfo(:,2) = arrayfun( @(x) str2num(cellMetSessInfo_extracted{x}{2}), 1:numUnits); %#ok<ST2NM>

%variables
stimInt = [0.5, 1.6, 0, 0];
stimNames = {'Low','High','None', 'None'};

%create FR data structure for plotting
clear data; fdNames = {'qFR_raw','qFR_norm'};
for f = 1:length(fdNames)
    data.(fdNames{f}) = [];
    for sInt = 1:length(stimNames)
        for ee = 1:length(params.environments) %for fam and nov environments
            env = params.environments{ee};
            for ct = 1:length(iCT)
                data2add = AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}])(iCT{ct},:);
                unitID = AllMap.(env).unitID.(['Volt_' stimNames{sInt}])(iCT{ct},:);
                sessInfo = cellMetSessInfo(iCT{ct},:);
                cellT = ones(length(unitID),1) .* ct; 
                Environ = ones(length(unitID),1) .* ee;
                Intensity = ones(length(unitID),1) .* sInt;
                data2add = [data2add, sessInfo, unitID, cellT, Environ, Intensity]; %
                data.(fdNames{f}) = [data.(fdNames{f}); data2add];
            end
        end
    end
    temp = data.(fdNames{f});
    data.(fdNames{f}) = temp(~isnan(temp(:,1)),:);
end

fdNames = {'deltaFRraw','deltaFRnorm'};
for f = 1:length(fdNames)
    data.(fdNames{f}) = [];
    for sInt = 1:length(stimNames)
        for ee = 1:length(params.environments) %for fam and nov environments
            env = params.environments{ee};
            for ct = 1:length(iCT)
                sessInfo = [];
                data2add = AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}])(iCT{ct},:);
                unitID = AllMap.(env).unitID.(['Volt_' stimNames{sInt}])(iCT{ct},:);
                sessInfo = cellMetSessInfo(iCT{ct},:);
                cellT = ones(length(unitID),1) .* ct; 
                Environ = ones(length(unitID),1) .* ee;
                Intensity = ones(length(unitID),1) .* sInt;
                data2add = data2add(all(~isnan(data2add),2),:); %excl NaN
                data.(fdNames{f}) = [data.(fdNames{f}); data2add];
            end
        end
    end
end

%% PLOT STIM EFFECTS ON FIRING RATES PER INTENSITY: ALL SESSIONS COMBINED FROM ALL ANIMALS
%identify session info
binEdges_FR = AllMap.binEdges_firingrate;
xTime = mean(getBinEdges(binEdges_FR),2); 
xTime_startbin = find(xTime >= -1.0250, 1); xTime_endbin = find(xTime <= 3.0250, 1, 'last');
xTime = xTime(xTime_startbin:xTime_endbin);
stim_famgoal = unique(allindex(allindex(:,6) == 1 & allindex(:,9) == 1, [1:2,7]),'rows'); %includes all mice
stim_novgoal = unique(allindex(allindex(:,6) ~= 1 & allindex(:,9) == 1, [1:2,7]),'rows');
stim_novsham = unique(allindex(allindex(:,6) ~= 1 & allindex(:,9) == 2, [1:2,7]),'rows');
nostim_nov = unique(allindex(allindex(:,6) ~= 1 & allindex(:,9) == 0, [1:2,7]),'rows');
fdNames = {'norm'}; %fdNames = {'raw','norm'};

for f = 1:length(fdNames)
    datTrace = data.(['deltaFR' fdNames{f}]);
    datQ = data.(['qFR_' fdNames{f}])(:,1);
    sessInfo = data.(['qFR_' fdNames{f}])(:,2:end);
    iNovGoalStim = ismember(sessInfo(:,1), params.goalshamMice) & ismember(sessInfo(:,2), stim_novgoal(:,2)) & sessInfo(:,6) == 2;
    iNovNoStim = ismember(sessInfo(:,1), params.WTmice) & sessInfo(:,6) == 2;
    disp(['num of iNovGoalStim = ' num2str(length(unique(sessInfo(iNovGoalStim, [1:2, 4]), 'rows')))])
    disp(['num of iNovNoStim = ' num2str(length(unique(sessInfo(iNovNoStim, [1:2,4]), 'rows')))])
    figure('Unit', 'inches', 'Position', [1.0208 1.7292 7 6]); clear h; 
    t = tiledlayout(2,length(stimInt),'TileSpacing','compact');
    for ct = [1 3]
        for sInt = 1:length(stimInt) %exclude the highest stim intensity (3.2V)
            if sInt == 4
                sigbar_pos_goalstim = 0.04; sigbar_pos_wt = 0.05;
            else
                sigbar_pos_goalstim = 0.6; sigbar_pos_wt = 0.7;
            end
    
            h(sInt) = nexttile; ax = h(sInt); 
            ops.ax = h(sInt);
            ops.error = 'sem';
            ops.x_axis = xTime;
            ops.alpha = 0.2;
            ops.line_width = 2;
            %novgoalstim
            ops.color_area = params.colors_goalstim(3,:);
            ops.color_line = params.colors_goalstim(3,:);
            plot_areaerrorbar(datTrace(iNovGoalStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt, xTime_startbin:xTime_endbin), ops); hold on
            %t-test with Bonferroni correction for multiple comparisons
            temp = datTrace(iNovGoalStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt, xTime_startbin:xTime_endbin);
            colors = params.colors_goalstim(3,:);
            
            binSize = 0.1;
            numBins = floor((xTime(end) - xTime(1)) / binSize); % Number of 0.05s bins
            ht = nan(numBins, 1);
            p = nan(numBins, 1);
            stats = cell(numBins, 1);
    
            % Perform t-tests for each bin
            for iBin = 1:numBins
                binStart = (iBin - 1) * binSize + xTime(1); % Bin start time in seconds
                binEnd = binStart + binSize; % Bin end time in seconds
                binIndices = find(xTime >= binStart & xTime < binEnd); % Indices of xTime within the bin
    
                if ~isempty(binIndices)
                    binData = temp(:, binIndices);
    
                    % Aggregate data within the bin (e.g., take the mean across the bin)
                    binMean = mean(binData, 2);
    
                    % Perform t-test on the aggregated data
                    [ht(iBin), p(iBin), ~, stats{iBin}] = ttest(binMean);
                end
            end
            
            sig = find(p < 0.05 / numel(p));
            arrayfun( @(ii) scatter(ax, (sig(ii) - 1) * binSize + xTime(1), sigbar_pos_goalstim,...
                'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
            box off; hold on; % added 07/14/2024
    
            %nostim, novel, added 07/14/2024
            ops.color_area = params.colors_nov(3,:);
            ops.color_line = params.colors_nov(3,:);
            plot_areaerrorbar(datTrace(iNovNoStim & sessInfo(:,5) == ct + 1, xTime_startbin:xTime_endbin), ops); hold on
            %t-test with Bonferroni correction for multiple comparisons
            temp = datTrace(iNovNoStim & sessInfo(:,5) == 2, xTime_startbin:xTime_endbin);
            colors = params.colors_nov(3,:);
            
            binSize = 0.1;
            numBins = floor((xTime(end) - xTime(1)) / binSize); % Number of 0.05s bins
            ht = nan(numBins, 1);
            p = nan(numBins, 1);
            stats = cell(numBins, 1);
    
            % Perform t-tests for each bin
            for iBin = 1:numBins
                binStart = (iBin - 1) * binSize + xTime(1); % Bin start time in seconds
                binEnd = binStart + binSize; % Bin end time in seconds
                binIndices = find(xTime >= binStart & xTime < binEnd); % Indices of xTime within the bin
    
                if ~isempty(binIndices)
                    binData = temp(:, binIndices);
    
                    % Aggregate data within the bin (e.g., take the mean across the bin)
                    binMean = mean(binData, 2);
    
                    % Perform t-test on the aggregated data
                    [ht(iBin), p(iBin), ~, stats{iBin}] = ttest(binMean);
                end
            end
            
            sig = find(p < 0.05 / numel(p));
            arrayfun( @(ii) scatter(ax, (sig(ii) - 1) * binSize + xTime(1), sigbar_pos_wt,...
                'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
            box off; hold on;
            axis square
            xline(h(sInt),0,'k:');yline(h(sInt),0,'k:');set(h(sInt), 'TickDir', 'out')
            title(h(sInt),stimNames(sInt));
            xlabel('Time to stim onset (s)')
            ylabel(['Norm FR change of ' celltypes{ct}])
            xlim([-1 3])
            if sInt == 4; ylim(ax, [-0.15 0.06]); end 
        end
    end
end
%% save figure
makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure06_CD';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')