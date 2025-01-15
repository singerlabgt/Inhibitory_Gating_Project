%% Jeong et al. 2023 MANUSCRIPT - FIGURE02_D,E
% NJeong 03/28/2023

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

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
iBad = cell_metrics.tags.Bad; %bad units to exclude
celltypes = {'PV Interneuron', 'Pyramidal Cell'};
for ct = 1:length(celltypes)
    tmp = setdiff(find(strcmp(cell_metrics.putativeCellType,celltypes{ct})), iBad)';
    iCT{ct} = tmp(tmp<=3465);
    if ct==1
        lightsensitive = cell_metrics.groundTruthClassification.lightsensitive;
        narrow = setdiff(find(startsWith(cell_metrics.putativeCellType,'Narrow')), iBad)';
        tmp = intersect(lightsensitive, narrow);
        iCT{ct} = tmp(tmp<=3465);
    end
end

%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);
sessions = unique(allindex(:,[1:2,7]), 'rows'); %define session as one date

%load data 
load(fullfile(dirs.data2load, 'AllMap'));

%determine animalID, recDay, novDay associated with each unit for indexing
numUnits = length(cell_metrics.sessionName);
cellMetSessInfo = nan(numUnits,3);
cellMetSessInfo_extracted = arrayfun( @(x) extract(cell_metrics.sessionName{x}, digitsPattern), 1:numUnits, 'UniformOutput', false);
cellMetSessInfo(:,1) = arrayfun( @(x) str2num(cellMetSessInfo_extracted{x}{1}), 1:numUnits); %#ok<ST2NM>
cellMetSessInfo(:,2) = arrayfun( @(x) str2num(cellMetSessInfo_extracted{x}{2}), 1:numUnits); %#ok<ST2NM>

%variables
stimInt = [0, 0.5, 1.6, 3.2];
stimNames = {'None','Low','High','Highest'};


%% Fig. 2D-E: firing rate comparison between goal stim and sham stim 

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
                cellT = ones(length(unitID),1) .* ct; %CellT of 1-4 indicate PV, narrow, wide, pyr, in that order
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
                cellT = ones(length(unitID),1) .* ct; %CellT of 1-4 indicate PV, narrow, wide, pyr, in that order
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

stim_famgoal = unique(allindex(allindex(:,6) == 1 & allindex(:,9) == 1, [1:2,7]),'rows'); %includes all mice
stim_novgoal = unique(allindex(allindex(:,6) ~= 1 & allindex(:,9) == 1, [1:2,7]),'rows');
stim_novsham = unique(allindex(allindex(:,6) ~= 1 & allindex(:,9) == 2, [1:2,7]),'rows');

fdNames = {'norm'}; %fdNames = {'raw','norm'};
pv_sham = [];
pv_goal = [];
pyr_sham = [];
pyr_goal = [];
for f = 1:length(fdNames)
    datTrace = data.(['deltaFR' fdNames{f}]);
    datQ = data.(['qFR_' fdNames{f}])(:,1);
    sessInfo = data.(['qFR_' fdNames{f}])(:,2:end);
    iNovGoalStim = ismember(sessInfo(:,1), params.goalshamMice) & ismember(sessInfo(:,2), stim_novgoal(:,2)) & sessInfo(:,6) == 2;
    iNovShamStim = ismember(sessInfo(:,1), params.goalshamMice) & ismember(sessInfo(:,2), stim_novsham(:,2)) & sessInfo(:,6) == 2;
    % 1 is PV, 2 is pyr
    for ct = 1:length(celltypes)
        clear h
        figure;
        t = tiledlayout(2,length(stimInt)-1,'TileSpacing','compact');

        for sInt = 1:length(stimInt)-1 %exclude the highest stim intensity (3.2V)
            if ct == 1 % pv
                pv_goal = [pv_goal; sum(iNovGoalStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt)];
                pv_sham = [pv_sham; sum(iNovShamStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt)];
            elseif ct == 2 % pyr
                pyr_goal = [pyr_goal; sum(iNovGoalStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt)];
                pyr_sham = [pyr_sham; sum(iNovShamStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt)];
            end
            h(sInt) = nexttile;
            ops.ax = h(sInt);
            ops.error = 'sem';
            ops.x_axis = xTime;
            ops.alpha = 0.2;
            ops.line_width = 2;
            ops.color_area = params.colors_shamstim(3,:);
            ops.color_line = params.colors_shamstim(3,:);
            plot_areaerrorbar(datTrace(iNovShamStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt, :), ops);
            box off; hold on;
            
            %novgoalstim
            ops.color_area = params.colors_goalstim(3,:);
            ops.color_line = params.colors_goalstim(3,:);
            plot_areaerrorbar(datTrace(iNovGoalStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt, :), ops);
            axis square
            xline(h(sInt),0,'k:');set(h(sInt), 'TickDir', 'out')
            title(h(sInt),stimNames(sInt));
        end
        
        for sInt = 1:length(stimInt)-1 %exclude the highest stim intensity (3.2V)
            h(sInt+3) = nexttile; box off; hold on;
            novSS = datQ(iNovShamStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt, :);
            scatter(h(sInt+3), ones(length(novSS),1), novSS, 30, 'Jitter', 'on',...
                'JitterAmount', 0.1, 'MarkerEdgeColor', params.colors_shamstim(3,:),...
                'MarkerFaceColor', params.colors_shamstim(3,:),'MarkerFaceAlpha', 0.2 )
            
            novGS = datQ(iNovGoalStim & sessInfo(:,5) == ct & sessInfo(:,7) == sInt, :);
            scatter(h(sInt+3), ones(length(novGS),1).* 2, novGS , 30, 'Jitter', 'on',...
                'JitterAmount', 0.1, 'MarkerEdgeColor', params.colors_goalstim(3,:),...
                'MarkerFaceColor', params.colors_goalstim(3,:),'MarkerFaceAlpha', 0.2 )
            errorbar(h(sInt+3), 1:2, [mean(novSS), mean(novGS)], [std(novSS)./sqrt(length(novSS)), std(novGS)./sqrt(length(novGS))],'k-');
            xticks(h(sInt+3), 1:2); xticklabels(h(sInt+3),{'Sham','Goal'})
            set(h(sInt+3), 'TickDir', 'out')
%             title(h(sInt+3), ['p=' num2str(ranksum(novSS,novGS))]);            
        end
        xlabel(h(2), 'Time to stim onset (s)')
        ylabel(t, [fdNames{f} '  change in firing rate from baseline'])
        xlim(h(1:3), [-1,2]); linkaxes(h(1:3),'xy'); linkaxes(h(4:6),'xy');
        xlim(h(4:6), [0.5 2.5]);
        suptitle([celltypes{ct} ': NovShamStim vs NovGoalStim'])
        saveas(t, fullfile(figdir,['Figure02_DE_' celltypes{ct} '.svg']));
    end
    
end
%% save figure
makefigurepretty(gcf)
figname = 'Figure02_DE';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')