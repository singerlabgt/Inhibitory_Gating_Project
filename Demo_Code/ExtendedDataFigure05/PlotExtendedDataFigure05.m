clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end
%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

sec_before = 1; % as duration
sec_after = 1;
binSize= 0.05; % 50ms bin
edges = -1* sec_before : binSize : sec_after;
light_onset_bin = find(edges == 0);


%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);

celltypes = {'Narrow Interneuron','Wide Interneuron','Pyramidal Cell', 'light-sensitive Pyramidal Cell', 'PV Interneuron'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Int.','WS Int.','Pyr.', 'Light-sensitive Pyr.', 'Optotagged'}; %names used in our manuscript
% load data containing trial-averaged firing rates of all cell types
load(fullfile(maindir, 'Demo_Data', 'spikemat_allsess_allCT'), 'spikemat_allsess_allCT')
%% EDFig. 5BF: plot heatmap
fig = figure('units','inch','position',[0,0,7,8]);
t = tiledlayout(3, 2);

for iCT = [3 5] % start from pyramidal
    ax = nexttile; 
    % normalize
    temp = spikemat_allsess_allCT{iCT};
    if isempty(temp)
        continue
    end
    vMax = max(temp,[],2);
    vMin = min(temp,[],2);
    temp = (temp - vMin) ./ (vMax - vMin);
    
    if strcmp(celltypes{iCT}, 'Pyramidal Cell') || strcmp(celltypes{iCT}, 'light-sensitive Pyramidal Cell')
        [~,maxLoc] = max(temp,[],2,'omitnan'); [~,b2] = sort(maxLoc); %sort pyramidal cells by peak location
    else
        [~,minLoc] = min(temp,[],2,'omitnan'); [~,b2] = sort(minLoc); %sort interneurons by min location
    end
    
    imagesc(temp(b2,:), [0 0.5]);box off; hold on;
    if iCT == length(celltypes); cb = colorbar('Location','eastoutside'); cb.Label.String = 'Normalized firing rate'; end
    if iCT == 4; ylabel('Unit index'); end
    if iCT == 5; xlabel('Time to Light Onset (seconds)'); end
    set(gca, 'XTick',[1,ceil(length(edges)/2-1),length(edges)-1],...
        'XTickLabel',edges([1,ceil(length(edges)/2)-1,length(edges)-1]));
    title(celltypeNames{iCT});
end
%% EDFig. 5CG: plot firing rate scatter plot: 
% Y axis - after-stim mean fr in hz. X axis - before-stim mean fr in hz
baseBins = 1:2;
for iCT = [3 5] % start from pyramidal
    ax = nexttile; 
    % normalize
    temp = spikemat_allsess_allCT{iCT};
    if isempty(temp)
        continue
    end
    preFR = arrayfun(@(iCell) mean(temp(iCell, 1 : (light_onset_bin - 1))) / binSize, 1:size(temp, 1));
    postFR = arrayfun(@(iCell) mean(temp(iCell, light_onset_bin : end)) / binSize, 1:size(temp, 1));
    % ft = fittype('a*x', 'independent', 'x', 'dependent', 'y');
    % f = fit(preFR', postFR', ft);
    % p = plot(f, preFR, postFR); 
    if iCT ==3
        col = params.colors_pyr(end,:);
    else
        col = params.colors_narrowInt(end,:);
    end
    p = scatter(preFR, postFR, 'filled', 'MarkerFaceColor', col);
    xlim(ax,[0 max([preFR postFR])]);
    ylim(ax,[0 max([preFR postFR])]);
    hold on
    plot(ax, xlim, ylim, '--k'); hold off
    if iCT == 3; ylabel('Post-stim firing rate (Hz)'); end
    if iCT == 4; xlabel('Pre-stim firing rate (Hz)'); end
    set(legend(gca), 'Visible', 'off')
end
%% EDFig. 5DH: plot traces
for iCT = [3 5]
    ax = nexttile;
    temp = spikemat_allsess_allCT{iCT};
    if isempty(temp)
        continue
    end
    vMax = max(temp,[],2);
    vMin = min(temp,[],2);
    temp = (temp - vMin) ./ (vMax - vMin);
    temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;
    if iCT ==3
        col = params.colors_pyr(end,:);
        ttesttail = 'left';
    else
        col = params.colors_narrowInt(end,:);
        ttesttail = 'right';
    end
    ops.ax     = ax;
    ops.x_axis = edges(1:end-1) + binSize/2;
    ops.color_area = col;
    ops.color_line = col;
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;
    xline(0); yline(0); hold on

    baseline = nanmean(temp(:, 1 : (light_onset_bin - 1)), 2);
    h = nan(size(temp,2),1);
    p = nan(size(temp,2),1);
    stats = cell(size(temp,2),1);
    for iT = 1:size(temp,2)
        [h(iT),p(iT),~,stats{iT}] = ttest2(temp(:,iT), baseline, 'Tail', ttesttail);
    end
    sig = find(p < 0.05 / numel(p));
    arrayfun( @(ii) scatter(ax, edges(sig(ii)), 100,...
                        'Marker','|','MarkerEdgeColor', col),1:length(sig));
end
%% save figure
makefigurepretty(gcf)
figname = 'ExtendedDataFigure05';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')