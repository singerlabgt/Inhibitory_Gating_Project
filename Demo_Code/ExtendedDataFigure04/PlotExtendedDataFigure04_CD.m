%% compare the same NS neurons' activity in reward trials w/o SWR
%% NOTE: use my matlab data path. need processing raw data

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
figdir = fullfile(maindir, 'Demo_Figures');
load(fullfile(maindir, 'Demo_Data', 'firingrate_wo_SWR')) % per trial NS interneuron firing rate w/o SWR events
% plot
position_binEdges = -60:2:70;
baseBins = 1:2;
colors = hex2rgb({'#d2dbec','#7894c5','#36489b'});

figure('Position',[440 486 644 212])
t = tiledlayout(1, 2);
ax = nexttile;
vMax = max(FR_noSWR, [], 2);
vMin = min(FR_noSWR, [], 2);
temp = (FR_noSWR - vMin) ./ (vMax - vMin); 
temp = (temp - nanmean(temp(:, baseBins), 2)) .* 100;
ops.ax     = ax;
ops.x_axis = mean(getBinEdges(position_binEdges),2);
ops.color_area = colors(end,:);
ops.color_line = colors(end,:);
ops.alpha      = 0.2;
ops.line_width = 2;
ops.error      = 'sem';
plot_areaerrorbar(temp, ops); hold on; box off;
h = nan(size(temp,2),1);
p = nan(size(temp,2),1);
stats = cell(size(temp,2),1);
for iT = 1:size(temp,2)
    [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT));
end
sig = find(p < 0.05 / numel(p));
arrayfun( @(ii) scatter(ax, position_binEdges(sig(ii)), 10,...
    'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
title('Trials without SWRs')
xlim(ax,[-40 40])
set(ax, 'XTick',position_binEdges(find(ismember(position_binEdges,[-40,0,10,40]))),...
    'XTickLabel',{'-40','0','10','40'});
xline(ax, 0, 'k:'); xline(ax, 10, 'k:');
xlabel(ax, 'Distance to familiar/novel RZ (deg)');
ylabel(ax, {'Change in residual firing','rate from baseline (%)'})
ylim([-25 20])

ax = nexttile;
vMax = max(FR_withSWR, [], 2);
vMin = min(FR_withSWR, [], 2);
temp = (FR_withSWR - vMin) ./ (vMax - vMin); 
temp = (temp - nanmean(temp(:, baseBins), 2)) .* 100;
ops.ax     = ax;
ops.x_axis = mean(getBinEdges(position_binEdges),2);
ops.color_area = colors(end,:);
ops.color_line = colors(end,:);
ops.alpha      = 0.2;
ops.line_width = 2;
ops.error      = 'sem';
plot_areaerrorbar(temp, ops); hold on; box off;
h = nan(size(temp,2),1);
p = nan(size(temp,2),1);
stats = cell(size(temp,2),1);
for iT = 1:size(temp,2)
    [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT));
end
sig = find(p < 0.05 / numel(p));
arrayfun( @(ii) scatter(ax, position_binEdges(sig(ii)), 10,...
    'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
title('Trials with SWRs')
xlim(ax,[-40 40])
set(ax, 'XTick',position_binEdges(find(ismember(position_binEdges,[-40,0,10,40]))),...
    'XTickLabel',{'-40','0','10','40'});
xline(ax, 0, 'k:'); xline(ax, 10, 'k:');
xlabel(ax, 'Distance to familiar/novel RZ (deg)');
ylabel(ax, {'Change in residual firing','rate from baseline (%)'})
ylim([-25 20])
yticks([-25 2])
makefigurepretty(gcf,1)

%% save figure
figname = 'ExtendedDataFigure04_CD';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

