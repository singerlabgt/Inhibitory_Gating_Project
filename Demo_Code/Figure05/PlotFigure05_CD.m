%% Jeong et al. 2023 MANUSCRIPT - FIGURE05_CD
% XZheng 10/01/2024

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

load(fullfile(dirs.data2load, 'allsess_mean_residual_time2rz.mat'));
%% plot 
fig = figure('units','inch','position',[1.3229 0.8229 6.1875 4.4896]);
t = tiledlayout(2, 2,'TileSpacing','compact');

colors = hex2rgb({'#d2dbec','#7894c5','#36489b'});
ax = nexttile;
baseBins = 1:2;
vMax = max(allsess_mean_residual_hit, [], 2);
vMin = min(allsess_mean_residual_hit, [], 2);
temp = (allsess_mean_residual_hit - vMin) ./ (vMax - vMin); 
normMap_RZ = temp;
%%% heatmap
temp = normMap_RZ(all(~isnan(normMap_RZ),2),:); %remove NaN cells
colormap("gray")
[~,minLoc] = min(temp,[],2); [~,b2] = sort(minLoc); %sort interneurons by min location
imagesc(timebin_edges, 1:size(temp,1), temp(b2,:));box off;
colorbar
set(gca, 'XTick', [-6 0 2 6]);
xlabel('Time to RZ (s)');
title('PV Interneurons')


ax = nexttile; 
temp = (temp - nanmean(temp(:, baseBins), 2)) .* 100;
ops.ax     = ax;
ops.x_axis = timebin_edges(1:end-1);
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
arrayfun( @(ii) scatter(ax, timebin_edges(sig(ii)), 10,...
    'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
title('Update task')
xlim(ax,[-6 6])
set(ax, 'XTick',timebin_edges(find(ismember(timebin_edges,[-6,0,2,6]))),...
    'XTickLabel',{'-6','0','2','6'});
xline(ax, 0, 'k:'); xline(ax, 10, 'k:');
xlabel(ax, 'Time to RZ (s)');
ylabel(ax, {'Change in residual firing','rate from baseline (%)'})
ylim([-30 20])

colors = hex2rgb({'#d2dbec','#7894c5','#36489b'});
baseBins = 1:2;
vMax = max(allsess_mean_residual_updatecue, [], 2);
vMin = min(allsess_mean_residual_updatecue, [], 2);
temp = (allsess_mean_residual_updatecue - vMin) ./ (vMax - vMin); 
normMap_update = temp;
%%% heatmap
ax = nexttile; 
temp = normMap_update(all(~isnan(normMap_update),2),:); %remove NaN cells
colormap("gray")
[~,minLoc] = min(temp,[],2); [~,b2] = sort(minLoc); %sort interneurons by min location
imagesc(timebin_edges, 1:size(temp,1), temp(b2,:));box off;
colorbar
set(gca, 'XTick', [-6 0 2 6]);
xlabel('Time to Update Cue (s)');
title('PV Interneurons')

ax = nexttile;
temp = (temp - nanmean(temp(:, baseBins), 2)) .* 100;
ops.ax     = ax;
ops.x_axis = timebin_edges(1:end-1);
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
arrayfun( @(ii) scatter(ax, timebin_edges(sig(ii)), 10,...
    'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
title('Update task')
xlim(ax,[-6 6])
set(ax, 'XTick',timebin_edges(find(ismember(timebin_edges,[-6,0,2,6]))),...
    'XTickLabel',{'-6','0','2','6'});
xline(ax, 0, 'k:'); xline(ax, 10, 'k:');
xlabel(ax, 'Time to Update Cue (s)');
ylabel(ax, {'Change in residual firing','rate from baseline (%)'})
ylim([-30 20])

%% save figure
makefigurepretty(gcf)
figname = 'Figure05_CD';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

