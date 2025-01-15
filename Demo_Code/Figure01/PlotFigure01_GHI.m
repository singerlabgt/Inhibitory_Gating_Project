clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameterss
[dirs, params] = getDefaultParameters(maindir); 
baseBins = 1:2; %first 2 bins to average across; to be used as baseline firing rate

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','Pyramidal Cell'}; %names used in our manuscript

 
fig = figure('units','inch','position',[0 0 6.5 4]);
t = tiledlayout(6,4,'TileSpacing','compact','Units','inches','OuterPosition',[0 0 6.5 4]);


%% Fig. 1G
%load the most updated distance-based residual firing rate map across
%distance relative to reward zone (nUnits x nBins)
filename = getlatestfile_with_string(dirs.data2load, 'allsess_raw_vs_residuals_distance2RZ.mat');
load(fullfile(dirs.data2load, filename));

%define celltype
cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
    1:length(celltypes),'UniformOutput',false);

colormap(gray);
populationAvg = cell(length(celltypes));
for ct = 1:length(celltypes)
    %min-max normalization per unit
    vMax = max(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
    vMin = min(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
    normMap_fam = (allsess_mean_fam_residual_hit(cellT{ct},:) - vMin) ./ (vMax - vMin);
    
    nexttile([3 2]);
    temp = normMap_fam(ismember(allsess_sessinfo(cellT{ct},1), params.WTmice),:); %select WT units only
    temp = temp(all(~isnan(temp),2),:);
    populationAvg{ct} = temp;
    if strcmp(celltypes{ct}, 'Pyramidal Cell')
        [~,maxLoc] = max(temp,[],2,'omitnan'); [~,b2] = sort(maxLoc); %sort pyramidal cells by peak location
    else
        [~,minLoc] = min(temp,[],2,'omitnan'); [~,b2] = sort(minLoc); %sort interneurons by min location
    end
    imagesc(temp(b2,:));box off;
    if ct == length(celltypes); cb = colorbar('Location','eastoutside'); cb.Label.String = 'Normalized residual firing rate'; end
    if ct == 1; ylabel('Unit index'); end
    if ct == 2; xlabel('Distance to familiar RZ (deg)'); end
    set(gca, 'XTick',[1,ceil(length(position_binEdges)/2-1),length(position_binEdges)-1],...
        'XTickLabel',position_binEdges([1,ceil(length(position_binEdges)/2)-1,length(position_binEdges)-1]));
    title(celltypeNames{ct});
end


%% Fig. 1I, left panel
ax = nexttile([3 2]);
for ct = 1:length(celltypes)
    if strcmp(celltypes{ct}, 'Pyramidal Cell')
        colors = params.colors_pyr;
    elseif strcmp(celltypes{ct}, 'Narrow Interneuron')
        colors = params.colors_narrowInt;
    end
    temp = populationAvg{ct};
    temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;
    ops.ax     = ax;
    ops.x_axis = mean(getBinEdges(position_binEdges),2);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;
    
    %t-test with Bonferroni correction for multiple comparisons
    h = nan(size(temp,2),1);
    p = nan(size(temp,2),1);
    stats = cell(size(temp,2),1);
    for iT = 1:size(temp,2)
        [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT));
    end
    sig = find(p < 0.05 / numel(p));
    arrayfun( @(ii) scatter(ax, position_binEdges(sig(ii)), 10+ct*2,...
        'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
end
xlim(ax,[-40 40])
set(ax, 'XTick',position_binEdges(find(ismember(position_binEdges,[-40,0,10,40]))),...
    'XTickLabel',{'-40','0','10','40'});
xline(ax, 0, 'k:'); xline(ax, 10, 'k:');
xlabel(ax, 'Distance to familiar RZ (deg)');
ylabel(ax, {'Change in residual firing','rate from baseline (%)'})


%% Fig. 1I, right panel
%load the most updated time-based residual firing rate map across time
%relative to reward zone entry (nUnits x nBins)
filename = getlatestfile_with_string(dirs.data2load, 'time2RZ.mat');
load(fullfile(dirs.data2load, filename));
time_binEdges = time_binEdges ./ params.samprate;

populationAvg = cell(length(celltypes));
for ct = 1:length(celltypes)
    units2incl = intersect(cellT{ct}, find(ismember(allsess_sessinfo(:,1), params.WTmice)));
    vMax = max(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
    vMin = min(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
    normMap_fam = (allsess_mean_fam_residual_hit(cellT{ct},:) - vMin) ./ (vMax - vMin);
    temp = normMap_fam(ismember(allsess_sessinfo(cellT{ct},1), params.WTmice),:); %select WT units only
    populationAvg{ct} = normMap_fam(all(~isnan(temp),2),:);
end

ax = nexttile([3 2]);
for ct = 1:length(celltypes)
    if strcmp(celltypes{ct}, 'Pyramidal Cell')
        colors = params.colors_pyr;
    elseif strcmp(celltypes{ct}, 'Narrow Interneuron')
        colors = params.colors_narrowInt;
    end
    temp = populationAvg{ct};
    temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;
    ops.ax     = ax;
    ops.x_axis = mean(getBinEdges(time_binEdges),2);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;
    
    %t-test with Bonferroni correction for multiple comparisons
    h = nan(size(temp,2),1);
    p = nan(size(temp,2),1);
    stats = cell(size(temp,2),1);
    for iT = 1:size(temp,2)
        [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT));
    end
    sig = find(p < 0.05 / numel(p));
    arrayfun( @(ii) scatter(ax, time_binEdges(sig(ii)), 10+ct*2,...
        'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
end
xlim(ax,[time_binEdges(1) time_binEdges(end)]);
set(ax, 'XTick',[time_binEdges(1),0, time_binEdges(end)],...
    'XTickLabel',time_binEdges([1,ceil(length(time_binEdges)/2)-1,length(time_binEdges)-1]));
xline(ax, 0, 'k:'); xline(ax, 10, 'k:');
xlabel(ax, 'Time to familiar RZ (s)');

% save figure
makefigurepretty(gcf)
figname = 'Figure01_GI';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
%% Fig. 1H
%load the most updated distance-based residual firing rate map across
%distance relative to reward zone (nUnits x nBins)
filename = getlatestfile_with_string(dirs.data2load, 'allsess_raw_vs_residuals_distance2RZ.mat');
load(fullfile(dirs.data2load, filename));

%define celltype
cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
    1:length(celltypes),'UniformOutput',false);

colormap(gray);
populationAvg = cell(length(celltypes));
climsfr = {[0 4], [5 25] [20 35], [20 70]};
%climsfr = {[0 4], [5 25], [20 50]};
for ct = [1]
    %min-max normalization per unit
    vMax = max(allsess_mean_fam_raw_hit(cellT{ct},:),[],2);
    vMin = min(allsess_mean_fam_raw_hit(cellT{ct},:),[],2);
    normMap_fam = (allsess_mean_fam_raw_hit(cellT{ct},:) - vMin) ./ (vMax - vMin);
    normMap_fam = allsess_mean_fam_raw_hit(cellT{ct},:);
    normMap_fam = normMap_fam(ismember(allsess_sessinfo(cellT{ct},1), params.WTmice),:); %select WT units only

    %find Inf values and replace with NaN
    normMap_fam(isinf(normMap_fam)) = NaN;
    avgMap = mean(normMap_fam,2,'omitnan'); 
    medMap = median(avgMap);
    quanMaps = quantile(avgMap,3);
    binMap = [0 5 25 max(avgMap)];
    binMap = [min(avgMap) quanMaps max(avgMap)];
    iM = discretize(avgMap,binMap);

    for b = 1:length(unique(iM))
        inclCell = iM == b;
        temp = normMap_fam(inclCell,:);
        populationAvg{ct} = temp(all(~isnan(temp),2),:);
    end
end

figure; ax = nexttile;
for ct = [1]
    if strcmp(celltypes{ct}, 'Pyramidal Cell')
        colors = params.colors_pyr;
    elseif strcmp(celltypes{ct}, 'Narrow Interneuron')
        colors = params.colors_narrowInt;
    elseif strcmp(celltypes{ct}, 'Wide Interneuron')
        colors = params.colors_wideInt;
    end
    temp = populationAvg{ct};
    ops.ax     = ax;
    ops.x_axis = mean(getBinEdges(position_binEdges),2);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;

end
xlim(ax,[-40 40])
set(ax, 'XTick',position_binEdges(find(ismember(position_binEdges,[-40,0,10,40]))), ...
    'XTickLabel',{'-40','0','10','40'});
xline(ax, 0, 'k:'); xline(ax, 10, 'k:');
xlabel(ax, 'Distance to familiar RZ (deg)');
ylabel(ax, {'Firing rate (Hz)'})

makefigurepretty(gcf)
figname = 'Figure01_H';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

