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

%% Fig. 1F, left (raw firing rate)
fig = figure('units','inch','position',[0 0 8 4]);
t = tiledlayout('flow');

%load the most updated distance-based residual firing rate map across
%distance relative to reward zone (nUnits x nBins)
filename = getlatestfile_with_string(dirs.data2load, 'allsess_raw_vs_residuals_lap.mat');
load(fullfile(dirs.data2load, filename));

%define celltype
cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
    1:length(celltypes),'UniformOutput',false);

colormap(gray);
for ct = 1
    %min-max normalization per unit
    vMax = max(allsess_mean_fam_raw_all(cellT{ct},:),[],2);
    vMin = min(allsess_mean_fam_raw_all(cellT{ct},:),[],2);
    normMap_fam = (allsess_mean_fam_raw_all(cellT{ct},:) - vMin) ./ (vMax - vMin);
    
    nexttile;
    % temp = normMap_fam(ismember(allsess_sessinfo(cellT{ct},1), params.WTmice),:); %select WT units only
    temp = normMap_fam;
    temp = temp(all(~isnan(temp),2),:);
    if strcmp(celltypes{ct}, 'Pyramidal Cell')
        [~,maxLoc] = max(temp,[],2,'omitnan'); [~,b2_raw] = sort(maxLoc); %sort pyramidal cells by peak location
    else
        [~,minLoc] = min(temp,[],2,'omitnan'); [~,b2_raw] = sort(minLoc); %sort interneurons by min location
    end
    imagesc(temp(b2_raw,:));box off;
    if ct == length(celltypes); cb = colorbar('Location','eastoutside'); cb.Label.String = 'Normalized residual firing rate'; end
    if ct == 1; ylabel('Unit index'); end
    if ct == 2; xlabel('Distance to familiar RZ (deg)'); end
    set(gca, 'XTick',[1,ceil(length(position_binEdges)/2-1),length(position_binEdges)-1],...
        'XTickLabel',position_binEdges([1,ceil(length(position_binEdges)/2)-1,length(position_binEdges)-1]));
    title({celltypeNames{ct}, 'Raw firing rate'});
end
%% Fig. 1F, right (residual firing rate)
for ct = 1
    %min-max normalization per unit
    vMax = max(allsess_mean_fam_residual_all(cellT{ct},:),[],2);
    vMin = min(allsess_mean_fam_residual_all(cellT{ct},:),[],2);
    normMap_fam = (allsess_mean_fam_residual_all(cellT{ct},:) - vMin) ./ (vMax - vMin);
    
    nexttile;
    temp = normMap_fam;
    temp = temp(all(~isnan(temp),2),:);

    imagesc(temp(b2_raw,:));box off;
    if ct == length(celltypes); cb = colorbar('Location','eastoutside'); cb.Label.String = 'Normalized residual firing rate'; end
    if ct == 1; ylabel('Unit index'); end
    if ct == 2; xlabel('Distance to familiar RZ (deg)'); end
    set(gca, 'XTick',[1,ceil(length(position_binEdges)/2-1),length(position_binEdges)-1],...
        'XTickLabel',position_binEdges([1,ceil(length(position_binEdges)/2)-1,length(position_binEdges)-1]));
    title({celltypeNames{ct}, 'Residual firing rate'});
end

%% save figure
makefigurepretty(gcf)
figname = 'Figure01_F';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

