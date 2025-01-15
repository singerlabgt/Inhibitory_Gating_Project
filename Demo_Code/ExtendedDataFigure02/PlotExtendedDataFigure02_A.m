%% Jeong et al. 2023 MANUSCRIPT - FIGURE 01_G,H,I
% NJeong 03/23/2023

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
celltypes = {'Narrow Interneuron','Wide Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','WS Interneuron','Pyramidal Cell'}; %names used in our manuscript

 
fig = figure('units','inch','position',[0 0 12 2]);
% hold on
tiledlayout(1,5,'TileSpacing','compact');


%% EDFig. 3A
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
    %inclFR = avgMap > 1; 
    %normMap_fam = normMap_fam(inclFR,:);
    %avgMap = avgMap(inclFR);

    medMap = median(avgMap);
    
    quanMaps = quantile(avgMap,3);
    binMap = [0 5 25 max(avgMap)];
    binMap = [min(avgMap) quanMaps max(avgMap)];
    iM = discretize(avgMap,binMap);

    for b = 1:length(unique(iM))
        inclCell = iM == b;
        temp = normMap_fam(inclCell,:);
        %temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;

        nexttile;
        populationAvg{ct} = temp(all(~isnan(temp),2),:);
        if strcmp(celltypes{ct}, 'Pyramidal Cell')
            [~,maxLoc] = max(temp,[],2,'omitnan'); [~,b2] = sort(maxLoc); %sort pyramidal cells by peak location
        else
            [~,minLoc] = min(temp,[],2,'omitnan'); [~,b2] = sort(minLoc); %sort interneurons by min location
        end
        imagesc(temp(b2,:), climsfr{b});box off;
        cb = colorbar('Location','eastoutside'); cb.Label.String = 'Firing rate (Hz)';
        if ct == length(celltypes); cb = colorbar('Location','eastoutside'); cb.Label.String = 'Normalized residual firing rate'; end
        if ct == 1; ylabel('Unit index'); end
        if ct == 2; xlabel('Distance to familiar RZ (deg)'); end
        set(gca, 'XTick',[1,ceil(length(position_binEdges)/2-1),length(position_binEdges)-1],...
            'XTickLabel',position_binEdges([1,ceil(length(position_binEdges)/2)-1,length(position_binEdges)-1]));
        title([ num2str(b) '-quartile']);
    end
    sgtitle(celltypeNames{ct});
end

makefigurepretty(gcf)
figname = 'ExtendedDataFigure02_A';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
