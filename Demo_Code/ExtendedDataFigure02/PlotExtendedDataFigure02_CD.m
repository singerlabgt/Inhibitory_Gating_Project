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

load(fullfile(maindir, "Demo_Data", 'allsess_raw_vs_residuals_lap.mat'))

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Wide Interneuron', 'Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron', 'WS Interneuron', 'Pyramidal Cell'}; %names used in our manuscript
cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
    1:length(celltypes),'UniformOutput',false);


%% compare raw vs residual per celltype
ax = arrayfun( @(x) subplot(2,3,x,'NextPlot','add','Box','off'), 1:6);
for ct = 1:length(celltypes)
    colormap(gray);
    vMax = max([allsess_mean_fam_raw_all(cellT{ct},:), allsess_mean_nov_raw_all(cellT{ct},:)],[],2);
    vMin = min([allsess_mean_fam_raw_all(cellT{ct},:), allsess_mean_nov_raw_all(cellT{ct},:)],[],2);
    
    temp = (allsess_mean_fam_raw_all(cellT{ct},:) - vMin) ./ (vMax - vMin);
    temp = temp(all(~isnan(temp),2),:); %remove NaN cells
    if strcmp(celltypes{ct}, 'Pyramidal Cell')
        [~,maxLoc] = max(temp,[],2); [~,b2] = sort(maxLoc); %sort pyramidal cells by peak location
    else
        [~,minLoc] = min(temp,[],2); [~,b2] = sort(minLoc); %sort interneurons by min location
    end
    
    imagesc(ax(ct),temp(b2,:));box off; 
    ylim(ax(ct), [1,length(b2)]);
    title(ax(ct), ['Raw: ' celltypes{ct}])
    set(ax(ct), 'XTick',[1,ceil(length(position_binEdges)/2-1),length(position_binEdges)-1],...
            'XTickLabel',position_binEdges([1,ceil(length(position_binEdges)/2)-1,length(position_binEdges)-1]));
    
    vMax = max([allsess_mean_fam_residual_all(cellT{ct},:), allsess_mean_nov_residual_all(cellT{ct},:)],[],2);
    vMin = min([allsess_mean_fam_residual_all(cellT{ct},:), allsess_mean_nov_residual_all(cellT{ct},:)],[],2);
    
    temp = (allsess_mean_fam_residual_all(cellT{ct},:) - vMin) ./ (vMax - vMin);    
    imagesc(ax(ct+3), temp(b2,:)); 
    ylim(ax(ct+3), [1,length(b2)]);
    title(ax(ct+3), ['Res: ' celltypes{ct}])
    set(ax(ct+3), 'XTick',[1,ceil(length(position_binEdges)/2-1),length(position_binEdges)-1],...
            'XTickLabel',position_binEdges([1,ceil(length(position_binEdges)/2)-1,length(position_binEdges)-1]));
        clearvars temp
end
set(gcf, 'Position', get(0, 'Screensize'),'PaperOrientation', 'landscape');
% save figure
makefigurepretty(gcf)
figname = 'ExtendedFigure02_C';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
%% ED Fig. D: variance explained by the model (new recording data added)
figure;
for ct = 1:length(celltypes)
    if startsWith(celltypes(ct), 'Narrow')
        colors = params.colors_narrowInt;
    elseif startsWith(celltypes(ct), 'Wide')
        colors = params.colors_wideInt;
    elseif startsWith(celltypes(ct), 'Pyr')
        colors = params.colors_pyr;
    end
    nexttile;
    h = histogram(allsess_stats(cellT{ct}, 1).*100,'BinEdges',0:5:55,'FaceColor',colors(3,:),'FaceAlpha',0.9);
    ylabel('Unit number')
    xlabel('Variance explained by the model (%)')
end
% save figure
makefigurepretty(gcf)
figname = 'ExtendedFigure02_D';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')