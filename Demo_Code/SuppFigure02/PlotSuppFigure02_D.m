% simply use matlab default boxplot for firing rate comparison 
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

%% plot
lightsensitive_units = setdiff(cell_metrics.groundTruthClassification.lightsensitive, cell_metrics.tags.Bad);
ns_units = setdiff(find(strcmp(cell_metrics.putativeCellType, 'Narrow Interneuron')),cell_metrics.tags.Bad);

celltype_names = [repelem("Light-sensitive", length(lightsensitive_units)), repelem("NS int.", length(ns_units))];
firingrates = [cell_metrics.firingRate(lightsensitive_units), cell_metrics.firingRate(ns_units)];
boxplot(firingrates, celltype_names)
ylabel('Firing rate (spikes/s)')

makefigurepretty(gcf)
figname = 'SuppFigure02_D';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
