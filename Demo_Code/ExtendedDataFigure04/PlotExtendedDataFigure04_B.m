
clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 
params.colors_pvbc = [86 194 175]./255;
params.colors_aac = [222 64 182]./255;
params.colors_aac = [222 64 182]./255;
params.colors_cck = [247 196 46]./255;
params.colors_bistratified = [126 86 31]./255;
params.colors_ungrouped_ns_pv = [122 122 122]./255;

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end
datadir = fullfile(maindir, 'Demo_Data'); 
load(fullfile(datadir, 'SpikeThetaPhase.mat')); % NeuronFile containing theta phase and TS.
load(fullfile(datadir, 'FileList.mat'));
load(fullfile(datadir, 'cell_metrics.mat'))
load(fullfile(datadir, 'allsess_shuffledSig_distance2RZ.mat'), 'allsess_sigDecrease_fam');
% load('Y:\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023\Data\allsess_shuffledSig_distance2RZ.mat')
NeuronFile = NeuronFile(1:42); % only load the WT animals 
FileList = FileList(1:42);

allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
[sessions, ind] = unique(allindex(:,1:2), 'rows', 'stable'); %define session as one date
identifier = repelem({'N'}, size(allindex, 1)); identifier(allindex(:, 1) == 4) = {'X'}; %identifier for mouse
ident_peranimal = repelem({'N'}, size(sessions, 1)); ident_peranimal(sessions(:, 1) == 4) = {'X'}; %identifier for mouse

SavePath = fullfile(maindir, 'Demo_Figures', 'ExtendedDataFigure04_B_pvsubgrouping'); if ~exist(SavePath, 'dir'); mkdir(SavePath); end
if ~exist(fullfile(SavePath, 'phase2theta'), 'dir')
    mkdir(fullfile(SavePath, 'phase2theta'))
end
if ~exist(fullfile(SavePath, 'phase2theta', 'diff_thr'), 'dir')
    mkdir(fullfile(SavePath, 'phase2theta', 'diff_thr'))
end
if ~exist(fullfile(SavePath, 'phase2theta', 'spike_thr'), 'dir')
    mkdir(fullfile(SavePath, 'phase2theta', 'spike_thr'))
end
if ~exist(fullfile(SavePath, 'fr2swr'), 'dir')
    mkdir(fullfile(SavePath, 'fr2swr'))
end



%% phase to theta -- NOW use t-stat method. Previously used curve fitting (phasehist_curvefitting.m)
if ~exist(fullfile(SavePath, 'phase2theta', 'tstat'))
    mkdir(fullfile(SavePath, 'phase2theta', 'tstat'))
end
universal_tol = 0;
scale_factor = 1;
binwidth = 18; 
% arguments used in phasehist_test
spike_thr = 30;
qval_thr = 0.05; 
tstat_thr = 0.5;
apply_tstat_thr = 1; % apply tstat threshold or not
params.colors_pvbc = [86 194 175]./255;
params.colors_aac = [222 64 182]./255;
params.colors_aac = [222 64 182]./255;
params.colors_cck = [247 196 46]./255;
params.colors_bistratified = [126 86 31]./255;
params.colors_ungrouped_ns_pv = [122 122 122]./255;
pyr_curation = 1;
if pyr_curation == 1
    phasehist_pyr;
end
%%
phasehist_ttest;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FR to ripple

load(fullfile(datadir, 'ThetaSWRCellType_pyr_corr.mat'));
% load('allsess_sigDecrease_fam');

%% plot FR change around RZ with heatmap and line plot
%% Dist2RZ
[allsess_unitType, allsess_mean_fam_residual_hit, allsess_sessinfo, position_binEdges] = multipleLinearRegression_Dist2RZ(SavePath, cell_metrics);

% include overlapped data only -- between neuron file and nuri's manuscript
% allsess data
uniqSess_neuronfile = [[FileList.Subj]',[FileList.Date]'];
idx_inNF = find(ismember(allsess_sessinfo(:, 1:2), uniqSess_neuronfile, 'rows'));
idx_NS = find(ismember(allsess_unitType, "Narrow Interneuron"));
idx_inNF_NS = intersect(idx_inNF, idx_NS);% in neuron file, narrow spiking interneurons. 
allsess_sessinfo = allsess_sessinfo(idx_inNF_NS, :);
allsess_unitType = allsess_unitType(idx_inNF_NS, :);
allsess_mean_fam_residual_hit = allsess_mean_fam_residual_hit(idx_inNF_NS, 5:end-5); % Xiao change on 10/25/2023. original: allsess_mean_fam_residual_hit = allsess_mean_fam_residual_hit(idx_inNF_NS, :); 
position_binEdges = position_binEdges(5:end-5); % Xiao added on 10/25/2023

% neuron type filtering in neuronfile -- not bad
neuronfile_unitType_theta = [NeuronFile.ThetaBasedCellType]';
neuronfile_unitType_theta = neuronfile_unitType_theta(arrayfun(@(x) ~isempty(neuronfile_unitType_theta{x}), 1:length(neuronfile_unitType_theta)));
neuronfile_unitType_swr = [NeuronFile.SWRBasedCellType]';
neuronfile_unitType_swr = neuronfile_unitType_swr(arrayfun(@(x) ~isempty(neuronfile_unitType_swr{x}), 1:length(neuronfile_unitType_swr)));
% for the ungrouped cell of each aspect, grab the cell type from another
% aspect. No need to run (if run also changes nothing in the result
% figures) if the reference aspect is theta while the ungrouped cells were
% reassigned by swr.
for iCell = 1:length(neuronfile_unitType_theta) 
    if neuronfile_unitType_swr{iCell} == ["ungrouped"]
        neuronfile_unitType_swr{iCell} = neuronfile_unitType_theta{iCell};
    end
end

swr2theta = 0; % first assign based on SWR then separate PVBC/BS based on theta; otherwise assign based on theta first then separate CCK/AAC/PVBC based on swr

if swr2theta == 1
    COUNT=0; 
    for iCell = 1:length(neuronfile_unitType_swr)
        if length(neuronfile_unitType_swr{iCell}) > 1
            if ~isempty(intersect(neuronfile_unitType_theta{iCell}, neuronfile_unitType_swr{iCell}))
                neuronfile_unitType_swr{iCell} = intersect(neuronfile_unitType_theta{iCell}, neuronfile_unitType_swr{iCell});
            else
                COUNT = COUNT+1;
            end
        end
    end
    neuronfile_unitType = neuronfile_unitType_swr;
else
    COUNT= 0;
    % for each theta multiassignment find the swr assignment
    for iCell = 1:length(neuronfile_unitType_swr)
        if length(neuronfile_unitType_theta{iCell}) > 1
            if ~isempty(intersect(neuronfile_unitType_theta{iCell}, neuronfile_unitType_swr{iCell}))
                neuronfile_unitType_theta{iCell} = intersect(neuronfile_unitType_theta{iCell}, neuronfile_unitType_swr{iCell});
            else
                COUNT = COUNT+1;
            end
        end
    end
    neuronfile_unitType = neuronfile_unitType_theta;
end
neuronfile_unitType = neuronfile_unitType(arrayfun(@(x) ~isempty(neuronfile_unitType{x}), 1:length(neuronfile_unitType)));
neuronfile_unitType_char = arrayfun(@(x) cell2mat(neuronfile_unitType{x}), 1:length(neuronfile_unitType), 'UniformOutput', false);
allsess_unitType = neuronfile_unitType_char;
summary(categorical(neuronfile_unitType_char))

% prepare data for plotting
baseBins = 1:2; %first 2 bins to average across; to be used as baseline firing rate

%load cell type info
celltypes = {'pvbc','aac','cck', 'bistratified', 'ungrouped_ns_pv'}; 
celltypeNames = {'PV+ Basket Cell','Axon-axonic Cell','CCK+ cells', 'Bistratified Cell', 'Ungrouped Cell'};

allsess_sigDecrease_fam_ref = allsess_sigDecrease_fam;
allsess_unitType_ref = allsess_unitType;
allsess_sessinfo_ref = allsess_sessinfo;
allsess_mean_fam_residual_hit_ref = allsess_mean_fam_residual_hit;

for sig = 0:1:2 % 1: only include the cells significantly decreasing fr around RZ; 
    if sig == 1 % only include cells with sig decrease
        allsess_sigDecrease_fam = allsess_sigDecrease_fam_ref(idx_inNF_NS, :);
        unit2incl = sum(allsess_sigDecrease_fam(:,11:16),2)>0; 
        disp([num2str(sum(unit2incl)) ' cells have sig decrease']) % 95 in 203 have sig decrease
        allsess_unitType = allsess_unitType_ref(:, find(unit2incl));
        allsess_sessinfo = allsess_sessinfo_ref(find(unit2incl), :);
        allsess_mean_fam_residual_hit = allsess_mean_fam_residual_hit_ref(find(unit2incl), :);
    elseif sig == 0 % only include cells without sig decrease
        allsess_sigDecrease_fam = allsess_sigDecrease_fam_ref(idx_inNF_NS, :);
        unit2incl = sum(allsess_sigDecrease_fam(:,11:16),2)<=0; 
        disp([num2str(sum(unit2incl)) ' cells dont have sig decrease']) % 108 in 203 don't have sig decrease
        allsess_unitType = allsess_unitType_ref(:, find(unit2incl));
        allsess_sessinfo = allsess_sessinfo_ref(find(unit2incl), :);
        allsess_mean_fam_residual_hit = allsess_mean_fam_residual_hit_ref(find(unit2incl), :);
    else % include all cells
        allsess_sigDecrease_fam = allsess_sigDecrease_fam_ref(idx_inNF_NS, :);
        unit2incl = sum(allsess_sigDecrease_fam(:,11:16),2)<=10000; 
        disp([num2str(sum(unit2incl)) ' cells have or dont have sig decrease']) % 108 in 203 don't have sig decrease
        allsess_unitType = allsess_unitType_ref(:, find(unit2incl));
        allsess_sessinfo = allsess_sessinfo_ref(find(unit2incl), :);
        allsess_mean_fam_residual_hit = allsess_mean_fam_residual_hit_ref(find(unit2incl), :);
    end

    for exact_match = 1%0:1
        for leave_one_out = 0%:1
            %define celltype
            if exact_match == 0 % including multi-assignment
                if leave_one_out == 0 % contain current celltype
                    cellT = arrayfun( @(x) find(contains(allsess_unitType, celltypes{x})),...
                        1:length(celltypes),'UniformOutput',false); 
                else 
                    cellT = arrayfun( @(x) find(~contains(allsess_unitType, celltypes{x})),...
                        1:length(celltypes),'UniformOutput',false);
                end
            else
                if leave_one_out == 0
                    cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
                        1:length(celltypes),'UniformOutput',false); % put into manuscript
                else
                    cellT = arrayfun( @(x) find(~strcmp(allsess_unitType, celltypes{x})),...
                        1:length(celltypes),'UniformOutput',false);
                end
            end
            t1 = PV_subgrouping_fig1G_heatmap(params, celltypes, celltypeNames, cellT, allsess_sessinfo, allsess_mean_fam_residual_hit, position_binEdges);
            if ~exist(fullfile(SavePath, 'dist2rz_fam'), 'dir'); mkdir(fullfile(SavePath, 'dist2rz_fam')); end
            makefigurepretty(gcf, 1)
            figname = ['HeatmapAroundRZ_sig' num2str(sig) 'exact_match' num2str(exact_match) 'leave_one_out' num2str(leave_one_out)];
            savefigALP(fullfile(SavePath, 'dist2rz_fam/'), figname, 'filetype', 'pdf')
            exportgraphics(t1, fullfile(SavePath, 'dist2rz_fam', ['HeatmapAroundRZ_sig' num2str(sig) 'exact_match' num2str(exact_match) 'leave_one_out' num2str(leave_one_out) '.png']), 'Resolution',300)
            
            populationAvg = cell(length(celltypes));
            for ct = 1:length(celltypes)
                %min-max normalization per unit
                vMax = max(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
                vMin = min(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
                normMap_fam = (allsess_mean_fam_residual_hit(cellT{ct},:) - vMin) ./ (vMax - vMin);
                temp = normMap_fam(ismember(allsess_sessinfo(cellT{ct},1), params.WTmice),:); %select WT units only
                populationAvg{ct} = temp(all(~isnan(temp),2),:);
            end
            t2 = PV_subgrouping_fig1H_residFRaroundRZ('dist', params, celltypes, celltypeNames, populationAvg, baseBins, position_binEdges);
            makefigurepretty(gcf, 1);
            figname = ['FRcurveAroundRZ_sig' num2str(sig) 'exact_match' num2str(exact_match) 'leave_one_out' num2str(leave_one_out)];
            savefigALP(fullfile(SavePath, 'dist2rz_fam/'), figname, 'filetype', 'pdf')
            exportgraphics(t2, fullfile(SavePath, 'dist2rz_fam', ['FRcurveAroundRZ_sig' num2str(sig) 'exact_match' num2str(exact_match) 'leave_one_out' num2str(leave_one_out) '.png']), 'Resolution',300)
        end
    end
end
% close all

