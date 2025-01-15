%% Jeong et al. 2023 MANUSCRIPT - FIGURE 03_A,B
% NJeong 03/23/2023

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path

%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%load cell type info
% load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','Pyramidal Cell'}; %names used in our manuscript

%% load the most updated distance-based residual firing rate map across
%distance relative to reward zone (nUnits x nBins)
filename = getlatestfile_with_string(dirs.data2load, 'allsess_raw_vs_residuals_stim_time2RZ.mat');
load(fullfile(dirs.data2load, filename));

%identify cell types for all units 
cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
    1:length(celltypes),'UniformOutput',false);

%load shuffled residual output data
filename = getlatestfile_with_string(dirs.data2load, 'shuffledSig');
load(fullfile(dirs.data2load, filename));


min_time = -4;
max_time = 4;
min_time_ind = find(time_binEdges == min_time*params.samprate);
max_time_ind = find(time_binEdges == max_time*params.samprate);
allsess_mean_nov_residual_hit_copy = allsess_mean_nov_residual_hit;
allsess_mean_nov_residual_hit = allsess_mean_nov_residual_hit(:, min_time_ind : max_time_ind);
time_binEdges_copy = time_binEdges;
time_binEdges = time_binEdges(min_time_ind : max_time_ind);
%identify cell types for all units 
% cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
%     1:length(celltypes),'UniformOutput',false);

%load shuffled residual output data
filename = getlatestfile_with_string(dirs.data2load, 'shuffledSig');
load(fullfile(dirs.data2load, filename));

% NOVEL DAY 1-3: plot WT sessions only, split into 3 days, per row=cell type
fig = figure('units','inch','position',[0,0,6.5,2.8]);
t = tiledlayout(1,3,'TileSpacing','compact');
for ii = 1:3
    ax(ii) = nexttile; hold on; box off;
for ct = 1:length(celltypes)
    
    %define color scheme per celltype
    if startsWith(celltypes(ct), 'Narrow')
        colors = params.colors_narrowInt;
        iInclUnits = allsess_sigDecrease_nov;
    elseif startsWith(celltypes(ct), 'Pyr')
        colors = params.colors_pyr;
        iInclUnits = allsess_sigIncrease_nov;
    end
    
    %define population min and max for normalization purposes
    vMax = max(allsess_mean_nov_residual_hit(cellT{ct},:),[],2);
    vMin = min(allsess_mean_nov_residual_hit(cellT{ct},:),[],2);
    normMap_nov = (allsess_mean_nov_residual_hit(cellT{ct},:) - vMin) ./ (vMax - vMin);
    colormap(gray);
    % for ii = 1:3
        sess2incl = ismember(allsess_sessinfo(cellT{ct},1), params.WTmice) & allsess_sessinfo(cellT{ct},3)==ii;
        temp = normMap_nov(sess2incl,:);
        sess2incl = sess2incl(all(~isnan(temp),2));
        temp = temp(all(~isnan(temp),2),:); %remove NaN cells
        
        %collect map from units with significant pos/neg slopes only
        unit2incl = sum(iInclUnits(cellT{ct},11:16),2)>0 & ismember(allsess_sessinfo(cellT{ct},1), params.WTmice) & allsess_sessinfo(cellT{ct},3)==ii;
        temp2 = normMap_nov(unit2incl,:);
        temp2 = temp2(all(~isnan(temp2),2),:); %remove NaN cells
        wt_nov{ct,ii} = temp2;
       
                
    % end
    % cb = colorbar;
    % cb.Layout.Tile = 'east';
     clear temp2;
    % for ii = 1:3
        temp = wt_nov{ct,ii};
        temp = (temp - nanmean(temp(:,1:2),2)) .* 100;
        wt_nov_delta{ct,ii} = (temp - nanmean(temp(:,1:2),2)) .* 100;
        temp2{ii} = nanmean(temp(:,11:12),2);
        
        ops.ax     = ax(ii);
        ops.x_axis = time_binEdges ./ params.samprate; %mean(getBinEdges(time_binEdges),2);
        ops.color_area = colors(3,:);
        ops.color_line = colors(3,:);
        ops.alpha      = 0.2;
        ops.line_width = 2;
        ops.error      = 'sem';
        plot_areaerrorbar(temp, ops); hold on;

        % add triangles to indicate persistent decrease: first bin followed by
        % 2.5 seconds
        for iT = 1 : size(temp, 2) - 24
            if sum(nanmean(temp(:, iT:iT+24), 1) < 0) == 25 && ct == 1
                scatter(ax(ii), time_binEdges(iT)/params.samprate, 50 - 10*ii, 30, ...
                    'Marker','^','MarkerEdgeColor',colors(ii,:), 'MarkerFaceColor',colors(ii, :))
                break 
            end
            if sum(nanmean(temp(:, iT:iT+24), 1) > 0) == 25 && ct == 2
                scatter(ax(ii), time_binEdges(iT)/params.samprate, 50 - 10*ii, 30, ...
                    'Marker','^','MarkerEdgeColor',colors(ii,:), 'MarkerFaceColor',colors(ii, :))
                break 
            end
        end
    % end
    xlim(ax(ii),[-4 4])
    ylim(ax(ii), [-25 40])
    set(ax(ii), 'XTick',[min_time,-2,0,2,max_time],...
        'XTickLabel',{'-4','-2','0','2', '4'}, 'TickDir', 'out');
    xline(ax(ii), 0, 'k:'); yline(ax(ii),0, 'k:');
end
title(['Novel Day ', num2str(ii)])
end
xlabel(ax(2),'Time to Novel RZ (sec)');
ylabel(ax(2),'Change in residual firing rate from baseline (%)');
linkaxes(ax,'xy');

makefigurepretty(gcf)
figname = 'ExtendedDataFigure07_A';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
