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
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron', 'Narrow Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'PV Interneuron', 'NS Interneuron','Pyramidal Cell'}; %names used in our manuscript

%load the most updated distance-based residual firing rate map across
%distance relative to reward zone (nUnits x nBins)
filename = getlatestfile_with_string(dirs.data2load, 'allsess_raw_vs_residuals_distance2RZ.mat');
load(fullfile(dirs.data2load, filename));

%identify cell types for all units 
cellT = arrayfun( @(x) find(strcmp(allsess_unitType, celltypes{x})),...
    1:length(celltypes),'UniformOutput',false);
% check if a cell is optotagged from all cells
total_metricID = 1:length(cell_metrics.sessionName);
lightsensitiveUnits = []; 
for iCell = 1:size(allsess_sessinfo, 1)
    iden = 'N'; if allsess_sessinfo(iCell, 1) == 4; iden = 'X'; end
    sessstr = [iden num2str(allsess_sessinfo(iCell, 1)) '_' num2str(allsess_sessinfo(iCell, 2)) '_CA3'];
    cluID = allsess_unitID(iCell);
    curr_metricID = total_metricID(find(ismember(cell_metrics.sessionName, sessstr) & ismember(cell_metrics.cluID, cluID)));
    if length(curr_metricID)~=1
        disp('error')
    end
    if ismember(curr_metricID, cell_metrics.groundTruthClassification.lightsensitive)
        lightsensitiveUnits = [lightsensitiveUnits, iCell];
    end
end
cellT{1} = intersect(cellT{1}, lightsensitiveUnits);
%load shuffled residual output data
filename = getlatestfile_with_string(dirs.data2load, 'shuffledSig');
load(fullfile(dirs.data2load, filename));

noRZstimsess = (allsess_sessinfo(:, 1) == 52 & allsess_sessinfo(:, 2) < 211228) | ...
    (allsess_sessinfo(:, 1) == 53 & allsess_sessinfo(:, 2) < 220217) | ...
    (allsess_sessinfo(:, 1) == 54 & allsess_sessinfo(:, 2) < 220211) | ...
    (allsess_sessinfo(:, 1) == 57 & allsess_sessinfo(:, 2) < 220405) | ...
    (allsess_sessinfo(:, 1) == 61 & allsess_sessinfo(:, 2) < 220908) | ...
    (allsess_sessinfo(:, 1) == 63 & allsess_sessinfo(:, 2) < 221009); % stim, but not in fam RZ (can include NRZ stim)


%% plot heatmap for PV cells only
fig = figure('units','inch','position',[0,0,6.5 2.5]);
t = tiledlayout(1, 2,'TileSpacing','compact');

% sig decreased cells
nexttile;
colormap(gray);
ct = 1;

vMax = max(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
vMin = min(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
normMap_fam = (allsess_mean_fam_residual_hit(cellT{ct},:) - vMin) ./ (vMax - vMin);

sess2incl = ismember(cellT{ct}, find(noRZstimsess)) ;
temp = normMap_fam(sess2incl,:);
sess2incl = sess2incl(all(~isnan(temp),2));
temp = temp(all(~isnan(temp),2),:); %remove NaN cells


if strcmp(celltypes{ct}, 'Pyramidal Cell')
    [~,maxLoc] = max(temp,[],2); [~,b2] = sort(maxLoc); %sort pyramidal cells by peak location
else
    [~,minLoc] = min(temp,[],2); [~,b2] = sort(minLoc); %sort interneurons by min location
end
imagesc(temp(b2,:));box off;
set(gca, 'XTick',[1,ceil(length(position_binEdges)/2-1),length(position_binEdges)-1],...
    'XTickLabel',position_binEdges([1,ceil(length(position_binEdges)/2)-1,length(position_binEdges)-1]));
xlabel('Distance to FAM RZ (deg)');
ylabel(celltypeNames{ct});
title('PV Interneurons')


% all opto-tagged cells
colors = params.colors_narrowInt;
ct = 1; % only analyze PV cells here
%define population min and max for normalization purposes
vMax = max(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
vMin = min(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
normMap_fam = (allsess_mean_fam_residual_hit(cellT{ct},:) - vMin) ./ (vMax - vMin);
    

%collect map from units with significant pos/neg slopes only
if startsWith(celltypeNames(ct), 'PV')
    % only include PV cells decreasing at RZ
    unit2incl = ismember(cellT{ct}, find(noRZstimsess));
end
temp2 = normMap_fam(unit2incl,:);
temp2 = temp2(all(~isnan(temp2),2),:); %remove NaN cells
temp = temp2;
temp = (temp - nanmean(temp(:,1:2),2)) .* 100;
ax = nexttile;
ops.ax     = ax(ct);
ops.x_axis = mean(getBinEdges(position_binEdges),2);
ops.color_area = colors(end,:);
ops.color_line = colors(end,:);
ops.alpha      = 0.2;
ops.line_width = 2;
ops.error      = 'sem';
plot_areaerrorbar(temp, ops); hold on;
disp(num2str(size(temp)))

%t-test with Bonferroni correction for multiple comparisons
h = nan(size(temp,2),1);
p = nan(size(temp,2),1);
stats = cell(size(temp,2),1);
for iT = 1:size(temp,2)
    [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT));
end
sig = find(p < 0.05 / numel(p));
arrayfun( @(ii) scatter(ax(ct), position_binEdges(sig(ii)), 25,...
    'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
xlim(ax(ct),[-40 40])
set(ax(ct), 'XTick',position_binEdges(find(ismember(position_binEdges,[-40,0,10,40]))),...
    'XTickLabel',{'-40','0','10','40'});
xline(ax(ct), 0, 'k:'); xline(ax(ct), 10, 'k:'); yline(ax(ct),0, 'k:');

ylabel(ax(ct),{'Change in residual firing', 'rate from baseline (%)'});
xlabel(ax(ct), 'Distance to FAM RZ (deg)');
linkaxes(ax,'xy');
title('all PV cells')
makefigurepretty(gcf)
figname = 'Figure01_L';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')