%% Jeong et al. 2023 MANUSCRIPT - FIGURE01_J
% NJeong 03/28/2023

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%
%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%% Fig. 1D. ROC example
%get all indices
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.WTmice);
allindex = allindex(animal_idx,:); %filter based on animals to include 
allindex = allindex(animal_idx,:); %filter based on animals to include 
[sessions, sessID] = unique(allindex(:,[1:2, 7]),'rows'); %session = [animalID, recording date, novelty day]

% exp = [0; find(diff(allindex(:,7))==-2); length(allindex)];
% for ii = 1:length(exp)-1
%     sets(ii,:) = unique(allindex(exp(ii)+1:exp(ii+1),2));
% end
load(fullfile(dirs.data2load, 'getNovelBehaviorROC_Speed_XZ_250105.mat'))

%% ROC curves of an example animal
animal = 11;
recDates = [200130 200131 200201];
ind_uniqSess = find(uniqSess(:, 1) == animal & ismember(uniqSess(:, 2),recDates));
legend_lbl = [arrayfun( @(x) {['Fam Day' num2str(x)]}, 1:length(recDates)), ...
    arrayfun( @(x) {['Nov Day' num2str(x)]}, 1:length(recDates))];

ax = axes('NextPlot','add','Box','off','XLim',[0 1], 'YLim',[0 1]);
plot_per_set = arrayfun( @(x) ...
    plot(ROC.fam_all.X{x}, ROC.fam_all.Y{x}, ...
    'Color', params.colors_fam(x-3,:), 'LineWidth', 2), ...
    4:6 );

plot_per_set = arrayfun( @(x) ...
    plot(ROC.nov_all.X{x}, ROC.nov_all.Y{x}, ...
    'Color', params.colors_nov(x-3,:), 'LineWidth', 2), ...
    4:6 );
plot([0 1], [0 1], 'k--')

title(ax, [params.iden num2str(animal) ' ' num2str(recDates(1)) ' to ' num2str(recDates(end))]);
xlabel(ax, 'False positive rate')
ylabel(ax, 'True positive rate')
legend(ax, legend_lbl, 'box', 'off', 'Location', 'southeast')

%% save figure
makefigurepretty(gcf);
figname = 'Figure01_D';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% AUC across all WT mice
wtFam = uniqSess(ismember(uniqSess(:,1), params.WTmice),:);
wtNov = uniqSess(ismember(uniqSess(:,1), params.WTmice),:);
auc_fam_wt = [ROC.fam_all.AUC];
auc_fam_wt = auc_fam_wt(ismember(uniqSess(:,1), params.WTmice))';
wtFam(:,4) = auc_fam_wt; wtFamLbl = repmat({'Fam'},size(wtFam,1),1);
auc_fam_wt = cell2mat(arrayfun( @(x) auc_fam_wt(uniqSess(1:length(auc_fam_wt),3)==x), 1:3,'UniformOutput',false));

auc_nov_wt = [ROC.nov_all.AUC];
auc_nov_wt = auc_nov_wt(ismember(uniqSess(:,1), params.WTmice))';
wtNov(:,4) = auc_nov_wt; wtNovLbl = repmat({'Nov'},size(wtFam,1),1);
auc_nov_wt = cell2mat(arrayfun( @(x) auc_nov_wt(uniqSess(1:length(auc_nov_wt),3)==x), 1:3,'UniformOutput',false));
auc_nov_wt(9,:) = []; %remove clear outlier 

fig = figure('units','inch','position',[0 0 3 4]);
t = tiledlayout(4,2,'TileSpacing','compact','Units','inches','OuterPosition',[0 0 2.5 4]);
ax(1) = nexttile([2,1]);
hold on; box off;
errorbar(ax(1), nanmean(auc_nov_wt),nanstd(auc_nov_wt) ./ sqrt(size(auc_nov_wt,1)),'Color',params.colors_nov(2,:),'LineWidth',2)
errorbar(ax(1), nanmean(auc_fam_wt),nanstd(auc_fam_wt) ./ sqrt(size(auc_fam_wt,1)),'Color',params.colors_fam(3,:),'LineWidth',2)

%plot in % change in AUC from day 1
auc_fam_wt = (auc_fam_wt - auc_fam_wt(:,1)) ./ auc_fam_wt(:,1) .*100;
auc_nov_wt = (auc_nov_wt - auc_nov_wt(:,1)) ./ auc_nov_wt(:,1) .*100;     

ax(2) = nexttile([2,1]);
hold on; box off;
errorbar(ax(2), nanmean(auc_nov_wt),nanstd(auc_nov_wt) ./ sqrt(size(auc_nov_wt,1)),'Color',params.colors_nov(2,:),'LineWidth',2)
errorbar(ax(2), nanmean(auc_fam_wt),nanstd(auc_fam_wt) ./ sqrt(size(auc_fam_wt,1)),'Color',params.colors_fam(3,:),'LineWidth',2)

hold off;
title(t, {'Speed ROC', 'familiar (black) vs novel (green)'})
ylabel(ax(1),'Area under the curve'); ylabel(ax(2),'Change from Day 1 (%)'); 
xlabel(ax(1),'Day in novel'); 
linkaxes(ax,'x'); xlim(ax, [0.5 3.5]); xticks(ax, 1:3);
%% save figure
makefigurepretty(gcf);
figname = 'Figure01_E';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')