%% Jeong et al. 2023 MANUSCRIPT - FIGURE02_B,C
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

%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);
sessions = unique(allindex(:,1:2),'rows');


%% Fig. 2B-C: Speed-based performance, goalstim vs shamstim per stim intensity
fig = figure('units','inch','position',[0 0 6.5 4]);
t = tiledlayout(4,7,'TileSpacing','compact','Units','inches','OuterPosition',[0 0 6.5 4]);

behType = {'Speed','Lickrate','DeltaLickrate','LickLatency'};
bt = 1;
currBT = [behType{bt} 'ROC'];
fname = getlatestfile_with_string(dirs.data2load, 'getNovelBehaviorROC_Speed');
load(fullfile(dirs.data2load, fname));

clear ax
for ii = 1:6    
    if ii == 1
        ax(ii) = nexttile([2,2]); 
        auc_nov_stim = [ROC.nov_all.AUC]'; titleLbl = 'All trials';
    elseif ii == 2
        ax(ii) = nexttile([2,1]);
        auc_nov_stim = [ROC.nov_stim_none.AUC]'; titleLbl = 'No stim';
    elseif ii == 3
        ax(ii) = nexttile([2,1]);
        auc_nov_stim = [ROC.nov_stim_low.AUC]'; titleLbl = 'Low stim';
    elseif ii == 4
        ax(ii) = nexttile([2,1]);
        auc_nov_stim = [ROC.nov_stim_high.AUC]'; titleLbl = 'High stim';
    elseif ii == 5
        ax(ii) = nexttile([2,1]);
        auc_nov_stim = [ROC.nov_stim_all.AUC]'; titleLbl = 'L+H stim';
    else
        ax(ii) = nexttile([2,1]);
        auc_nov_stim = [ROC.nov_stim_short.AUC]'; titleLbl = 'Short stim';
    end
    
    %identify goal stim and sham stim session info: [animalID, recDay, novelDay]
    famgoalstim = unique(allindex(ismember(allindex(:,1), params.goalshamMice)...
        & allindex(:,6) == 1 & allindex(:,8) == 1 & allindex(:,9) == 1, [1:2,7]),'rows');
    novgoalstim = unique(allindex(ismember(allindex(:,1), params.goalshamMice)...
        & allindex(:,6) ~= 1 & allindex(:,8) == 1 & allindex(:,9) == 1, [1:2,7]),'rows');
    novshamstim = unique(allindex(ismember(allindex(:,1), params.goalshamMice)...
        & allindex(:,6) ~= 1 & allindex(:,8) == 1 & allindex(:,9) == 2, [1:2,7]),'rows');
    
    auc_nov_goalstim = auc_nov_stim(...
        ismember(uniqSess(:,1), novgoalstim(:,1)) & ismember(uniqSess(:,2), novgoalstim(:,2))...
        & ismember(uniqSess(:,3), novgoalstim(:,3)));
    % auc_nov_goalstim(13:15) = nan(3,1);
    auc_nov_shamstim = auc_nov_stim(...
        ismember(uniqSess(:,1), novshamstim(:,1)) & ismember(uniqSess(:,2), novshamstim(:,2))...
        & ismember(uniqSess(:,3), novshamstim(:,3)));
    % novgoalstim(13:15,:) = [61,220903,1;61,220904,2;61,220905,3]; %these three sessions were excluded bc off-target probe insertion  
    novgoalstim(:,4) = auc_nov_goalstim; novgoalstim = sortrows(novgoalstim,2); goalstimLbl = repmat({'Goal'},size(novgoalstim,1),1);
    novshamstim(:,4) = auc_nov_shamstim; shamstimLbl = repmat({'Sham'},size(novshamstim,1),1);
    auc_nov_goalstim = cell2mat(arrayfun( @(x) auc_nov_goalstim(novgoalstim(:,3)==x), 1:3,'UniformOutput',false));
    auc_nov_shamstim = cell2mat(arrayfun( @(x) auc_nov_shamstim(novshamstim(:,3)==x), 1:3,'UniformOutput',false));

    %plot in terms of percent change in AUC score from day 1
    auc_nov_goalstim = (auc_nov_goalstim - auc_nov_goalstim(:,1)) ./ auc_nov_goalstim(:,1) .*100;
    auc_nov_shamstim = (auc_nov_shamstim - auc_nov_shamstim(:,1)) ./ auc_nov_shamstim(:,1) .*100;

    hold on; box off; 
    jitter = 0.05; x_values_jittered = repmat([1:3], 5, 1) + jitter * randn(size(auc_nov_goalstim));
    scatter(ax(ii), reshape(x_values_jittered, [], 1), reshape(auc_nov_goalstim, [], 1), 15, params.colors_goalstim(3,:), 'filled', 'MarkerFaceAlpha', 0.5);

    jitter = 0.05; x_values_jittered = repmat([1:3], 5, 1) + jitter * randn(size(auc_nov_shamstim));
    scatter(ax(ii), reshape(x_values_jittered, [], 1), reshape(auc_nov_shamstim, [], 1), 15, params.colors_shamstim(3,:), 'filled', 'MarkerFaceAlpha', 0.5);

    errorbar(ax(ii), nanmean(auc_nov_goalstim),...
        nanstd(auc_nov_goalstim) ./ sqrt(size(auc_nov_goalstim,1)),...
        'Color',params.colors_goalstim(3,:),'LineWidth',2);
    errorbar(ax(ii), nanmean(auc_nov_shamstim),nanstd(auc_nov_shamstim) ./ sqrt(size(auc_nov_shamstim,1)),'Color',params.colors_shamstim(3,:),'LineWidth',2)
    xticks(ax(ii),1:3)
    title(ax(ii), titleLbl)
    hold off;     
end
ylabel(ax(1),'Area under the curve'); xlabel(ax(1),'Day in novel'); 
linkaxes(ax,'xy'); xlim(ax, [0.5 3.5]); xticks(ax, 1:3); 


%% save figure
makefigurepretty(gcf)
figname = 'Figure02_BC';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')