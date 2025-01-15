clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%all behavioral types to run 
behType = {'Speed','Lickrate','LickLatency'};

% get all indices 
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions

for bt = 1:length(behType)  
    fname = getlatestfile_with_string(fullfile(maindir, "Demo_Data/"), ['getNovelBehaviorROC_' behType{bt}]);
    load(fullfile(maindir, "Demo_Data/", fname));
    %% goalstim vs shamstim per stim intensity    
    for ii = 1:6
        if ii == 1
            ax(ii) = nexttile([2,2]);
            auc_nov_stim = [ROC.nov_all.AUC]'; titleLbl = 'AllTrials';        
        elseif ii == 2
            ax(ii) = nexttile([2,1]);
            auc_nov_stim = [ROC.nov_stim_low.AUC]'; titleLbl = 'LowStim';
        elseif ii == 3
            ax(ii) = nexttile([2,1]);
            auc_nov_stim = [ROC.nov_stim_high.AUC]'; titleLbl = 'HighStim';
        elseif ii == 4
            ax(ii) = nexttile([2,1]);
            auc_nov_stim = [ROC.nov_stim_none.AUC]'; titleLbl = 'NoStim';
        elseif ii == 5
            ax(ii) = nexttile([2,1]);
            auc_nov_stim = [ROC.nov_stim_all.AUC]'; titleLbl = 'BothStim';
        else
            ax(ii) = nexttile([2,1]);
            auc_nov_stim = [ROC.nov_stim_short.AUC]'; titleLbl = 'ShortStim'; %trials with stim <=5 seconds long only 
        end
        
        
        %identify goal stim and sham stim session info: [animalID, recDay, novelDay]
        goalshamMice = [57 61 63 62 65];
        famgoalstim = unique(allindex(ismember(allindex(:,1), goalshamMice)...
            & allindex(:,6) == 1 & allindex(:,8) == 1 & allindex(:,9) == 1, [1:2,7]),'rows');
        novgoalstim = unique(allindex(ismember(allindex(:,1), goalshamMice)...
            & allindex(:,6) ~= 1 & allindex(:,8) == 1 & allindex(:,9) == 1, [1:2,7]),'rows');
        novshamstim = unique(allindex(ismember(allindex(:,1), goalshamMice)...
            & allindex(:,6) ~= 1 & allindex(:,8) == 1 & allindex(:,9) == 2, [1:2,7]),'rows');
        

    end
    %% AUC in FAM with stim in PVxAi32 mice 
    auc_noVR = [ROC.noVR_all.AUC];
    auc_fam_stim_none = [ROC.fam_stim_none.AUC];
    auc_fam_stim_low = [ROC.fam_stim_low.AUC];
    auc_fam_stim_high = [ROC.fam_stim_high.AUC];
    auc_fam_stim_all = [ROC.fam_stim_all.AUC];
    auc_fam_all = [ROC.fam_stim_all.AUC];
    
    auc_fam_stim_low = auc_fam_stim_low(...
        ismember(uniqSess(:,1), famgoalstim(:,1)) & ismember(uniqSess(:,2), famgoalstim(:,2))...
        & ismember(uniqSess(:,3), famgoalstim(:,3)))';
    auc_fam_stim_high = auc_fam_stim_high(...
        ismember(uniqSess(:,1), famgoalstim(:,1)) & ismember(uniqSess(:,2), famgoalstim(:,2))...
        & ismember(uniqSess(:,3), famgoalstim(:,3)))';
    auc_fam_stim_all = auc_fam_stim_all(...
        ismember(uniqSess(:,1), famgoalstim(:,1)) & ismember(uniqSess(:,2), famgoalstim(:,2))...
        & ismember(uniqSess(:,3), famgoalstim(:,3)))';
    auc_fam_stim_none = auc_fam_stim_none(...
        ismember(uniqSess(:,1), famgoalstim(:,1)) & ismember(uniqSess(:,2), famgoalstim(:,2))...
        & ismember(uniqSess(:,3), famgoalstim(:,3)))';
    auc_fam_all = auc_fam_all(...
        ismember(uniqSess(:,1), famgoalstim(:,1)) & ismember(uniqSess(:,2), famgoalstim(:,2))...
        & ismember(uniqSess(:,3), famgoalstim(:,3)))';
    auc_noVR = auc_noVR(...
        ismember(uniqSess(:,1), famgoalstim(:,1)) & ismember(uniqSess(:,2), famgoalstim(:,2))...
        & ismember(uniqSess(:,3), famgoalstim(:,3)))';
    
    figure; ax = axes('NextPlot','add','Box','off','XLim',[0.5,6.5]);
    y = [auc_noVR,auc_fam_stim_none,auc_fam_stim_low,auc_fam_stim_high,auc_fam_stim_all,auc_fam_all];
    plot(ax, y' ,'Color',params.colors_fam(1,:))
    errorbar(ax, nanmean(y),nanstd(y) ./ sqrt(size(y,1)),'Color',params.colors_fam(3,:),'LineWidth',2)
    xticks(ax,1:length(y));
    xticklabels(ax,{'NoVR','NoStim','LowStim','HighStim','AllStim','AllTrials'})
    xlabel(ax,'Stim condition in Familiar')
    ylabel(ax,'Area under the curve')
    title(ax, [behType{bt}, ' ROC in Familiar'])
    %% save figure
    makefigurepretty(gcf);
    figname = 'SuppFigure05_C';
    savefigALP([figdir '/'], figname, 'filetype', 'pdf')
   
end

