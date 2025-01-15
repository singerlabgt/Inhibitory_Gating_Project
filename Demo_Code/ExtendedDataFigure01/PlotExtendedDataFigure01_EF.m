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
%get all indices
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.animals);
allindex = allindex(animal_idx,:); %filter based on animals to include 
allindex = allindex(animal_idx,:); %filter based on animals to include 
sessions = unique(allindex(:,[1:2,7]), 'rows'); %define session as one date


%% bar plot speed ROC across days comparison between goal stim vs sham stim
behType = {'Speed', 'LickLatency', 'Lickrate'};
for bt = 1:length(behType)
    currBT = behType{bt};
    fname = getlatestfile_with_string(fullfile(maindir, 'Demo_Data'), ['getNovelBehaviorROC_' currBT]);
    load(fullfile(maindir, 'Demo_Data', fname));
    
    %% goalstim vs shamstim per stim intensity
    %identify goal stim and sham stim session info: [animalID, recDay, novelDay]
    famgoalstim = unique(allindex(ismember(allindex(:,1), params.goalshamMice)...
        & allindex(:,6) == 1 & allindex(:,8) == 1 & allindex(:,9) == 1, [1:2,7]),'rows');
    novgoalstim = unique(allindex(ismember(allindex(:,1), params.goalshamMice)...
        & allindex(:,6) ~= 1 & allindex(:,8) == 1 & allindex(:,9) == 1, [1:2,7]),'rows');
    novshamstim = unique(allindex(ismember(allindex(:,1), params.goalshamMice)...
        & allindex(:,6) ~= 1 & allindex(:,8) == 1 & allindex(:,9) == 2, [1:2,7]),'rows');
        
    auc_nov_stim_none = [ROC.nov_stim_none.AUC]';
    auc_nov_stim = [ROC.nov_stim_low.AUC]';
    auc_nov_stim_high = [ROC.nov_stim_high.AUC]';
    auc_nov_stim_all = [ROC.nov_stim_all.AUC]';
    auc_nov_all = [ROC.nov_all.AUC]';
    
    auc_nov_goalstim = auc_nov_stim(...
        ismember(uniqSess(:,1), novgoalstim(:,1)) & ismember(uniqSess(:,2), novgoalstim(:,2))...
        & ismember(uniqSess(:,3), novgoalstim(:,3)));
    auc_nov_shamstim = auc_nov_stim(...
        ismember(uniqSess(:,1), novshamstim(:,1)) & ismember(uniqSess(:,2), novshamstim(:,2))...
        & ismember(uniqSess(:,3), novshamstim(:,3)));
    novgoalstim(:,4) = auc_nov_goalstim; novgoalstim = sortrows(novgoalstim,2); goalstimLbl = repmat({'Goal'},size(novgoalstim,1),1);
    novshamstim(:,4) = auc_nov_shamstim; shamstimLbl = repmat({'Sham'},size(novshamstim,1),1);
    auc_nov_goalstim = cell2mat(arrayfun( @(x) auc_nov_goalstim(novgoalstim(:,3)==x), 1:3,'UniformOutput',false));
    auc_nov_shamstim = cell2mat(arrayfun( @(x) auc_nov_shamstim(novshamstim(:,3)==x), 1:3,'UniformOutput',false));
    
    figure; ax = arrayfun( @(x) subplot(2,2,x,'NextPlot','add','Box','off','XLim',[0.5,3.5]),1:4);
    set(gcf,'units','inches','position',[5,5,3,4])
    
    iPlot = 1;
    for ii = 4:5
        
        if ii == 1
            auc_nov_stim = [ROC.nov_stim_none.AUC]'; titleLbl = 'NoStim';
            ylabel(ax(ii),'Area under the curve')
        elseif ii == 2
            auc_nov_stim = [ROC.nov_stim_low.AUC]'; titleLbl = 'LowStim';
        elseif ii == 3
            auc_nov_stim = [ROC.nov_stim_high.AUC]'; titleLbl = 'HighStim';
            xlabel(ax(ii),'Day in Novel')
        elseif ii == 4
            auc_nov_stim = [ROC.nov_stim_all.AUC]'; titleLbl = 'BothStim'; % referred as 'BothStim' in Nuri's normalized AUC code
        elseif ii == 5
            auc_nov_stim = [ROC.nov_stim_short.AUC]'; titleLbl = 'ShortStim';
        else
            auc_nov_stim = [ROC.nov_all.AUC]'; titleLbl = 'AllTrials';
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
        auc_nov_shamstim = auc_nov_stim(...
            ismember(uniqSess(:,1), novshamstim(:,1)) & ismember(uniqSess(:,2), novshamstim(:,2))...
            & ismember(uniqSess(:,3), novshamstim(:,3)));
        novgoalstim(:,4) = auc_nov_goalstim; novgoalstim = sortrows(novgoalstim,2); goalstimLbl = repmat({'Goal'},size(novgoalstim,1),1);
        novshamstim(:,4) = auc_nov_shamstim; shamstimLbl = repmat({'Sham'},size(novshamstim,1),1);
        auc_nov_goalstim = cell2mat(arrayfun( @(x) auc_nov_goalstim(novgoalstim(:,3)==x), 1:3,'UniformOutput',false));
        auc_nov_shamstim = cell2mat(arrayfun( @(x) auc_nov_shamstim(novshamstim(:,3)==x), 1:3,'UniformOutput',false));
        
        errorbar(ax(iPlot), nanmean(auc_nov_goalstim),nanstd(auc_nov_goalstim) ./ sqrt(size(auc_nov_goalstim,1)),'Color',params.colors_goalstim(3,:),'LineWidth',2)
        errorbar(ax(iPlot), nanmean(auc_nov_shamstim),nanstd(auc_nov_shamstim) ./ sqrt(size(auc_nov_shamstim,1)),'Color',params.colors_shamstim(3,:),'LineWidth',2)
        xticks(ax(iPlot),1:3)
        title(ax(iPlot), titleLbl)
                

        iPlot = iPlot+1;
        auc_nov_goalstim = (auc_nov_goalstim - auc_nov_goalstim(:,1)) ./ auc_nov_goalstim(:,1) .*100;
        auc_nov_shamstim = (auc_nov_shamstim - auc_nov_shamstim(:,1)) ./ auc_nov_shamstim(:,1) .*100;
        errorbar(ax(iPlot), nanmean(auc_nov_goalstim),nanstd(auc_nov_goalstim) ./ sqrt(size(auc_nov_goalstim,1)),'Color',params.colors_goalstim(3,:),'LineWidth',2)
        errorbar(ax(iPlot), nanmean(auc_nov_shamstim),nanstd(auc_nov_shamstim) ./ sqrt(size(auc_nov_shamstim,1)),'Color',params.colors_shamstim(3,:),'LineWidth',2)
        xticks(ax(iPlot),1:3)
        title(ax(iPlot), [titleLbl, ', change from day1'])
                

        iPlot = iPlot+1;
    end
    suptitle([behType{bt}, ' ROC: goal stim (blue) vs sham stim (orange)'])
    %set(gcf,'units','inches','position',[5,5,8,3])
    makefigurepretty(gcf)
    figname = fullfile(figdir, ['ExtendedDataFigure01_EF_' 'all_Nov_AUC_' behType{bt}]);
    print(gcf, figname, '-dpdf','-r300'); clf;
 
    
    %% wt mice only
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
    
    figure; ax = arrayfun( @(x) subplot(1,2,x,'NextPlot','add','Box','off','XLim',[0.5,3.5]),1:2); 
    errorbar(ax(1), nanmean(auc_nov_wt),nanstd(auc_nov_wt) ./ sqrt(size(auc_nov_wt,1)),'Color',params.colors_nov(3,:),'LineWidth',2)
    errorbar(ax(1), nanmean(auc_fam_wt),nanstd(auc_fam_wt) ./ sqrt(size(auc_fam_wt,1)),'Color',params.colors_fam(3,:),'LineWidth',2)
    xticks(ax(1),1:3)
    xlabel(ax(1),'Day')
    ylabel(ax(1),'Area under the curve')
    title(ax(1), [behType{bt}, ' ROC, familiar (black) vs novel (green)'])

    %plot in % change in AUC from day 1
    auc_fam_wt = (auc_fam_wt - auc_fam_wt(:,1)) ./ auc_fam_wt(:,1) .*100;
    auc_nov_wt = (auc_nov_wt - auc_nov_wt(:,1)) ./ auc_nov_wt(:,1) .*100;     
    errorbar(ax(2), nanmean(auc_nov_wt),nanstd(auc_nov_wt) ./ sqrt(size(auc_nov_wt,1)),'Color',params.colors_nov(3,:),'LineWidth',2)
    errorbar(ax(2), nanmean(auc_fam_wt),nanstd(auc_fam_wt) ./ sqrt(size(auc_fam_wt,1)),'Color',params.colors_fam(3,:),'LineWidth',2)
    xticks(ax(2),1:3)
    xlabel(ax(2),'Day')
    ylabel(ax(2),'Area under the curve')
    title(ax(2), [behType{bt}, ' ROC change from day 1, familiar (black) vs novel (green)'])

    set(gcf,'units','inches','position',[5,5,3,2])
    makefigurepretty(gcf)
    figname = fullfile(figdir, ['ExtendedDataFigure01_EF_' 'all_WT_AUC_' behType{bt}]);
    print(gcf, figname, '-dpdf','-r300'); clf;
    
    
    
end
close all 

