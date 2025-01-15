%% plot_behaviorHistROC
%plot behavioral performance metric distribution comparing AZ/RZ vs NRZ along with ROC curves across days
%plots speed/lickrate/licklatency histograms across days along with ROC curves
%NJ created 12/31/2022
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
savefigdir = fullfile(figdir, 'ExtendedDatatFigure01_BCD_BehaviorHistROC'); if ~isfolder(savefigdir); mkdir(savefigdir); end
load(fullfile(maindir, 'Demo_Data', 'getNovelBehaviorROC_Speed_XZ_250105.mat'))
structdir = fullfile(maindir, 'Demo_Data', 'singlesess_behavior_speeddist/');
%% EDfig. 1BCD example speed distributions and ROC
%get all indices
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.animals);
allindex = allindex(animal_idx,:); %filter based on animals to include 
allindex = allindex(animal_idx,:); %filter based on animals to include 
[sets, ~, ~] = splitSessions2Set(allindex);

behType = {'Speed'}; bt = 1;
xlbls = {'Speed (deg/s)'};
currBT = [behType{bt} 'ROC'];
roc_sessinfo = cell2mat( arrayfun( @(x)...
    unique(ROC.sessInfo{x}(:,[1:2,7]),'rows'),...
    1:length(ROC.sessInfo),'UniformOutput',false)' );

for ss = 1:length(sets)-2
    animal = unique( sets{ss}(:, 1) );
    recDates = sets{ss}(:, 2);
    if animal == 4
        iden = 'X';
    else
        iden = 'N';
    end
    ax = arrayfun( @(x) subplot(3,4,x,'NextPlot','add','Box','off'), 1:12);
    for d = 1:length(recDates)
        recDay = recDates(d);
        novelDay = sets{ss}(d, 3);
        idx = find(roc_sessinfo(:,1) == animal & roc_sessinfo(:,2) == recDay & roc_sessinfo(:,3) == novelDay);
        stringname = [iden num2str(animal) '_' num2str(recDay)];
        
        %get latest file with name that contains the specified string
        fname = getlatestfile_with_string(structdir, stringname);
        load(fullfile(structdir, fname));
        
        temp = [nanmean(data.noVR.az,2); nanmean(data.noVR.cz,2); nanmean(data.fam.az,2); nanmean(data.fam.cz,2); nanmean(data.nov.az,2); nanmean(data.nov.cz,2)];
        edges = min(temp) : 0.5 : prctile(temp,95); 
        
        %NoVR
        histogram(ax(d), nanmean(data.noVR.az,2),'EdgeColor',params.colors_fam(3,:),'BinEdges',edges,'Normalization','probability','DisplayName','NoVR AZ','DisplayStyle','stairs','LineWidth',2);
        histogram(ax(d), nanmean(data.noVR.cz,2),'EdgeColor',[0.7 0.5 0.7],'BinEdges',edges,'Normalization','probability','DisplayName','NoVR NRZ','DisplayStyle','stairs','LineWidth',2);
        axis square
        plot(ax(4), ROC.noVR_all.X{idx}, ROC.noVR_all.Y{idx}, 'Color', params.colors_fam(d,:), 'LineWidth',2, 'DisplayName',['Day ' num2str(d)])
        axis square
        
        %Fam
        histogram(ax(d+4), nanmean(data.fam.az,2),'EdgeColor',params.colors_fam(3,:),'BinEdges',edges,'Normalization','probability','DisplayName','Fam AZ','DisplayStyle','stairs','LineWidth',2);
        histogram(ax(d+4), nanmean(data.fam.cz,2),'EdgeColor',[0.7 0.5 0.7],'BinEdges',edges,'Normalization','probability','DisplayName','Fam NRZ','DisplayStyle','stairs','LineWidth',2);
        axis square
        plot(ax(8), ROC.fam_all.X{idx}, ROC.fam_all.Y{idx}, 'Color', params.colors_fam(d,:), 'LineWidth',2, 'DisplayName',['Day ' num2str(d)])
        axis square
        
        %Nov
        histogram(ax(d+8), nanmean(data.nov.az,2),'EdgeColor',params.colors_nov(3,:),'BinEdges',edges,'Normalization','probability','DisplayName','Nov AZ','DisplayStyle','stairs','LineWidth',2);
        histogram(ax(d+8), nanmean(data.nov.cz,2),'EdgeColor',[0.7 0.5 0.7],'BinEdges',edges,'Normalization','probability','DisplayName','Nov NRZ','DisplayStyle','stairs','LineWidth',2);
        axis square            
        plot(ax(12), ROC.nov_all.X{idx}, ROC.nov_all.Y{idx}, 'Color', params.colors_nov(d,:), 'LineWidth',2, 'DisplayName',['Day ' num2str(d)])
        axis square
        xlabel(ax(10), xlbls(bt))
        ylabel(ax(5), 'Proportion of trials')
        xlabel(ax(12), 'False positive rate')
        ylabel(ax(12), 'True positive rate')
        clear data
    end
    axis(ax, 'square')
    linkaxes(ax(1:3), 'xy')
    linkaxes(ax(5:7), 'xy')
    linkaxes(ax(9:11), 'xy')
    linkaxes(ax([4,8,12]), 'xy')
    arrayfun( @(x) title(ax(x), ['Day ' num2str(x)]), 1:3)
    arrayfun( @(x) plot(ax(x), [0 1], [0 1], 'k--'), [4,8,12])
    suptitle([iden num2str(animal) ' - ' currBT ' noVR/Fam/Nov'])
    set(gcf, 'Position', get(0, 'Screensize'));
    orient(gcf,'landscape')
    currentsavefigdir = fullfile(savefigdir, ['Set', num2str(ss) '_' iden num2str(animal) ' - ' currBT]);
    saveas(gcf, currentsavefigdir, 'png')
    print(gcf, currentsavefigdir,'-dpdf','-r300','-bestfit')
    clf; 
end