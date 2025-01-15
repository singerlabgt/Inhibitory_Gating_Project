%% plot residual firing rate hit vs miss trials per celltype
clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

sec_before = 1; % as duration
sec_after = 1;
binSize= 0.05; % 50ms bin
edges = -1* sec_before : binSize : sec_after;
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end
fileID = fopen(fullfile(figdir, 'ExtendedDataFigure06_stats.txt'), 'w'); fclose(fileID); %create stats info file 
%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);

bins2incl = 11:13; %find peak or trough 10 degrees around RZ (-10 to 20) 
% bins2incl = 11:16; %old version from Nuri
ctNames = {'NS Int','Pyr'}; %excl WS 

%% load residual distance to RZ data file
load(fullfile(maindir, "Demo_Data/", "allsess_raw_vs_residuals_distance2RZ.mat"));
celltype = zeros(length(allsess_unitID),1);
celltype( arrayfun( @(x) strcmp(allsess_unitType{x},'Narrow Interneuron'), 1:length(allsess_unitType))) = 1;
celltype( arrayfun( @(x) strcmp(allsess_unitType{x},'Wide Interneuron'), 1:length(allsess_unitType))) = 3;
celltype( arrayfun( @(x) strcmp(allsess_unitType{x},'Pyramidal Cell'), 1:length(allsess_unitType))) = 2;

%% plot mean +/- SEM residual change from baseline hit vs miss trials in WT mice
fig = figure('units','inch','position',[0,0,5,4]);
ax = arrayfun( @(x) subplot(4,2,x,'NextPlot','add','Box','off'), 1:8);
arrayfun( @(x) yline(ax(x),0,'k:'),1:length(ax));
arrayfun( @(x) xline(ax(x),0,'k:'),1:length(ax));

normMap_famHit = (allsess_mean_fam_residual_hit - nanmin(allsess_mean_fam_residual_hit,[],2))...
    ./ (nanmax(allsess_mean_fam_residual_hit,[],2) - nanmin(allsess_mean_fam_residual_hit,[],2));
normMap_famMiss = (allsess_mean_fam_residual_miss - nanmin(allsess_mean_fam_residual_miss,[],2))...
    ./ (nanmax(allsess_mean_fam_residual_miss,[],2) - nanmin(allsess_mean_fam_residual_miss,[],2));
normMap_novHit = (allsess_mean_nov_residual_hit - nanmin(allsess_mean_nov_residual_hit,[],2))...
    ./ (nanmax(allsess_mean_nov_residual_hit,[],2) - nanmin(allsess_mean_nov_residual_hit,[],2));
normMap_novMiss = (allsess_mean_nov_residual_miss - nanmin(allsess_mean_nov_residual_miss,[],2))...
    ./ (nanmax(allsess_mean_nov_residual_miss,[],2) - nanmin(allsess_mean_nov_residual_miss,[],2));
    
data4R_fam = []; data4R_nov = []; %data structure long form saved for LMM in R
for ct = 1:length(ctNames)
   
    i2plot = ( celltype == ct & ismember(allsess_sessinfo(:,1), params.WTmice));
    
    ops.ax = ax(ct);
    ops.x_axis = mean(getBinEdges(position_binEdges),2);
    ops.color_area = [0 0 0];
    ops.color_line = [0 0 0];
    ops.alpha = 0.2;
    ops.line_width = 2;
    ops.error = 'sem';
    
    %% fam hit vs miss
    axes(ax(ct))
    hold on
    map2plotH = normMap_famHit(i2plot,:);
    map2plotH = map2plotH - nanmean(map2plotH(:,1:2),2);    
    sessinfoH = [allsess_sessinfo(i2plot,:), allsess_unitID(i2plot), celltype(i2plot)];
    sessinfoH = sessinfoH(all(~isnan(map2plotH),2),:);
    map2plotH = map2plotH(all(~isnan(map2plotH),2),:); %remove NaN cells
    plot_areaerrorbar(map2plotH, ops); hold on
    
    clear trialLbl 
    ops.color_area = [0.5 0 0];
    ops.color_line = [0.5 0 0];
    map2plotF = normMap_famMiss(i2plot,:);
    map2plotF = map2plotF - nanmean(map2plotF(:,1:2),2);
    sessinfoF = [allsess_sessinfo(i2plot,:), allsess_unitID(i2plot), celltype(i2plot)];
    sessinfoF = sessinfoF(all(~isnan(map2plotF),2),:);
    map2plotF = map2plotF(all(~isnan(map2plotF),2),:); %remove NaN cells
    plot_areaerrorbar(map2plotF, ops)
    title(ax(ct), ['Fam ' ctNames{ct}])
    
    %quantificaiton correct vs incorrect
    if ct == 1 %find trough within specified bin positions for NS units
        plotHits = nanmin(map2plotH(:,bins2incl),[],2);
        plotFail = nanmin(map2plotF(:,bins2incl),[],2);
    else %find peak if Pyr cell 
        plotHits = nanmax(map2plotH(:,bins2incl),[],2);
        plotFail = nanmax(map2plotF(:,bins2incl),[],2);
    end
    scatter(ax(ct+2), ones(size(map2plotH,1),1) .* 1, plotHits,...
        'MarkerFaceColor', [0 0 0], 'MarkerFaceAlpha', 0.1,...
        'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.1,...
        'jitter', 'on', 'jitterAmount', 0.1);
    scatter(ax(ct+2), ones(size(map2plotF,1),1) .*2, plotFail,...
        'MarkerFaceColor', [0.5 0 0], 'MarkerFaceAlpha', 0.1,...
        'MarkerEdgeColor', [0.5 0 0], 'MarkerEdgeAlpha', 0.1,...
        'jitter', 'on', 'jitterAmount', 0.1);
    errorbar(ax(ct+2), [1.3, 1.7], [nanmean(plotHits), nanmean(plotFail)],...
        [nanstd(plotHits), nanstd(plotFail)],...
        'Color',[0 0 0],'Marker','o','MarkerEdgeColor',[0 0 0],...
        'LineWidth',2,'MarkerFaceColor','w');
    title(ax(ct+2),['p=' num2str(ranksum(plotHits, plotFail))])
    xlim(ax(ct+2), [0.5,2.5]); xticks(ax(ct+2), [1, 2]); xticklabels(ax(ct+2), []);
    
    %stats info
    hits = nanmean(plotHits,2); trialLbl = ones(length(hits),1);
    fails = nanmean(plotFail,2); trialLbl = [trialLbl; zeros(length(fails),1)];
    
    data4R_fam = [data4R_fam;  [hits; fails], trialLbl, [sessinfoH; sessinfoF] ]; %[residualFR, correct or not, sessInfo, cellID, celltype]
    fileID = fopen(fullfile(figdir, 'stats.txt'), 'a');
    fprintf(fileID,'%s\r\n', deblank(['WT Fam Correct ' ctNames{ct} ' [min, 25th, 50th, 75th, max] = ' num2str( prctile(hits, [0,25,50,75,100])) ]));
    fprintf(fileID,'%s\r\n', deblank(['WT Fam Correct ' ctNames{ct} ' n = ' num2str(length(hits))...
        ', mean +/- SEM = ' num2str(mean(hits)) ' +/- ' num2str(std(hits) ./ sqrt(length(hits)))]));
    fprintf(fileID,'%s\r\n', deblank(['WT Fam Incorrect ' ctNames{ct} ' [min, 25th, 50th, 75th, max] = ' num2str( prctile(fails, [0,25,50,75,100])) ]));
    fprintf(fileID,'%s\r\n', deblank(['WT Fam Incorrect ' ctNames{ct} ' n = ' num2str(length(fails))...
        ', mean +/- SEM = ' num2str(mean(fails)) ' +/- ' num2str(std(fails) ./ sqrt(length(fails)))]));

    
    %% nov hit vs miss
    clear trialLbl sessinfoF sessinfoH
    ops.ax = ax(ct+4);
    axes(ax(ct+4))
    hold on
    ops.color_area = [0 0 0];
    ops.color_line = [0 0 0];
    map2plotH = normMap_novHit(i2plot,:);
    map2plotH = map2plotH - mean(map2plotH(:,1:2),2);
    sessinfoH = [allsess_sessinfo(i2plot,:), allsess_unitID(i2plot), celltype(i2plot)];
    sessinfoH = sessinfoH(all(~isnan(map2plotH),2),:);
    map2plotH = map2plotH(all(~isnan(map2plotH),2),:); %remove NaN cells
    
    plot_areaerrorbar(map2plotH, ops); hold on
    ops.color_area = [0.5 0 0];
    ops.color_line = [0.5 0 0];
    map2plotF = normMap_novMiss(i2plot,:);
    map2plotF = map2plotF - mean(map2plotF(:,1:2),2);
    sessinfoF = [allsess_sessinfo(i2plot,:), allsess_unitID(i2plot), celltype(i2plot)];
    sessinfoF = sessinfoF(all(~isnan(map2plotF),2),:);
    map2plotF = map2plotF(all(~isnan(map2plotF),2),:); %remove NaN cells
    plot_areaerrorbar(map2plotF, ops)
    title(ax(ct+4), ['Nov ' ctNames{ct}])
    
    
    if ct == 1 %find trough within specified bin positions for NS units
        plotHits = nanmin(map2plotH(:,bins2incl),[],2);
        plotFail = nanmin(map2plotF(:,bins2incl),[],2);
    else %find peak if Pyr cell 
        plotHits = nanmax(map2plotH(:,bins2incl),[],2);
        plotFail = nanmax(map2plotF(:,bins2incl),[],2);
    end
    scatter(ax(ct+6), ones(size(map2plotH,1),1) .* 1, plotHits,...
        'MarkerFaceColor', [0 0 0], 'MarkerFaceAlpha', 0.1,...
        'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.1,...
        'jitter', 'on', 'jitterAmount', 0.1); hold on; 
    scatter(ax(ct+6), ones(size(map2plotF,1),1) .*2, plotFail,...
        'MarkerFaceColor', [0.5 0 0], 'MarkerFaceAlpha', 0.1,...
        'MarkerEdgeColor', [0.5 0 0], 'MarkerEdgeAlpha', 0.1,...
        'jitter', 'on', 'jitterAmount', 0.1);
    errorbar(ax(ct+6), [1.3, 1.7], [nanmean(plotHits), nanmean(plotFail)],...
        [nanstd(plotHits), nanstd(plotFail)],...
        'Color',[0 0 0],'Marker','o','MarkerEdgeColor',[0 0 0],...
        'LineWidth',2,'MarkerFaceColor','w');    
    xlim(ax(ct+6), [0.5,2.5]); xticks(ax(ct+6), [1, 2]);  xticklabels(ax(ct+6), []);
    title(ax(ct+6),['p=' num2str(ranksum(plotHits, nanmin(plotFail)))])
    
    %stats info
    hits = nanmean(plotHits,2); trialLbl = ones(length(hits),1);
    fails = nanmean(plotFail,2); trialLbl = [trialLbl; zeros(length(fails),1)];
    
    data4R_nov = [data4R_nov;  [hits; fails], trialLbl, [sessinfoH; sessinfoF] ]; %[residualFR, correct or not, sessInfo, cellID, celltype]
    
    fileID = fopen(fullfile(figdir, 'stats.txt'), 'a');
    fprintf(fileID,'%s\r\n', deblank(['WT Nov Correct ' ctNames{ct} ' [min, 25th, 50th, 75th, max] = ' num2str( prctile(hits, [0,25,50,75,100])) ]));
    fprintf(fileID,'%s\r\n', deblank(['WT Nov Correct ' ctNames{ct} ' n = ' num2str(length(hits))...
        ', mean +/- SEM = ' num2str(mean(hits)) ' +/- ' num2str(std(hits) ./ sqrt(length(hits)))]));
    fprintf(fileID,'%s\r\n', deblank(['WT Nov Incorrect ' ctNames{ct} ' [min, 25th, 50th, 75th, max] = ' num2str( prctile(fails, [0,25,50,75,100])) ]));
    fprintf(fileID,'%s\r\n', deblank(['WT Nov Incorrect ' ctNames{ct} ' n = ' num2str(length(fails))...
        ', mean +/- SEM = ' num2str(mean(fails)) ' +/- ' num2str(std(fails) ./ sqrt(length(fails)))]));


end
xlim(ax([1:2, 5:6]), [-40 20])
linkaxes(ax([1:2, 5:6]),'xy')
linkaxes(ax([3:4, 7:8]),'xy')
xlabel(ax(7),'Distance to RZ (deg)')
ylabel(ax(7),'Norm. \Delta residual FR from baseline')
for a = 1:8
    ax(a);
    set(ax(a), 'TickDir', 'out')
end
makefigurepretty(gcf);
figname = 'ExtendedDataFigure06_AB';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')