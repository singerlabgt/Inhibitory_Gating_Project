clc; clear; close all;
%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';

addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end
%get all indices for baseline and behavioral files
[dirs, params] = getDefaultParameters(maindir); 
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.animals);
allindex = allindex(animal_idx,:); %filter based on animals to include 
allindex = allindex(animal_idx,:); %filter based on animals to include 
[sessions, ~] = unique(allindex(:,[1:2,7]), 'rows'); %define session as one date

load(fullfile(maindir, 'Demo_Data', 'firingrate_perSpeedQuartile'))
%% ED Fig. 3A1: use all narrow interneurons, FAM
adjusted_p = arrayfun( @(x)...
    signrank(all_rate.fam(all_info.fam(:,5)==1,x),all_rate.noVR(all_info.fam(:,5)==1,x)),...
    1:size(all_rate.fam,2)) .* size(all_rate.fam,2);
% plot one line per NS interneuron, compare FR between Fam and noVR per speed quartile
figure('Units', 'inches', 'Position', [0 0 2 2]); ax = axes('NextPlot','add','Box','off','XLim',[0.5 7],'XTick',1.5:1.5:7,'XTickLabel',{'1','2','3','4'}); 
prev = 1;
for sQ = 1:4
    temp = [all_rate.fam(all_info.fam(:,5)==1,sQ), all_rate.noVR(all_info.fam(:,5)==1,sQ)]; 
    temp = temp(all(~isnan(temp),2),:); 
    plot(ax, (ones(size(temp,1),2) .* [prev, prev+1])', temp','Color',params.colors_narrowInt(1,:));
    errorbar(ax, (ones(1,2) .* [prev, prev+1])', mean(temp), std(temp) ./ sqrt(size(temp,1)),'LineWidth',2,'Color',params.colors_narrowInt(3,:));
    text(ax, prev-0.5, 90,['p=' num2str(signrank(temp(:,1), temp(:,2)))])
    disp(['Fam NS Interneurons - SQ ', num2str(sQ), ' - P = ', num2str(signrank(temp(:,1), temp(:,2)))]) 
    prev = prev +1.5;    
end
title('Fam')
xlabel('Speed quartile'); ylabel('Firing rate (spikes/s)');
makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure03_A1';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% ED Fig. 3A2: use all narrow interneurons, NOV
adjusted_p = arrayfun( @(x)...
    signrank(all_rate.nov(all_info.nov(:,5)==1,x),all_rate.noVR(all_info.nov(:,5)==1,x)),...
    1:size(all_rate.nov,2)) .* size(all_rate.nov,2);
% plot one line per NS interneuron, compare FR between nov and noVR per speed quartile
figure('Units', 'inches', 'Position', [0 0 2 2]); ax = axes('NextPlot','add','Box','off','XLim',[0.5 7],'XTick',1.5:1.5:7,'XTickLabel',{'1','2','3','4'}); 
prev = 1;
for sQ = 1:4
    temp = [all_rate.nov(all_info.nov(:,5)==1,sQ), all_rate.noVR(all_info.nov(:,5)==1,sQ)]; 
    temp = temp(all(~isnan(temp),2),:); 
    plot(ax, (ones(size(temp,1),2) .* [prev, prev+1])', temp','Color',params.colors_narrowInt(1,:));
    errorbar(ax, (ones(1,2) .* [prev, prev+1])', mean(temp), std(temp) ./ sqrt(size(temp,1)),'LineWidth',2,'Color',params.colors_narrowInt(3,:));
    text(ax, prev-0.5, 90,['p=' num2str(signrank(temp(:,1), temp(:,2)))])
        disp(['Nov NS Interneurons - SQ ', num2str(sQ), ' - P = ', num2str(signrank(temp(:,1), temp(:,2)))]) 
    prev = prev +1.5;    
end
title('Nov')
xlabel('Speed quartile'); ylabel('Firing rate (spikes/s)');
makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure03_A2';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% ED Fig. 3A3: abs speed values at each quartile
% unsure which color to use...select FAM color temporarily.
figure('Units', 'inches', 'Position', [0 0 2 1.8]); ax = axes('NextPlot','add','Box','off','XLim',[0.5 4.5],'XTick',1:4,'XTickLabel',{'1','2','3','4'}); 
plot(speed_alldaymat', 'Color', params.colors_fam(1,:), 'LineWidth',1); hold on
errorbar(1:4, mean(speed_alldaymat), std(speed_alldaymat), 'Color', params.colors_fam(end,:), 'LineWidth',2); hold off
xlabel('Speed quartile'); ylabel('Speed (deg/s)');
title('Fig S4A, absolute speed distribution')
makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure03_A3';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% ED Fig. 3B1
celltypes = {'PV Interneuron','Narrow Interneuron','Wide Interneuron','Pyramidal Cell'};
figure;
ax = arrayfun( @(x) subplot(1,3,x,'NextPlot','add','Box','off','XLim',[0.5 4.5],'xtick',1:4),1:3);
for iCT = 1:3
    tempN = rateQ(rateQ(:,end)==iCT,:);
    if iCT == 1
        colors = params.colors_narrowInt;
    elseif iCT == 2
        colors = params.colors_wideInt;
    elseif iCT == 3
        colors = params.colors_pyr;
    end
    
    plot(ax(iCT), tempN(:,1:4)','color',colors(1,:))
    errorbar(ax(iCT), 1:4, arrayfun( @(x) mean(tempN(:,x)),1:4),...
        arrayfun( @(x) std(tempN(:,x)) ./ sqrt(length(tempN(:,x))),1:4),'Color',colors(3,:),'LineWidth',2)
    title(ax(iCT),[celltypes{iCT+1} ' (1 line = 1 unit)']);
    ylabel(ax(iCT),'AZ trial-averaged firing rate (spikes/s)');
    xlabel(ax(iCT),'Speed quartile')
end
axis(ax,'square');
makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure03_B1';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% EDFig. 3B2: calculate by RZ-centered correct trials.
prct2use = 0:25:100; %equally sized speed percentiles
speed_alldaymat = nan([length(sessions), length(prct2use) - 1]);
for iSess = 1:length(sessions)
    filedir = fullfile(maindir, "Demo_Data/", "singlesess_ratemap_distance2RZ/");
    filename = getlatestfile_with_string(filedir, [params.iden num2str(sessions(iSess,1)) '_' num2str(sessions(iSess,2))]);
    load(fullfile(filedir, filename));
    temp = mean(outmap.bin5.fam.smoothSpeed(outmap.bin5.fam.labels(:,2)==1 & outmap.bin5.fam.labels(:,3)==1,11:12),2);
    speedQ = prctile(temp,[0 25 50 75 100]);
    abs_speedDist_avg = arrayfun(@(x) mean(temp(temp >= speedQ(x) & temp <= speedQ(x + 1))), 1 : length(speedQ) - 1); % average speed within that quartile
    speed_alldaymat(iSess, :) = abs_speedDist_avg; 
end
figure('Units', 'inches', 'Position', [0 0 2 1.8]); ax = axes('NextPlot','add','Box','off','XLim',[0.5 4.5],'XTick',1:4,'XTickLabel',{'1','2','3','4'}); 
plot(speed_alldaymat', 'Color', params.colors_fam(1,:), 'LineWidth',1); hold on
errorbar(1:4, mean(speed_alldaymat), std(speed_alldaymat), 'Color', params.colors_fam(end,:), 'LineWidth',2); hold off
xlabel('Speed quartile'); ylabel('Speed (deg/s)');
title('Fig S4B, absolute speed distribution')
makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure03_B2';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
