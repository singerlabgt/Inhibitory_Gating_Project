%% Plot SWR power grouped comparison
clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 
baseBins = 1:2; %first 2 bins to average across; to be used as baseline firing rate

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

load(fullfile(maindir, "Demo_Data", "ripplePSTH.mat"))
%%
%time points to be used as x-axis
times = edges ./ ephys_samprate; 

colors = [ 0 0 1; 1 0 0];
%WT Fam: compare between celltypes
figure;
ax = axes('NextPlot','add','Box','off');
yline(ax,0,'k--')
ops.ax = ax;
ops.x_axis = mean(getBinEdges(times),2);
ops.alpha      = 0.5;
ops.line_width = 2;
ops.error      = 'sem';
tempBins = lookup2( [-1, 1], times);%find 1s before and after mid-SWR
for ct = 1:2          
    if ct == 1
        tempR = rate_int_narrow;
        tempSess = sessinfo_int_narrow;
    elseif ct == 2
        tempR = rate_pyr;
        tempSess = sessinfo_pyr;
    end
    tempR = tempR ./ nanmax(tempR,[],2);
    tempSess = table2array(tempSess);
    rate2plot = tempR(ismember(tempSess(:,1), params.WTmice) &...
        tempSess(:,6) == 1 & tempSess(:,8) == 0 & tempSess(:,9) == 0,:);
    rate2plot = (rate2plot - nanmean(rate2plot(:,1:21),2));

    ops.color_area = colors(ct,:);
    ops.color_line = colors(ct,:);
    plot_areaerrorbar(rate2plot, ops); hold on;
end
ax2 = axes('Position',[.6 .6 .3 .3],'Box','off','NextPlot','add');
xline(ax2,0,'k--')
ops.ax = ax2;
ops.x_axis = ops.x_axis(tempBins(1):tempBins(2)-1);
for ct = 1:2          
    if ct == 1
        tempR = rate_int_narrow;
        tempSess = sessinfo_int_narrow;
    elseif ct == 2
        tempR = rate_pyr;
        tempSess = sessinfo_pyr;
    end
    tempR = tempR ./ nanmax(tempR,[],2);
    tempSess = table2array(tempSess);
    rate2plot = tempR(ismember(tempSess(:,1), params.WTmice) &...
        tempSess(:,6) == 1 & tempSess(:,8) == 0 & tempSess(:,9) == 0,:);
    rate2plot = (rate2plot - nanmean(rate2plot(:,1:21),2));

    ops.color_area = colors(ct,:);
    ops.color_line = colors(ct,:);
    plot_areaerrorbar(rate2plot(:,tempBins(1):tempBins(2)-1), ops); hold on;
end

xlabel(ax, 'Time to SWR midpoint (s)')
ylabel(ax, 'Normalized \Delta firing from baseline')
title(ax, 'WT in Familiar')
text(ax, -8, 0.19, 'Pyr','Color','r');
text(ax, -8, 0.16, 'Narrow int','Color','b');
% save figure
figname = 'SuppFigure08_F_left';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%%
%WT Nov: compare between celltypes 
figure;
ax = axes('NextPlot','add','Box','off');
yline(ax,0,'k--')
ops.ax = ax;
ops.x_axis = mean(getBinEdges(times),2);
ops.alpha      = 0.5;
ops.line_width = 2;
ops.error      = 'sem';
for ct = 1:2         
    if ct == 1
        tempR = rate_int_narrow;
        tempSess = sessinfo_int_narrow;
    elseif ct == 2
        tempR = rate_pyr;
        tempSess = sessinfo_pyr;
    end
    tempR = tempR ./ nanmax(tempR,[],2);
    tempSess = table2array(tempSess);
    rate2plot = tempR(ismember(tempSess(:,1), params.WTmice) &...
        tempSess(:,6) ~= 1 & tempSess(:,8) == 0 & tempSess(:,9) == 0,:);
    rate2plot = (rate2plot - nanmean(rate2plot(:,1:21),2));

    ops.color_area = colors(ct,:);
    ops.color_line = colors(ct,:);
    plot_areaerrorbar(rate2plot, ops); hold on;
end
ax2 = axes('Position',[.6 .6 .3 .3],'Box','off','NextPlot','add');
xline(ax2,0,'k--')
ops.ax = ax2;
ops.x_axis = ops.x_axis(tempBins(1):tempBins(2)-1);
for ct = 1:2          
    if ct == 1
        tempR = rate_int_narrow;
        tempSess = sessinfo_int_narrow;
    elseif ct == 2
        tempR = rate_pyr;
        tempSess = sessinfo_pyr;
    end
    tempR = tempR ./ nanmax(tempR,[],2);
    tempSess = table2array(tempSess);
    rate2plot = tempR(ismember(tempSess(:,1), params.WTmice) &...
        tempSess(:,6) ~= 1 & tempSess(:,8) == 0 & tempSess(:,9) == 0,:);
    rate2plot = (rate2plot - nanmean(rate2plot(:,1:21),2));

    ops.color_area = colors(ct,:);
    ops.color_line = colors(ct,:);
    plot_areaerrorbar(rate2plot(:,tempBins(1):tempBins(2)-1), ops); hold on;
end

xlabel(ax, 'Time to SWR midpoint (s)')
ylabel(ax, 'Normalized \Delta firing from baseline')
title(ax, 'WT in Novel')
text(ax, -8, 0.19, 'Pyr','Color','r');
text(ax, -8, 0.16, 'Narrow int','Color','b');
% save figure
figname = 'SuppFigure08_F_right';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% one panel per cell type: Fam vs Nov for WT and GoalStim vs ShamStim for PVxAi32
celltypes = {'Narrow interneuron','Pyramidal cell'};
figure;
ax = arrayfun( @(x) subplot(1,2,x,'NextPlot','add','Box','off'), 1:2);
ops.x_axis = mean(getBinEdges(times),2);
ops.x_axis = ops.x_axis(tempBins(1):tempBins(2)-1);
ops.alpha      = 0.5;
ops.line_width = 2;
ops.error      = 'sem';
for ct = 1:2          
    xline(ax(ct),0,'k--')

    if ct == 1
        tempR = rate_int_narrow;
        tempSess = sessinfo_int_narrow;
    elseif ct == 2
        tempR = rate_pyr;
        tempSess = sessinfo_pyr;
    end
    
    %WT mice 
    ops.ax = ax(ct);
    tempR = tempR ./ max(tempR,[],2);
    tempSess = table2array(tempSess);
    rate2plot = tempR(ismember(tempSess(:,1), params.WTmice) &...
        tempSess(:,6) == 1 & tempSess(:,8) == 0 & tempSess(:,9) == 0,:);
    rate2plot = (rate2plot - nanmean(rate2plot(:,1:21),2));
    ops.color_area = params.colors_fam(3,:);
    ops.color_line = params.colors_fam(3,:);
    plot_areaerrorbar(rate2plot(:,tempBins(1):tempBins(2)-1), ops); hold on;
    
    rate2plot = tempR(ismember(tempSess(:,1), params.WTmice) &...
        tempSess(:,6) ~= 1 & tempSess(:,8) == 0 & tempSess(:,9) == 0,:);
    rate2plot = (rate2plot - nanmean(rate2plot(:,1:21),2));    
    ops.color_area = params.colors_nov(2,:);
    ops.color_line = params.colors_nov(2,:);
    plot_areaerrorbar(rate2plot(:,tempBins(1):tempBins(2)-1), ops); hold on;
    title(ax(ct), celltypes{ct})    
    xlabel(ax(ct), 'Time to SWR midpoint (s)')
    ylabel(ax(ct), 'Normalized \Delta firing from baseline')
    
end
linkaxes(ax, 'xy');
orient(gcf,'landscape')
set(gcf, 'Position', get(0, 'Screensize'));
% save figure
makefigurepretty(gcf)
figname = 'SuppFigure08_G';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

