%% plot distribution of ripple rate in Hz per animal/day/environment

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

load(fullfile(maindir, "Demo_Data", "riprate_for_R.mat"))
load(fullfile(maindir, "Demo_Data", "allindex_ripplerate.mat"))
load(fullfile(maindir, "Demo_Data", "RipData_250ms.mat"))

riprate = allindex.ripsInStopped./ allindex.StoppedDuration;
numTh = 10;
%% ExtendedDataFigure08_B: all three days combined, compare fam vs nov ripple rate
figure;
clear v_WT
v_WT.fam = riprate(ismember(allindex.Animal,params.WTmice)...
    & ismember(allindex.Date, sessions(ripples_per_session >= numTh,2))...
    & allindex.VR == 1 );
v_WT.nov = riprate(ismember(allindex.Animal,params.WTmice)...
    & ismember(allindex.Date, sessions(ripples_per_session >= numTh,2))...
    & allindex.VR ~= 1 );
maxL = max(length(v_WT.fam), length(v_WT.nov));
temp = nan(maxL,2);
temp(1:length(v_WT.fam), 1) = v_WT.fam;
temp(1:length(v_WT.nov), 2) = v_WT.nov;
v1 = violinplot_half(temp);
colormat = [params.colors_fam(end,:); params.colors_nov(2,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 100;
end
xlabel('Fam (black) vs Nov (green)')
ylabel('Ripple rate (Hz)')
title(['WT mice with at least ' num2str(numTh) ' ripples']);
% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure08_B';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
