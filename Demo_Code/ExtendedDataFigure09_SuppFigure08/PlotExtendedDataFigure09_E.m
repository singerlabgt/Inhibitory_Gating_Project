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

numTh = 10;

riprate = allindex.ripsInStopped./ allindex.StoppedDuration;
%% SuppFigure08_E: PV animals with goal and sham stim conditions
%compare overall difference between novGoalStim vs novShamStim
figure;
clear v_PV
v_PV.novGoalStim = riprate(ismember(allindex.Animal,params.goalshamMice)...
    & ismember(allindex.Date, sessions(ripples_per_session >= numTh,2))...
    & allindex.VR ~= 1 & allindex.Stimulation == 1 & allindex.StimLocation == 1 );
v_PV.novShamStim = riprate(ismember(allindex.Animal,params.goalshamMice)...
    & ismember(allindex.Date, sessions(ripples_per_session >= numTh,2))...
    & allindex.VR ~= 1 & allindex.Stimulation == 1 & allindex.StimLocation == 2 );
maxL = max(length(v_PV.novGoalStim), length(v_PV.novShamStim));
temp = nan(maxL,2);
temp(1:length(v_PV.novGoalStim), 1) = v_PV.novGoalStim;
temp(1:length(v_PV.novShamStim), 2) = v_PV.novShamStim;
v1 = violinplot_half(temp);
colormat = [params.colors_goalstim(2,:); params.colors_shamstim(2,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 100;
end
ylabel('Ripple rate (Hz)')
title(['PVxAi32 mice with at least ' num2str(numTh) ' ripples']);
% save figure
makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure09_E';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
