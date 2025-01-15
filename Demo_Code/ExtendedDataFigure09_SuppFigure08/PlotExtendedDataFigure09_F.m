%% Plot Co-activation Proability during SWR

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

load(fullfile(maindir, "Demo_Data/", "coactivationProb.mat"))

%% PVxAi32 mice only 
clear temp 
figure;
v_PV.novGoalStim = table2array(coactivePairs_nov(...
    ismember(coactivePairs_nov.Animal, params.goalshamMice)...
    & coactivePairs_nov.VR ~= 1 & coactivePairs_nov.Stimulation == 1 ... 
    & coactivePairs_nov.StimLocation == 1, end));
v_PV.novShamStim = table2array(coactivePairs_nov(...
    ismember(coactivePairs_nov.Animal, params.goalshamMice)...
    & coactivePairs_nov.VR ~= 1 & coactivePairs_nov.Stimulation == 1 ...
    & coactivePairs_nov.StimLocation == 2, end));

maxL = max(height(v_PV.novGoalStim), height(v_PV.novShamStim));
temp = nan(maxL,2);
temp(1:height(v_PV.novGoalStim), 1) = v_PV.novGoalStim;
temp(1:height(v_PV.novShamStim), 2) = v_PV.novShamStim;
v1 = violinplot_half(temp);
colormat = [params.colors_goalstim(3,:); params.colors_shamstim(3,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 50;
end
text(0.6,0.5,['ranksum p=' num2str(ranksum(v_PV.novGoalStim,v_PV.novShamStim))])
ylabel('Coactivation probability during SWR')
xlabel('Goal stim (blue) vs Sham stim (orange)')
title('PVxAi32 mice')
% save figure
makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure09_F';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

