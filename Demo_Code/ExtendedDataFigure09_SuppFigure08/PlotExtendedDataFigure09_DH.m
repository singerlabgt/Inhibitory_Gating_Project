%% Plot SWR duration grouped comparison
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

load(fullfile(maindir, "Demo_Data", "allripples_out.mat"))

%% PVxAi32 mice using RZ ripples only
clear temp v_PV
out.rip_duration_s(out.rip_duration_s>100) = nan; %exclude one peculiar one with 500+ sec long ripple 
isRZripples = abs(out.animalDistance2RZ) <= 10; %within 10 degrees around RZ

figure;

v_PV.novGoalStim = out.rip_duration_s(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==1);
v_PV.novShamStim = out.rip_duration_s(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==2);

maxL = max(length(v_PV.novGoalStim), length(v_PV.novShamStim));
temp = nan(maxL,2);
temp(1:length(v_PV.novGoalStim), 1) = v_PV.novGoalStim;
temp(1:length(v_PV.novShamStim), 2) = v_PV.novShamStim;
v1 = violinplot_half(temp); ylim([0,1.05]);
colormat = [params.colors_goalstim(3,:); params.colors_shamstim(3,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 50;
end
text(0.6,9,['ranksum p=' num2str(ranksum(v_PV.novGoalStim,v_PV.novShamStim))])
hold on; breakyaxis([0.25,0.95]); 
ylabel('Ripple duration')
xlabel('Goal stim (blue) vs Sham stim (orange)')
title('PVxAi32 mice using RZ ripples only')
% save figure
figname = 'ExntendedDataFigure09_H';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')


clear temp v_PV
figure;
v_PV.famGoalStim = out.rip_duration_s(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)==1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==1);
v_PV.famNoStim = out.rip_duration_s(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)==1 & out.sessindex(:,8)==0 & out.rip_duringStopped == 1 & out.sessindex(:,9)==0);

maxL = max(length(v_PV.famGoalStim), length(v_PV.famNoStim));
temp = nan(maxL,2);
temp(1:length(v_PV.famGoalStim), 1) = v_PV.famGoalStim;
temp(1:length(v_PV.famNoStim), 2) = v_PV.famNoStim;
v1 = violinplot_half(temp); ylim([0,1.9]);
colormat = [params.colors_goalstim(3,:); params.colors_fam(3,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 50;
end
text(0.6,.18,['ranksum p=' num2str(ranksum(v_PV.famGoalStim,v_PV.famNoStim))])
hold on; breakyaxis([0.4, 1.8])
ylabel('Ripple duration')
xlabel('Fam Goal stim (blue) vs No stim (black)')
title('PVxAi32 mice using all ripples')
% save figure
figname = 'ExtendedDataFigure09_DH';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')