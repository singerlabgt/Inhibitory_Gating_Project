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

load(fullfile(maindir, "Demo_Data", "allripples_out.mat"))


%% WT mice using RZ rippels only
%plot ripple power 
clear temp v_WT
figure;
isRZripples = abs(out.animalDistance2RZ) <= 10; %within 10 degrees around RZ
v_WT.Fam = out.rip_maxthreshold(isRZripples & ismember(out.sessindex(:,1), params.WTmice) &  out.rip_duringStopped==1 & out.sessindex(:,6)==1);
v_WT.Nov = out.rip_maxthreshold(isRZripples & ismember(out.sessindex(:,1), params.WTmice) &  out.rip_duringStopped==1 & out.sessindex(:,6)~=1);
maxL = max(length(v_WT.Fam), length(v_WT.Nov));
temp = nan(maxL,2);
temp(1:length(v_WT.Fam), 1) = v_WT.Fam;
temp(1:length(v_WT.Nov), 2) = v_WT.Nov;
v1 = violinplot_half(temp);
colormat = [params.colors_fam(3,:); params.colors_nov(2,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 50;
end
text(0.5,9,['ranksum p=' num2str(ranksum(v_WT.Fam,v_WT.Nov))])
ylabel('Ripple power')
xlabel('Fam (black) vs Nov (green)')
title('WT mice')
% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure08_D';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')