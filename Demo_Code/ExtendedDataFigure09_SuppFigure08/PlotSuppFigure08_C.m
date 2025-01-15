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
%% WT mice only 
%plot coactivity probability
clear temp
figure;
v_WT.Fam = table2array(coactivePairs_fam(ismember(coactivePairs_fam.Animal, params.WTmice) & coactivePairs_fam.VR == 1, end));
v_WT.Nov = table2array(coactivePairs_nov(ismember(coactivePairs_nov.Animal, params.WTmice) & coactivePairs_nov.VR ~= 1, end));
maxL = max(height(v_WT.Fam), height(v_WT.Fam));
temp = nan(maxL,2);
temp(1:height(v_WT.Fam), 1) = v_WT.Fam;
temp(1:height(v_WT.Nov), 2) = v_WT.Nov;
v1 = violinplot_half(temp);
colormat = [params.colors_fam(3,:); params.colors_nov(2,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 50;
end
text(0.5,0.5,['ranksum p=' num2str(ranksum(v_WT.Fam,v_WT.Nov))])
ylabel('Coactivation probability during SWR')
xlabel('Fam (black) vs Nov (green)')
title('WT mice')
% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure08_C';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
