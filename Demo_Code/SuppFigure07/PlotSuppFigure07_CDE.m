%% circular stats analysis
clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%%
%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

load(fullfile(maindir, 'Demo_Data', 'thetaphase_stim'))
load(fullfile(maindir, 'Demo_Data', 'thetapower_stim'))
load(fullfile(maindir, 'Demo_Data', 'thetatrace_example.mat'))

%% Supp Fig. 7C
% plot goal stim theta trace example
figure; t = tiledlayout(2, 1);
nexttile;
x_time = x_time_goalstim;
plot(x_time, filtered_theta_goalstim)
xlabel('Time (Sec)'); xlim([0 1]); ylim([-500 500])
ylabel('Amplitude (uV)')
title('filtered theta trace')
nexttile;
plot(x_time, raw_lfptrace_goalstim)
vel = round(mean(velocity_goalstim), 2);
xlabel('Time'); xlim([0 1]); ylim([-500 500])
ylabel('Amplitude (uV)')
title('raw LFP trace')
sgtitle(['Goal stim example, vel=' num2str(vel) ' deg/s'])
makefigurepretty(gcf);
figname = ['SuppFigure07_C_thetatrace_goalstim'];
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

% plot sham stim theta trace example
figure; t = tiledlayout(2, 1);
nexttile;
x_time = x_time_shamstim;
plot(x_time, filtered_theta_shamstim)
xlabel('Time (Sec)'); xlim([0 1]); ylim([-500 500])
ylabel('Amplitude (uV)')
title('filtered theta trace')
nexttile;
plot(x_time, raw_lfptrace_shamstim)
vel = round(mean(velocity_shamstim), 2);
xlabel('Time'); xlim([0 1]); ylim([-500 500])
ylabel('Amplitude (uV)')
title('raw LFP trace')
sgtitle(['sham stim example, vel=' num2str(vel) ' deg/s'])
makefigurepretty(gcf);
figname = ['SuppFigure07_C_thetatrace_shamstim'];
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

% plot theta power distribution. Theta is extracted during per running
% block when the animal runs at speed greater than 2 degrees/sec
figure('Units', 'inches', 'Position', [5.4167 3.1771 2.0521 1.5833]);
maxL = max(length(thetapower_per_running_block{1}), length(thetapower_per_running_block{2}));
temp = nan(maxL,2);
temp(1:length(thetapower_per_running_block{1}), 1) = thetapower_per_running_block{1};
temp(1:length(thetapower_per_running_block{2}), 2) = thetapower_per_running_block{2};
v1 = violinplot_half(temp, [], 'CenterSpace', 0.05, 'MedianSize', 20);
% v1 = violinplot_half(temp, 'MedianSize', 50);
colormat = [params.colors_goalstim(2,:); params.colors_shamstim(2,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 15;
end
ylabel(['Theta power (uV^2)']) 
xticks([])
ylim([0 3.5*1e5])
makefigurepretty(gcf);
figname = ['SuppFigure07_C_thetapower'];
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
%% Supp Fig. 7D
% Calculate vector strength
vec_strength_goalstim = []; vec_strength_shamstim = [];
figure;
t = tiledlayout(2, 2);
for condition = 1:2
    if condition == 1
        angles = thetaphase_goalstim;
        preferred = deg2rad(preferredthetaphase_goalstim);
        colors = params.colors_goalstim;
        cond = 'goalstim';
    else
        angles = thetaphase_shamstim;
        preferred = deg2rad(preferredthetaphase_shamstim);
        colors = params.colors_shamstim;
        cond = 'shamstim';
    end

    % calculate vector strength for each cell, then plot histgram
    vec_strength = nan(length(angles), 1);
    for iCell = 1:length(angles)
        R = circ_r(angles{iCell});
        vec_strength(iCell) = R;
    end
    nexttile;
    %histogram(vec_strength, 'Normalization', 'probability', 'FaceColor', colors(end,:));
    vs_edges = 0:0.1:1;
    counts = histcounts(vec_strength, vs_edges);
    nh = counts./sum(counts);
    plot(vs_edges(1:end-1), nh, 'Color', colors(end,:), 'LineWidth', 1)
    ylabel('Fraction of cells')
    xlabel('Vector strength')
    title('vector strength')
    xticks([0 0.5 1])
    ylim([0 0.25])
    yticks([0 0.1 0.2])

    if cond == 'goalstim'
        vec_strength_goalstim = [vec_strength_goalstim; vec_strength];
    else
        vec_strength_shamstim = [vec_strength_shamstim; vec_strength];
    end
    % bins = 0:0.02:0.14;
    % h = histcounts(vec_strength, bins);
    % nh = h./sum(h)
    % plot(bins(1:end-1),nh)

    % Compute histogram data
    num_bins = 36;
    [counts, bin_centers] = histcounts(preferred, num_bins);
    
    % Convert bin centers to polar coordinates (radians)
    bin_centers = bin_centers(1:end-1) + diff(bin_centers)/2;
    
    % Plot the polar histogram
    nexttile;
    polarhistogram(counts, num_bins)
    polarhistogram('BinEdges', [bin_centers bin_centers(end)+2*pi/num_bins], 'BinCounts', counts, 'FaceColor', colors(end,:));
    % polarhistogram('BinEdges', [bin_centers bin_centers(end)+2*180/num_bins], 'BinCounts', counts, 'FaceColor', colors(end,:), 'Normalization','probability');
    title('preferred theta phase per cell');
    
    % Optional: Adjust the appearance of the plot
    ax = gca;
    ax.RLim = [0 max(counts)];
end

makefigurepretty(gcf,1)
savefigALP([figdir '/'], 'SuppFigure07_D', 'filetype', 'pdf')

%% Supp Fig 7E
figure('Units', 'inches', 'Position', [5.4167 3.1771 2.0521 1.5833]);
% compare vector strength using violins
maxL = max(length(vec_strength_goalstim), length(vec_strength_shamstim));
temp = nan(maxL,2);
temp(1:length(vec_strength_goalstim), 1) = vec_strength_goalstim;
temp(1:length(vec_strength_shamstim), 2) = vec_strength_shamstim;
v1 = violinplot_half(temp);
disp(['there are ', num2str(length(vec_strength_goalstim)), ' goal place cells for vector strength'])
disp(['there are ', num2str(length(vec_strength_shamstim)), ' sham place cells for vector strength'])

colormat = [params.colors_goalstim(2,:); params.colors_shamstim(2,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 25;
end
ylabel('Vector strength')
xlabel('Goal stim (blue) vs Sham stim (orange)')
yticks([0 0.5 1])

makefigurepretty(gcf,1)
savefigALP([figdir '/'], 'SuppFigure07_E', 'filetype', 'pdf')