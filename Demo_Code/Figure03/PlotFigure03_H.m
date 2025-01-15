clear; close all;
%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

% load Neuron information (bad, celltype, sessioninfo, etc. for each day)
load(fullfile(maindir, "Demo_Data/", "NeuroInfo.mat"));
% load neurons with significant increase or decrease in AZRZ (for PV
% selection)
load(fullfile(maindir, "Demo_Data/", "allsess_shuffledSig_distance2RZ.mat"));
% load neuron firing rates per day (for PV selection)
load(fullfile(maindir, "Demo_Data/", "allsess_raw_vs_residuals_distance2RZ.mat"));
% load ratemap for each day 
outmapdir = fullfile(maindir, "Demo_Data/", "singlesess_ratemap_distance2RZ/");
%%
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.animals);
allindex = allindex(animal_idx,:); %filter based on animals to include 
allindex = allindex(animal_idx,:); %filter based on animals to include 

[sessions, sessID] = unique(allindex(:,[1:2, 7]),'rows'); %session = [animalID, recording date, novelty day]

%%
correct_only = 1;
decrease_only = 0; % select neurons hvagin significant decrease at goals
separatezone = 0;
scale = 1;

allday_speed = zeros(0, 65);
allday_lickrate = zeros(0, 65);
allday_trialidx = zeros(0, 1);
allday_novenv = zeros(0, 1);
allday_novelexposureday = zeros(0, 1);
allday_iZone = zeros(0, 1);
allday_zonetrial = zeros(0, 1); % trial num on each zone

for iDay = 1:length(NeuronInfo)
    disp(['processing day' num2str(iDay)])
    fn = getlatestfile_with_string(outmapdir, NeuronInfo(iDay).SessionName);
    if length(fn) < 1
        continue
    end
    load(fullfile(outmapdir, fn))
    % identify the bad cells from NeuronInfo
    notbadidx = find(ismember(NeuronInfo(iDay).TSIDtotal, NeuronInfo(iDay).TSID_notbad));
    notbadCluID = NeuronInfo(iDay).CluIDtotal(notbadidx); % check if match allsess
    % identify session related cells in allsess
    sessionname = NeuronInfo(iDay).SessionName;
    sessionname_split = strsplit(sessionname, '_');
    animalnum = str2num(sessionname_split{1}(2:end));
    day = str2num(sessionname_split{2});
    
    % only include WT animals
    if ~ismember(animalnum, params.WTmice)
        continue
    end

    currsess_idx = find(allsess_sessinfo_sh(:, 1)==animalnum & allsess_sessinfo_sh(:, 2)==day);
    if isempty(currsess_idx)
        continue
    end
    if min(allsess_unitID_sh(currsess_idx)' == notbadCluID) == 0 
        disp(['NeuronInfo notbadcluID and allsess_unitID_sh does not match on iDay ', num2str(iDay)])
    end
    % select neurons:
    % not tagged bad, nov decreasing -15-20 deg, high firing (max > 30), NS interneurons
    NS_units = find(ismember(allsess_unitType_sh(currsess_idx), 'Narrow Interneuron'));
    
    if decrease_only == 1
        AZRZ_decrease_units = find(sum(allsess_sigDecrease_fam(currsess_idx, 11:16), 2) > 0);
        unit2incl =  intersect(AZRZ_decrease_units, NS_units);
    else
        unit2incl = NS_units;
    end
    highfiring_units = find(mean(allsess_mean_fam_raw_all(currsess_idx,:), 2) > 20.1099);
    unit2incl = intersect(unit2incl, highfiring_units);
    % if isempty(unit2incl)
    %     continue
    % end
    % select trials:
    if correct_only == 1
        select_trials = find(outmap.bin2.nov.labels(:, 2) == 1 & outmap.bin2.nov.labels(:, 3) == 1); % RZ-centered and receive reward
    else
        select_trials = find(outmap.bin2.nov.labels(:, 2) == 1 & outmap.bin2.nov.labels(:, 3) == 0); % only select RZ_centered
    end
    if isempty(select_trials) || length(select_trials) < 5
        continue
    end
    
    nov_speed = outmap.bin2.nov.smoothSpeed(select_trials, :);
    % normalize
    vMax = nanmax(nov_speed, [], 2);
    vMin = nanmin(nov_speed, [], 2)+0.0001; 
    nov_speed = (nov_speed - vMin) ./ (vMax - vMin);
    nov_speed_az = nanmean(nov_speed(:, ismember(outmap.bin2.binEdges, -10:0)), 2);
    nov_speed_cz = nanmean(nov_speed(:, ismember(outmap.bin2.binEdges, 40:50)), 2);

    nov_lickrate = outmap.bin2.nov.smoothLickrate(select_trials, :);
    % normalize
    vMax = nanmax(nov_lickrate, [], 2);
    vMin = nanmin(nov_lickrate, [], 2)+0.0001;
    nov_lickrate = (nov_lickrate - vMin) ./ (vMax - vMin);
    nov_lickrate_az = nanmean(nov_lickrate(:, ismember(outmap.bin2.binEdges, -10:0)), 2);
    nov_lickrate_cz = nanmean(nov_lickrate(:, ismember(outmap.bin2.binEdges, 40:50)), 2);
    
    
    %% append data to all day table
    allday_speed = [allday_speed; [nov_speed_az, nov_speed_cz]];
    allday_lickrate = [allday_lickrate; [nov_lickrate_az, nov_lickrate_cz]];
    allday_trialidx = [allday_trialidx; (1:size(nov_speed_az, 1))']; 
    allday_novelexposureday = [allday_novelexposureday; repelem(unique(allsess_sessinfo_sh(currsess_idx, 3)), size(nov_speed_az, 1), 1)];
    % get nov env
    index = allindex(allindex(:, 1) == animalnum & allindex(:, 2) == day, :);
    env = unique(index(:, 6));
    if ismember(2, env)
        allday_novenv = [allday_novenv; repelem(2, size(nov_speed_az, 1), 1)];
    elseif ismember(3, env)
        allday_novenv = [allday_novenv; repelem(3, size(nov_speed_az, 1), 1)];
    end
    zones = outmap.bin2.nov.labels(select_trials, 4);
    allday_iZone = [allday_iZone; zones];
    zonetrial = zeros(size(zones));
    for iZone = 1:3
        iZone_trialidx = find(zones == iZone);
        zonetrial(iZone_trialidx) = 1:length(iZone_trialidx);
    end
    allday_zonetrial = [allday_zonetrial; zonetrial];

end

%% ROC plot for first and last trial block. all data pooled together
env = 2; % the env the animal is exposed to novel track for the first time
select_exposureday = 1; % novel day 1
trialblock = 25;
    
day1_select_trials = allday_novelexposureday == select_exposureday & allday_novenv == env;
novelday1_speed = allday_speed(day1_select_trials, :); 
novelday1_lickrate = allday_lickrate(day1_select_trials, :); 
novelday1_zonetrial = allday_trialidx(day1_select_trials);

% first trial block
first_trials = 1:trialblock;
roc_first_trials_speed = calcBehaviorROC(novelday1_speed(ismember(novelday1_zonetrial, first_trials), 1), novelday1_speed(ismember(novelday1_zonetrial, first_trials), 2));
roc_first_trials_lickrate = calcBehaviorROC(novelday1_lickrate(ismember(novelday1_zonetrial, first_trials), 1), novelday1_lickrate(ismember(novelday1_zonetrial, first_trials), 2));

% last trial block
iBlock = 4;
last_trials = trialblock*(iBlock-1)+1 : trialblock*(iBlock);
roc_last_trials_speed = calcBehaviorROC(novelday1_speed(ismember(novelday1_zonetrial, last_trials), 1), novelday1_speed(ismember(novelday1_zonetrial, last_trials), 2));
roc_last_trials_lickrate = calcBehaviorROC(novelday1_lickrate(ismember(novelday1_zonetrial, last_trials), 1), novelday1_lickrate(ismember(novelday1_zonetrial, last_trials), 2));

fig = figure('units','inch','position',[0 0 4 2]);
t = tiledlayout(1, 2);

nexttile; hold on; box off; 
plot(roc_first_trials_speed.X, roc_first_trials_speed.Y, 'Color', params.colors_nov(1, :), 'LineWidth',2);
plot(roc_last_trials_speed.X, roc_last_trials_speed.Y, 'Color', params.colors_nov(3, :), 'LineWidth',2);
plot(0:1,0:1,'LineStyle','--', 'Color','k')
title('speed ROC')

nexttile; hold on; box off; 
plot(roc_first_trials_lickrate.X, roc_first_trials_lickrate.Y, 'Color', params.colors_nov(1, :), 'LineWidth',2);
plot(roc_last_trials_lickrate.X, roc_last_trials_lickrate.Y, 'Color', params.colors_nov(3, :), 'LineWidth',2);
plot(0:1,0:1,'LineStyle','--', 'Color','k')
title('lickrate ROC')

makefigurepretty(gcf,1)
figname = 'Figure03_H_ROC';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')