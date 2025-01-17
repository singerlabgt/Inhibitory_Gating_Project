%% Jeong et al. 2023 MANUSCRIPT - FIGURE 04_B,C,D,E
% NJeong 03/23/2023

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir);

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%load all session information associated with ripple rate data 
load(fullfile(dirs.data2load, 'allindex_ripplerate.mat'));

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','Pyramidal Cell'}; %names used in our manuscript
iBad = cell_metrics.tags.Bad; %bad cell index
iPV = setdiff(cell_metrics.groundTruthClassification.lightsensitive, iBad);
iNarrow = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Narrow')), iBad);
iPyr = setdiff(find(startsWith(cell_metrics.putativeCellType, 'Pyr')), iBad);

%set paramters
downsamplefactor = 15;
ephys_samprate = params.samprate;
binsize_ms = 10;                                        %bin size in ms, trying same binsize as Schlingloff et al., 2014, JNeurosci
binsize_samp = binsize_ms * ephys_samprate / 1000;      %bin size in the number of samples based on sampling rate
secBefore = 10;                                         %number of seconds before SWR-midpoint to look at for SWR-centered FR
secAfter = 10;                                          %number of seconds after SWR-midpoint to look at for SWR-centered FR
speedTh = 2;                                            %speed threshold in deg/s to detect candidate slowing periods
speed_durTh = 5;                                        %find periods that are at least 5 sec long
nontheta_durTh = 5;                                     %include nonthetas >5 sec only
minRips = 10;                                           %minimum number of ripples per session to be included in analysis


%% Ripple Properties Analysis: Ripple power, duration and grouped comparison for WT and PVxAi32 mice

%load SWR info
load(fullfile(dirs.data2load, 'allripples_out.mat'));

%include sessions with at least X ripple periods only
nRips = arrayfun( @(x) size(spikeseries{x},2), 1:length(spikeseries))';
sess2incl = find(nRips >= minRips);

%defined edges for SWR-midpoint-centered firing rates
edges = -secBefore *ephys_samprate : binsize_samp : secAfter*ephys_samprate;

%load place coding data structure 
load(fullfile(dirs.data2load, 'placecodingout_09212024.mat'));


%% create structure to identify place-modulated units in each environment
binEdges = getBinEdges(0:5:360);
nDeg = 10; %in degrees; distance before and after to count as 'goal'

%initialize 
concat_sessInfo = [];       %Nx3, concatenated list of [animal ID, date, novel day]
concat_rzBins.fam = [];     %Nx12, concatenated bin numbers for reward zone in familiar environment
concat_rzBins.nov = [];     %Nx12,concatenated bin numbers for reward zone in novel environment
spatialinfo_fam = [];       %Nx1, concatenated spatial information in familiar environment
spatialinfo_nov = [];       %Nx1, concatenated spatial information in novel environment
isPlaceMod_fam = [];        %Nx1 (binary), 1 if determined as spatially modulated in familiar environment, 0 if not
isPlaceMod_nov = [];        %Nx1 (binary), 1 if determined as spatially modulated in novel environment, 0 if not
famMap = [];                %Nx72, rate map in familiar environment across linearized track
novMap = [];                %Nx72, rate map in novel environment across linearized track

sessions = table2array(unique(allindex(:,1:2),'rows'));
for ii = 1:length(sessions)
    currSess = sessions(ii,:);
    iSess = find(sessions(:,1) == currSess(1) & sessions(:,2) == currSess(2));
    spatialinfo_fam = vertcat(spatialinfo_fam, s_ratemap(iSess).outmap.bin5.fam.spatialinfo_original');
    
    unitID = s_ratemap(iSess).outmap.unitID;
    for iUnit = 1:length(unitID)
        tempInfo = s_ratemap(iSess).outmap.bin5.fam.spatialinfo_original(iUnit);
        threshold_95th = prctile(s_ratemap(iSess).outmap.bin5.fam.spatialinfo_shuffled(:,iUnit), 95); %95th percentile shuffled ratemap
        isPlaceMod_fam = [isPlaceMod_fam; tempInfo > threshold_95th];
        
        if isempty(s_ratemap(iSess).outmap.bin5.nov.ratemap)
            isPlaceMod_nov = [isPlaceMod_nov; nan];
        else
            tempInfo = s_ratemap(iSess).outmap.bin5.nov.spatialinfo_original(iUnit);
            threshold_95th = prctile(s_ratemap(iSess).outmap.bin5.nov.spatialinfo_shuffled(:,iUnit), 95); %95th percentile shuffled ratemap
            isPlaceMod_nov = [isPlaceMod_nov; tempInfo > threshold_95th];
        end
    end
    
    concat_sessInfo = vertcat(concat_sessInfo, repmat(currSess, length(unitID), 1));
    
    famMap = vertcat(famMap, s_ratemap(iSess).outmap.bin5.fam.ratemap_original);
    if isempty(s_ratemap(iSess).outmap.bin5.nov.ratemap)
        s_ratemap(iSess).outmap.bin5.nov.ratemap_original =...
            nan(size(s_ratemap(iSess).outmap.bin5.fam.ratemap_original));
        spatialinfo_nov = vertcat(spatialinfo_nov, nan(size( s_ratemap(iSess).outmap.bin5.fam.spatialinfo_original')));
    else
        spatialinfo_nov = vertcat(spatialinfo_nov, s_ratemap(iSess).outmap.bin5.nov.spatialinfo_original');
    end
    novMap = vertcat(novMap, s_ratemap(iSess).outmap.bin5.nov.ratemap_original);    
    
    
    numUnits = length(s_ratemap(ii).outmap.unitID);
    rz_fam = thetaRZ_fam(ii,:);
    rz_fam = [wrapTo360(rz_fam - nDeg); wrapTo360(rz_fam + nDeg)]';
    temp = arrayfun(@ (x) find(binEdges(:,1)==rz_fam(x,1)),1:3);
    bins_fam = cell2mat( arrayfun(@ (x) temp(x):temp(x)+ (ceil(nDeg*2/5)-1), 1:3, 'UniformOutput', false));
    concat_rzBins.fam = vertcat(concat_rzBins.fam, repmat(bins_fam, numUnits, 1));
    
    rz_nov = thetaRZ_nov(ii,:);
    rz_nov = [wrapTo360(rz_nov - nDeg); wrapTo360(rz_nov + nDeg)]';
    if any(isnan(rz_nov)) %doing this mess bcuz of some fam-only animals
        bins_nov = nan(size(bins_fam));
    else
        temp = arrayfun(@ (x) find(binEdges(:,1)==rz_nov(x,1)),1:3);
        bins_nov = cell2mat( arrayfun(@ (x) temp(x):temp(x)+ (ceil(nDeg*2/5)-1), 1:3, 'UniformOutput', false));
    end
    concat_rzBins.nov = vertcat(concat_rzBins.nov, repmat(bins_nov, numUnits, 1));
end

%identify goal-representing place fields
%count multiple fiels per cell, if available; find units with peak rate
%within goal areas, defined as 10 deg before & after reward zone start angle (i.e. AZ+RZ)
[~, I] = nanmax(famMap,[], 2);
iPeakInFamGoal = find(ismember(I, concat_rzBins.fam)); %indices of PFs with a peak firing in AZ/RZ in at least once
iNoPeakInFamGoal = find(~ismember(I, concat_rzBins.fam)); %indices of PFs with a peak firing outside AZ/RZ
[~, I] = nanmax(novMap,[], 2);
iPeakInNovGoal = find(ismember(I, concat_rzBins.nov));
iNoPeakInNovGoal = find(~ismember(I, concat_rzBins.nov));


clear swr_centered_rate
for ii = 1:length(sess2incl)
    temp = nan(size( spikeseries{sess2incl(ii)},1), length(edges)-1);
    for iUnit = 1:size( spikeseries{sess2incl(ii)},1)
        temp(iUnit,:) = mean( (cell2mat(arrayfun( @(x) ...
            histcounts(spikeseries{sess2incl(ii)}{iUnit,x}, edges),...
            1:size(spikeseries{sess2incl(ii)},2), 'UniformOutput', false)') ./ (binsize_ms./1000)) );
    end
    swr_centered_rate{ii,1} = temp; %each line is a unit averaged over ripples
    swr_centered_rate{ii,2} = allindex(sess2incl(ii),:);
end

%% create structure for plotting 
rate_int_narrow = [];
rate_pyr = [];
sessinfo_int_narrow = [];
sessinfo_pyr = [];
coactivePairs_fam = [];
coactivePairs_nov = [];
%go thru cell array of swr psth and collect interneurons
for ii = 1:length(swr_centered_rate)
    currSess = swr_centered_rate{ii,2};
    sessidx = find(startsWith(cell_metrics.sessionName, [params.iden num2str(currSess.Animal) '_' num2str(currSess.Date)]));
    int_narrow = find(ismember(sessidx, iNarrow));
    pyr = find(ismember(sessidx, iPyr));
        
    %SWR-PSTH dataset
    rate_int_narrow = vertcat( rate_int_narrow, swr_centered_rate{ii,1}(int_narrow,:));
    sessinfo_int_narrow = vertcat( sessinfo_int_narrow, repmat(currSess, length(int_narrow),1));
    rate_pyr = vertcat( rate_pyr, swr_centered_rate{ii,1}(pyr,:));
    sessinfo_pyr = vertcat( sessinfo_pyr, repmat(currSess, length(pyr),1));    
    
    %calculate SWR coactivity probability = number of ripples where two units are coactive, divided by number of total ripples
    goalPyr_fam = intersect(find(ismember(sessidx,iPeakInFamGoal)), pyr);
    goalPyr_nov = intersect(find(ismember(sessidx,iPeakInNovGoal)), pyr);
    
    %in familiar environment, use Pyr. cells with familiar goal modulation
    if currSess.VR == 1 && length(goalPyr_fam) > 1 %had to add this bc nchoosek doesn't work for single number
        famCellPairs = nchoosek(1:length(goalPyr_fam), 2);
        for iPair = 1:size(famCellPairs,1)
            unit1 = famCellPairs(iPair, 1);
            unit2 = famCellPairs(iPair, 2);
            
            numCoactiveRips = length(intersect(...
                find(spikesduringSWR{sess2incl(ii)}( unit1, :)),...
                find(spikesduringSWR{sess2incl(ii)}( unit2, :))) );
            coactivePr_fam = numCoactiveRips ./ size( spikesduringSWR{sess2incl(ii)},2);
            coactivePr_fam = array2table([unit1, unit2, coactivePr_fam], 'VariableNames',{'unit1','unit2','CoactiveProb'});
            coactivePairs_fam = [coactivePairs_fam; [currSess, coactivePr_fam] ];
        end
    end
    
    %in novel environment, use Pyr. cells with novel goal modulation only
    if currSess.VR ~= 1 && length(goalPyr_nov) > 1
        novCellPairs = nchoosek(goalPyr_nov, 2);
        for iPair = 1:size(novCellPairs,1)
            unit1 = novCellPairs(iPair, 1);
            unit2 = novCellPairs(iPair, 2);
            
            numCoactiveRips = length(intersect(find...
                (spikesduringSWR{sess2incl(ii)}( unit1, :)),...
                find(spikesduringSWR{sess2incl(ii)}( unit2, :))) );
            coactivePr_nov = numCoactiveRips ./ size( spikesduringSWR{sess2incl(ii)},2);
            coactivePr_nov = array2table([unit1, unit2, coactivePr_nov], 'VariableNames',{'unit1','unit2','CoactiveProb'});
            coactivePairs_nov = [coactivePairs_nov; [currSess, coactivePr_nov] ];
        end
    end
end

%number of ripples per animal and per session
animals = unique(allindex.Animal); sessions = unique([allindex.Animal, allindex.Date],'rows');
ripples_per_animal = arrayfun( @(x) sum(allindex.ripsInStopped(allindex.Animal == animals(x))), 1:length(animals))';
ripples_per_session = arrayfun( @(x) sum(allindex.ripsInStopped(allindex.Animal == sessions(x,1) & allindex.Date == sessions(x,2))), 1:size(sessions,1))';



%% PLOTTING
fig = figure('units','inch','position',[0 0 7.2 2.7]);
%t = tiledlayout(1,2,'TileSpacing','compact','Units','inches','OuterPosition',[0 0 6.5 2]);
t = tiledlayout(1,6,'TileSpacing','compact');


%% Figure 5B: Sharp-wave ripple rate in novel environment (goal stim vs sham stim)
nexttile; hold on; box off; 
riprate = allindex.ripsInStopped ./ allindex.StoppedDuration; % all sessions included
v_PV.novGoalStim = riprate(ismember(allindex.Animal,params.goalshamMice)...
    & ismember(allindex.Date, sessions(ripples_per_session >= minRips ,2))...
    & allindex.VR ~= 1 & allindex.Stimulation == 1 & allindex.StimLocation == 1 );
v_PV.novShamStim = riprate(ismember(allindex.Animal,params.goalshamMice)...
    & ismember(allindex.Date, sessions(ripples_per_session >= minRips ,2))...
    & allindex.VR ~= 1 & allindex.Stimulation == 1 & allindex.StimLocation == 2 );
maxL = max(length(v_PV.novGoalStim), length(v_PV.novShamStim));
temp = nan(maxL,2);
temp(1:length(v_PV.novGoalStim), 1) = v_PV.novGoalStim;
temp(1:length(v_PV.novShamStim), 2) = v_PV.novShamStim;
v1 = violinplot_half(temp, [], 'CenterSpace', 0.05, 'MedianSize', 20);
colormat = [params.colors_goalstim(2,:); params.colors_shamstim(2,:)];
for ii = 1:length(v1)
    v1(1,ii).ViolinColor = colormat(ii,:);
    v1(1,ii).ScatterPlot.SizeData = 25;
end
ylabel('Ripple rate (Hz)')
xticks([])


%% Figure 5C: Coactivation probability during SWR in novel environment (goal stim vs sham stim)
nexttile; hold on; box off; 
coactive_novgoal = table2array(coactivePairs_nov(...
    ismember(coactivePairs_nov.Animal, params.goalshamMice)...
    & coactivePairs_nov.VR ~= 1 & coactivePairs_nov.Stimulation == 1 ... 
    & coactivePairs_nov.StimLocation == 1, end));
coactive_novsham = table2array(coactivePairs_nov(...
    ismember(coactivePairs_nov.Animal, params.goalshamMice)...
    & coactivePairs_nov.VR ~= 1 & coactivePairs_nov.Stimulation == 1 ...
    & coactivePairs_nov.StimLocation == 2, end));

% downsample the nov sham stim group to match the data size in goal stim
rng(0, "twister")
sampleidx = randi(length(coactive_novgoal), length(coactive_novgoal), 1);
coactive_novgoal = coactive_novgoal(sampleidx);
sampleidx = randi(length(coactive_novsham), length(coactive_novgoal), 1);
coactive_novsham = coactive_novsham(sampleidx);
total_A = length(coactive_novgoal); total_B = length(coactive_novsham);

% exclude 0's - added by XZ 05/19
coactive_novgoal_zeros = coactive_novgoal(coactive_novgoal == 0); coactive_novsham_zeros = coactive_novsham(coactive_novsham == 0); 
perc_zeros = (length(coactive_novsham_zeros) + length(coactive_novgoal_zeros)) /  (length(coactive_novsham) + length(coactive_novgoal_zeros)); %disp(num2str(perc_zeros))
coactive_novgoal_low = coactive_novgoal(coactive_novgoal <= 0.05); coactive_novsham_low = coactive_novsham(coactive_novsham <= 0.05);
perc_low = (length(coactive_novsham_low) + length(coactive_novgoal_low)) /  (length(coactive_novsham) + length(coactive_novgoal_zeros)); %disp(num2str(perc_low))

coactive_novgoal = coactive_novgoal(coactive_novgoal > 0.05);
coactive_novsham = coactive_novsham(coactive_novsham > 0.05);


v_PV.novGoalStim = coactive_novgoal;
v_PV.novShamStim = coactive_novsham; 
maxL = max(height(v_PV.novGoalStim), height(v_PV.novShamStim));
temp = nan(maxL,2);
temp(1:height(v_PV.novGoalStim), 1) = v_PV.novGoalStim;
temp(1:height(v_PV.novShamStim), 2) = v_PV.novShamStim;
v2 = violinplot_half(temp, [], 'CenterSpace', 0.05, 'MedianSize', 20);
colormat = [params.colors_goalstim(3,:); params.colors_shamstim(3,:)];
for ii = 1:length(v2)
    v2(1,ii).ViolinColor = colormat(ii,:);
    v2(1,ii).ScatterPlot.SizeData = 25;
end
ylabel('Coactivation probability')
xticks([])

%%% bar plot showing the count difference of [0 0.05] values for goal/sham
% groups
nexttile; hold on; box off; 
X = categorical({'Goal stim', 'Sham stim'});
X = reordercats(X,{'Goal stim', 'Sham stim'});
data = [length(coactive_novgoal_low), length(coactive_novsham_low)] ./ (length(coactive_novsham) + length(coactive_novgoal_zeros));
b = bar(X, data, 'FaceColor', 'flat', 'EdgeColor','none', 'FaceAlpha',0.8);
b.CData = colormat;
set(gca, 'TickDir', 'out')
ylabel('Proportion of pairs')

% chi-square test
% (https://www.mathworks.com/matlabcentral/answers/96572-how-can-i-perform-a-chi-square-test-to-determine-how-statistically-different-two-proportions-are-in)
% Data for Group A and Group B
positive_A = length(coactive_novgoal_low); positive_B = length(coactive_novsham_low);
% Observed data
n1 = positive_A; N1 = total_A;
n2 = positive_B; N2 = total_B;
% Calculate the observed negatives (O)
O1 = n1;  % Observed positives in Group A
O2 = N1 - n1;  % Observed negatives in Group A
O3 = n2;  % Observed positives in Group B
O4 = N2 - n2;  % Observed negatives in Group B
% Calculate expected counts (E) under the null hypothesis
total = N1 + N2;  % Total number of observations
% Expected counts for Group A (positives and negatives)
E1 = (N1 * (n1 + n2)) / total;  % Expected positives in Group A
E2 = N1 - E1;  % Expected negatives in Group A
% Expected counts for Group B (positives and negatives)
E3 = (N2 * (n1 + n2)) / total;  % Expected positives in Group B
E4 = N2 - E3;  % Expected negatives in Group B
% Calculate Chi-square statistic (𝜒²)
chi2stat = ((O1 - E1)^2 / E1) + ((O2 - E2)^2 / E2) + ((O3 - E3)^2 / E3) + ((O4 - E4)^2 / E4);
% Calculate degrees of freedom (df) for a 2x2 table
df = 1;  % For a 2x2 table, df = (rows - 1) * (columns - 1)
p = 1 - chi2cdf(chi2stat, df);  
phi = sqrt(chi2stat / total);

fprintf('Chi-square Statistic (𝜒²): %.4f\n', chi2stat);
fprintf('Degrees of Freedom (df): %d\n', df);
fprintf('p-value: %.4f\n', p);
fprintf('Effect Size (Phi coefficient): %.4f\n', phi);
fprintf('Expected Counts: Group A Positives = %.2f, Group A Negatives = %.2f\n', E1, E2);
fprintf('Expected Counts: Group B Positives = %.2f, Group B Negatives = %.2f\n', E3, E4);
%% Figure 5D: Ripple power in novel environment (goal stim vs sham stim)
clear temp v_PV
%figure; ax = axes('NextPlot','add','Box','off');
isRZripples = abs(out.animalDistance2RZ) <= 10; %animal position is within 10 degrees around RZ (i.e. AZ+RZ)

totalnovgoalsess = out.sessindex(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==1, :);
totalnovgoalsess_rips = out.rip_maxthreshold(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==1, :);
totalnovshamsess = out.sessindex(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==2, :);
totalnovshamsess_rips = out.rip_maxthreshold(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==2, :);

nexttile; hold on; box off; 

novgoalsess = unique(out.sessindex(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==1, :), 'rows');
novshamsess = unique(out.sessindex(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==2, :), 'rows');


power_novgoalsess = [];
for iSess = 1:size(novgoalsess, 1)
    sessrow = novgoalsess(iSess, :);
    rip_idx_insess = find(ismember(totalnovgoalsess, sessrow, 'rows'));
    power_novgoalsess = [power_novgoalsess, mean(totalnovgoalsess_rips(rip_idx_insess))];
end

power_novshamsess = [];
for iSess = 1:size(novshamsess, 1)
    sessrow = novshamsess(iSess, :);
    rip_idx_insess = find(ismember(totalnovshamsess, sessrow, 'rows'));
    power_novshamsess = [power_novshamsess, mean(totalnovshamsess_rips(rip_idx_insess))];
end

v_PV.novGoalStim = power_novgoalsess;
v_PV.novShamStim = power_novshamsess;
maxL = max(length(v_PV.novGoalStim), length(v_PV.novShamStim));
temp = nan(maxL,2);
temp(1:length(v_PV.novGoalStim), 1) = v_PV.novGoalStim;
temp(1:length(v_PV.novShamStim), 2) = v_PV.novShamStim;
v3 = violinplot_half(temp, [], 'CenterSpace', 0.05, 'MedianSize', 20);
colormat = [params.colors_goalstim(3,:); params.colors_shamstim(3,:)];
for ii = 1:length(v3)
    v3(1,ii).ViolinColor = colormat(ii,:);
    v3(1,ii).ScatterPlot.SizeData = 25;
end
ylabel('Ripple power')
title('PVxAi32 mice using RZ ripples only')
yticks([3 4 5 6])
xticks([])

%%% look at higher quartile -- added by XZ 08082024
nexttile; hold on; box off;
threshold = prctile([v_PV.novGoalStim, v_PV.novShamStim], 50);
power_novgoalsess = [];
for iSess = 1:size(novgoalsess, 1)
    sessrow = novgoalsess(iSess, :);
    rip_idx_insess = find(ismember(totalnovgoalsess, sessrow, 'rows'));
    if mean(totalnovgoalsess_rips(rip_idx_insess)) >= threshold
        power_novgoalsess = [power_novgoalsess, mean(totalnovgoalsess_rips(rip_idx_insess))];
    end
end

power_novshamsess = [];
for iSess = 1:size(novshamsess, 1)
    sessrow = novshamsess(iSess, :);
    rip_idx_insess = find(ismember(totalnovshamsess, sessrow, 'rows'));
    if mean(totalnovshamsess_rips(rip_idx_insess)) > threshold
        power_novshamsess = [power_novshamsess, mean(totalnovshamsess_rips(rip_idx_insess))];
    end
end

v_PV.novGoalStim = power_novgoalsess;
v_PV.novShamStim = power_novshamsess;
maxL = max(length(v_PV.novGoalStim), length(v_PV.novShamStim));
temp = nan(maxL,2);
temp(1:length(v_PV.novGoalStim), 1) = v_PV.novGoalStim;
temp(1:length(v_PV.novShamStim), 2) = v_PV.novShamStim;
v3 = violinplot_half(temp, [], 'CenterSpace', 0.05, 'MedianSize', 20);
colormat = [params.colors_goalstim(3,:); params.colors_shamstim(3,:)];
for ii = 1:length(v3)
    v3(1,ii).ViolinColor = colormat(ii,:);
    v3(1,ii).ScatterPlot.SizeData = 25;
end
ylabel('Ripple power')
% title('PVxAi32 mice using RZ ripples only')
ylim([3.4 4.8])
yticks([3.5 4 4.5])
xticks([])
%% Figure 5E: Ripple duration in novel environment (goal stim vs sham stim)
clear temp v_PV

totalnovgoalsess = out.sessindex(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==1, :);
totalnovgoalsess_rips = out.rip_duration_s(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==1, :);
totalnovshamsess = out.sessindex(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==2, :);
totalnovshamsess_rips = out.rip_duration_s(isRZripples & ...
    ismember(out.sessindex(:,1), params.goalshamMice)...
    & out.sessindex(:,6)~=1 & out.sessindex(:,8)==1 & out.rip_duringStopped == 1 & out.sessindex(:,9)==2, :);
nexttile; hold on; box off; 

ylabel('Ripple duration (s)')
duration_novgoalsess = [];
for iSess = 1:size(novgoalsess, 1)
    sessrow = novgoalsess(iSess, :);
    rip_idx_insess = find(ismember(totalnovgoalsess, sessrow, 'rows'));
    duration_novgoalsess = [duration_novgoalsess, mean(totalnovgoalsess_rips(rip_idx_insess))];
end

duration_novshamsess = [];
for iSess = 1:size(novshamsess, 1)
    sessrow = novshamsess(iSess, :);
    rip_idx_insess = find(ismember(totalnovshamsess, sessrow, 'rows'));
    duration_novshamsess = [duration_novshamsess, mean(totalnovshamsess_rips(rip_idx_insess))];
end
v_PV.novGoalStim = duration_novgoalsess;
v_PV.novShamStim = duration_novshamsess;
maxL = max(length(v_PV.novGoalStim), length(v_PV.novShamStim));
temp = nan(maxL,2);
temp(1:length(v_PV.novGoalStim), 1) = v_PV.novGoalStim;
temp(1:length(v_PV.novShamStim), 2) = v_PV.novShamStim;
v4 = violinplot_half(temp, [], 'CenterSpace', 0.05, 'MedianSize', 20);
ylim([0,1.05]);
colormat = [params.colors_goalstim(3,:); params.colors_shamstim(3,:)];
for ii = 1:length(v4)
    v4(1,ii).ViolinColor = colormat(ii,:);
    v4(1,ii).ScatterPlot.SizeData = 25;
end
ylim([0 0.14])
xticks([])
%% save figure
makefigurepretty(gcf)
figname = 'Figure04_BCDE';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')
