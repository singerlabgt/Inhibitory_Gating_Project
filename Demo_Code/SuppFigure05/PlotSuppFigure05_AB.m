%% Population Aaverage Time2AZ Ratemap across Stim Intensities
clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end


% get all indices 
[allindex, ~] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
[sessions, ind] = unique(allindex(:,[1:2,7]), 'rows'); %define session as one date



%define binsize for firing activity and behavior (in ms; same as fieldnames)
spikeBinSz = 50;
behavBinSz = 200;

% %get all directories and parameters
% [dirs, params] = getDirectoriesAndParams();
% datesincl = [];
% datesexcl = [];
% 
% allindexT = selectindextable(dirs.spreadsheet_rec, 'animal', params.animals, ...
%     'datesincluded', datesincl, 'datesexcluded', datesexcl);
% allindexT = allindexT(ismember(allindexT.VR, 1:3),:); %include active behavior only with no VR manipulations
% allindex = allindexT{:,{'Animal','Date','Recording','Include','Depth','VR','NovelExposure','Stimulation','StimLocation'}};
% [sessions, ind] = unique(allindex(:,[1:2,7]), 'rows'); %define session as one date
% identifier = allindexT{ind,'ID'};

%% create a population average ratemap output structure of N x M size;
%% N units x M time bins centered at AZ/stim onset (2sec before and after)

%load cell type info
load(fullfile(maindir, 'Demo_Data/', 'cell_metrics.mat'))
iBad = cell_metrics.tags.Bad; %bad units to exclude
celltypes = {'PV Interneuron','Narrow Interneuron','Wide Interneuron','Pyramidal Cell'};
for ct = 1:length(celltypes)
    iCT{ct} = setdiff(find(strcmp(cell_metrics.putativeCellType,celltypes{ct})), iBad)';
    if ct==1
        lightsensitive = cell_metrics.groundTruthClassification.lightsensitive;
        narrow = setdiff(find(startsWith(cell_metrics.putativeCellType,'Narrow')), iBad)';
        iCT{ct} = intersect(lightsensitive, narrow);
    end
end

%determine [animalID, recDay, novDay] associated with each unit for indexing
numUnits = length(cell_metrics.sessionName);
cellMetSessInfo = nan(numUnits,3);
splitStr = regexp(cell_metrics.sessionName,'_','split'); %doing it this way bc different identifier letters and single vs double digit animal number
cellMetSessInfo(:,1) = arrayfun( @(x) str2double(extractAfter(splitStr{x}{1}, 1)), 1:numUnits);
cellMetSessInfo(:,2) = arrayfun( @(x) str2double(splitStr{x}{2}), 1:numUnits); 
cellMetSessInfo(:,3) = arrayfun( @(x) sessions(cellMetSessInfo(x,1) == sessions(:,1) & cellMetSessInfo(x,2) == sessions(:,2), 3), 1:numUnits);

%variables
stimInt = [0, 0.5, 1.6, 3.2];
stimNames = {'None','Low','High','Highest'};
stimTypes = {'NovGoalStim','NovShamStim','FamGoalStim','NovAZStim','FamAZStim'};
frFname = ['bin' num2str(spikeBinSz)]; %spike fieldname based on chosen binsize
bFname = ['bin' num2str(behavBinSz)]; %behavior fieldname based on chosen binsize


%% make a giant output structure:
%% for FR/speed/lickrate across position bins: N units (or trials) x M time bins centered at AZ/stim onset (2sec before and after)
%% for quantification ("q"): N x 1 with N being the number of units or trials -- this is the mean FR/speed/lickR averaged over the 2-second AZ/stim bins

%output structure fieldnames
fdNames = {'FRrawxPos','deltaFRraw','deltaFRnorm','qFR_raw','qFR_norm','unitID',...
    'SpeedxPos_raw','SpeedxPos_norm','qSpeed_raw','qSpeed_norm',...
    'LickRxPos_raw','LickRxPos_norm','qLickR_raw','qLickR_norm','sessInfo'};
for sInt = 1:length(stimNames)
    voltName = ['Volt_' stimNames{sInt}];
    for ee = 1:length(params.environments)
        env = params.environments{ee};
        for f = 1:length(fdNames)
            AllMap.(env).(fdNames{f}).(voltName) = [];
        end
    end
end

%go thru each session
nov_none_day1 = 0;
nov_none_day2 = 0;
nov_none_day3 = 0;
nov_low_day1 = 0;
nov_low_day2 = 0;
nov_low_day3 = 0;
nov_high_day1 = 0;
nov_high_day2 = 0;
nov_high_day3 = 0;
fam_none_day3 = 0;
fam_low_day3 = 0;
fam_high_day3 = 0;
for iSess = 1:length(sessions)
    currSess = sessions(iSess, :);
    animal = currSess(1);
    recDay = currSess(2);
    index = allindex(allindex(:,1) == animal & allindex(:,2) == recDay, :);
    novDay = unique(index(:, 7));
    
    %figure out if any sort of stim was used
    novGoalStim = 0;     novShamStim = 0;     famGoalStim = 0; %zero if not
    if sum(index(:,6) ~= 1 & index(:,9)==1) > 0 %see if novgoalstim was used
        novGoalStim = 1;
    elseif sum(index(:,6) ~= 1 & index(:,9)==2) > 0
        novShamStim = 1;
    elseif sum(index(:,6) == 1 & index(:,9)==1) > 0
        famGoalStim = 1;
    end
    
    %load dist2RZ ratemap file for timewindow and stim info for each trial
    if animal == 4
        iden = 'X';
    else
        iden = 'N';
    end
    fdir = fullfile(maindir, 'Demo_Data', 'singlesess_ratemap_time2AZ/');
    time2AZfname = getlatestfile_with_string(...
        fdir,[iden num2str(animal) '_' num2str(recDay)]);
    timeMap = load(fullfile(fdir,time2AZfname));
    timeMap = timeMap.outmap;
    unitID = timeMap.unitID;
    
    %get binEdges based on bin sizes - same for both environments
    fdNames = fieldnames(timeMap); fdNames = fdNames(startsWith(fdNames,'bin'));
    binEdges_FR = timeMap.(frFname).fam.binEdges ./ params.samprate; %for firing rates
    binEdges_B = timeMap.(bFname).fam.binEdges ./ params.samprate; %for behavior, separating this bc usually larger binsize
    
    for ee = 1:length(params.environments) %for 1-fam and 2-nov environments
        env = params.environments{ee}; 
        
        %trialV = stim voltage used per trial, gives numeric value in V
        trialV = arrayfun( @(x) unique(...
            timeMap.(frFname).(env).stimVoltage{x}(find(timeMap.(frFname).(env).stimVoltage{x})),'stable'),...
            1:length(timeMap.(frFname).(env).stimVoltage), 'UniformOutput', false);
        trialV(find(arrayfun( @(x) isempty(trialV{x}), 1:length(timeMap.(frFname).(env).stimVoltage)))) = {0}; %put no stim as zero
        trialV = arrayfun( @(x) trialV{x}(1), 1:length(timeMap.(frFname).(env).stimVoltage))';
        
        %get 1 FR x Time2AZ averaged across trials per neuron per stim intensity
        for sInt = 1:length(stimInt)
            trialInt = lookup2(trialV,stimInt);
            
            tempFR = nan(length(timeMap.unitID), length(binEdges_FR)-1);
            tempFR_norm = nan(length(timeMap.unitID), length(binEdges_FR)-1);
            tempSpeed = []; speedRaw = []; speedNorm = [];
            tempLickR = []; lickRaw = []; lickNorm = [];
            
            %put NaN if there is no novel exposure; this is to keep the
            %matrix size the same as cell_metrics.mat
            if isempty(timeMap.(frFname).(env).ratemap)
                centerTrials = [];
            elseif ee==2 && novShamStim == 1
                %for sham stim, get NRZ-centered trials
                centerTrials = timeMap.(frFname).(env).labels(:,2) == 0; %zero here refers to NRZ-centered trials (i.e. Non-AZ centered trials)
            else
                %for all other stim or no-stim sessions, include AZ-centered trials only
                centerTrials = timeMap.(frFname).(env).labels(:,2) == 1;
            end
            
            %find right trials to include with correct stim type
            trials2incl = intersect(find(centerTrials), find(trialInt == sInt));
            
            %mean firing rate per unit per stim intensity
            if length(trials2incl) > 5 %only attempt to average if enough trials
                if ee == 1 && sInt == 1 && novDay == 3
                    fam_none_day3 = fam_none_day3 + length(trials2incl);
                elseif ee == 1 && sInt == 2 && novDay == 3
                    fam_low_day3 = fam_low_day3 + length(trials2incl);
                elseif ee == 1 && sInt == 3 && novDay == 3
                    fam_high_day3 = fam_high_day3 + length(trials2incl);
                elseif ee == 2 
                    if sInt == 1
                        if novDay == 1
                            nov_none_day1 = nov_none_day1 + length(trials2incl);
                        elseif novDay == 2
                            nov_none_day2 = nov_none_day2 + length(trials2incl);
                        elseif novDay == 3
                            nov_none_day3 = nov_none_day3 + length(trials2incl);
                        end
                    elseif sInt == 2
                        if novDay == 1
                            nov_low_day1 = nov_low_day1 + length(trials2incl);
                        elseif novDay == 2
                            nov_low_day2 = nov_low_day2 + length(trials2incl);
                        elseif novDay == 3
                            nov_low_day3 = nov_low_day3 + length(trials2incl);
                        end
                    elseif sInt == 3
                        if novDay == 1
                            nov_high_day1 = nov_high_day1 + length(trials2incl);
                        elseif novDay == 2
                            nov_high_day2 = nov_high_day2 + length(trials2incl);
                        elseif novDay == 3
                            nov_high_day3 = nov_high_day3 + length(trials2incl);
                        end
                    end
                end
                for iUnit = 1:length(timeMap.unitID)
                    tempFR(iUnit,:) = nanmean(timeMap.(frFname).(env).ratemap(trials2incl,:,iUnit));
                end
            end
            
            %specify FR bins to use for baseline and during-stimulation period
            temp = getBinEdges(binEdges_FR);
            baseT = find(ismember(temp(:,1), -1:0)); baseT = baseT(1) : baseT(end)-1; %baseline: avg of [-1 to 0s) prior to stim onset
            duringST = find(ismember(temp(:,1), 0:2)); duringST = duringST(1) : duringST(end)-1;
            
            %calcualte mean change in firing rate, averaged over 2s post-stim
            deltaRaw = tempFR - mean(tempFR(:,baseT),2); %raw change from baseline
            frMeanRaw = nanmean(deltaRaw(:,duringST),2); %one per unit
            
            %calcualte mean change in Normalized firing rate, averaged over 2s post-stim
            deltaNorm = tempFR ./ max(tempFR,[],2); %first normalize by peakFR per unit
            deltaNorm = deltaNorm - mean(deltaNorm(:,baseT),2); %normalized change from baseline
            frMeanNorm = nanmean(deltaNorm(:,duringST),2); %one per unit
            
            %behavior per stim intensity (not averaged across trials)
            %may need to use different bins for baseline and during stim
            %for behavior than spike data
            temp = getBinEdges(binEdges_B);
            baseT = find(ismember(temp(:,1), -1:0)); baseT = baseT(1) : baseT(end)-1; %baseline: avg of [-1 to 0s) prior to stim onset
            duringST = find(ismember(temp(:,1), 0:2)); duringST = duringST(1) : duringST(end)-1;
            
            if length(trials2incl) > 5
                tempSpeed = timeMap.(bFname).(env).smoothSpeed(trials2incl,:);
                deltaSpeedRaw = tempSpeed - mean(tempSpeed(:,baseT),2); %raw change from baseline
                speedRaw = nanmean(deltaSpeedRaw(:,duringST),2);
                deltaSpeedNorm = tempSpeed ./ max(tempSpeed,[],2); %normalize by peak
                deltaSpeedNorm = deltaSpeedNorm - mean(deltaSpeedNorm(:,baseT),2);
                speedNorm = nanmean(deltaSpeedNorm(:,duringST),2);
                
                tempLickR = timeMap.(bFname).(env).smoothLickrate(trials2incl,:);
                deltaLickRaw = tempLickR - mean(tempLickR(:,baseT),2); %raw change from baseline
                lickRaw = nanmean(deltaLickRaw(:,duringST),2);
                deltaLickNorm = tempLickR ./ max(tempLickR,[],2); %normalize by peak
                deltaLickNorm = deltaLickNorm - mean(deltaLickNorm(:,baseT),2);
                lickNorm = nanmean(deltaLickNorm(:,duringST),2);
                
                AllMap.(env).SpeedxPos_raw.(['Volt_' stimNames{sInt}]) = [AllMap.(env).SpeedxPos_raw.(['Volt_' stimNames{sInt}]); deltaSpeedRaw];
                AllMap.(env).SpeedxPos_norm.(['Volt_' stimNames{sInt}]) = [AllMap.(env).SpeedxPos_norm.(['Volt_' stimNames{sInt}]); deltaSpeedNorm];
                AllMap.(env).LickRxPos_raw.(['Volt_' stimNames{sInt}]) = [AllMap.(env).LickRxPos_raw.(['Volt_' stimNames{sInt}]); deltaLickRaw];
                AllMap.(env).LickRxPos_norm.(['Volt_' stimNames{sInt}]) = [AllMap.(env).LickRxPos_norm.(['Volt_' stimNames{sInt}]); deltaLickNorm];
                
                fdNames = {'qSpeed_raw','qSpeed_norm','qLickR_raw','qLickR_norm'}; %behavior fieldnames
                temp = [speedRaw, speedNorm, lickRaw, lickNorm];
                if ~isempty(temp)
                    for f = 1:length(fdNames)
                        AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}]) = [AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}]); temp(:,f)];
                    end
                end
            end
            %gather population data for each stim intensity
            AllMap.(env).FRrawxPos.(['Volt_' stimNames{sInt}]) = [AllMap.(env).FRrawxPos.(['Volt_' stimNames{sInt}]); tempFR];
            AllMap.(env).deltaFRraw.(['Volt_' stimNames{sInt}]) = [AllMap.(env).deltaFRraw.(['Volt_' stimNames{sInt}]); deltaRaw];
            AllMap.(env).deltaFRnorm.(['Volt_' stimNames{sInt}]) = [AllMap.(env).deltaFRnorm.(['Volt_' stimNames{sInt}]); deltaNorm];
            fdNames = {'qFR_raw','qFR_norm'}; %FR fieldnames
            temp = [frMeanRaw, frMeanNorm];
            for f = 1:length(fdNames)
                AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}]) = [AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}]); temp(:,f)];
            end
            AllMap.(env).unitID.(['Volt_' stimNames{sInt}]) = [AllMap.(env).unitID.(['Volt_' stimNames{sInt}]); unitID];
            
            
            AllMap.(env).sessInfo.(['Volt_' stimNames{sInt}]) = [AllMap.(env).sessInfo.(['Volt_' stimNames{sInt}]); repmat([animal,recDay,novDay],length(speedNorm),1)];
        end
    end
    clear timeMap
end
manuscript_N = [nov_none_day1, nov_low_day1, nov_high_day1, nov_none_day2, nov_low_day2, nov_high_day2...
    nov_none_day3, nov_low_day3, nov_high_day3, fam_none_day3, fam_low_day3, fam_high_day3];
%save output
AllMap.binsize_firingrate = spikeBinSz;
AllMap.binsize_behavior = behavBinSz;
AllMap.binEdges_firingrate = binEdges_FR; %in seconds, 0 = time of light onset
AllMap.binEdges_behavior = binEdges_B;  %in seconds, 0 = time of light onset


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% create FR data structure for plotting
% clear data; fdNames = {'qFR_raw','qFR_norm'};
% for f = 1:length(fdNames)
%     data.(fdNames{f}) = [];
%     for sInt = 1:length(stimNames)
%         for ee = 1:length(params.environments) %for fam and nov environments
%             env = params.environments{ee};
%             for ct = 1:length(iCT)
%                 data2add = AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}])(iCT{ct},:);
%                 unitID = AllMap.(env).unitID.(['Volt_' stimNames{sInt}])(iCT{ct},:);
%                 sessInfo = cellMetSessInfo(iCT{ct},:);
%                 cellT = ones(length(unitID),1) .* ct; %CellT of 1-4 indicate PV, narrow, wide, pyr, in that order
%                 Environ = ones(length(unitID),1) .* ee;
%                 Intensity = ones(length(unitID),1) .* sInt;
%                 data2add = [data2add, sessInfo, unitID, cellT, Environ, Intensity]; %
%                 data.(fdNames{f}) = [data.(fdNames{f}); data2add];
%             end
%         end
%     end
%     temp = data.(fdNames{f});
%     data.(fdNames{f}) = temp(~isnan(temp(:,1)),:);
% end
% 
% 
% fdNames = {'deltaFRraw','deltaFRnorm'};
% for f = 1:length(fdNames)
%     data.(fdNames{f}) = [];
%     for sInt = 1:length(stimNames)
%         for ee = 1:length(params.environments) %for fam and nov environments
%             env = params.environments{ee};
%             for ct = 1:length(iCT)
%                 sessInfo = [];
%                 data2add = AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}])(iCT{ct},:);
%                 unitID = AllMap.(env).unitID.(['Volt_' stimNames{sInt}])(iCT{ct},:);
%                 sessInfo = cellMetSessInfo(iCT{ct},:);
%                 cellT = ones(length(unitID),1) .* ct; %CellT of 1-4 indicate PV, narrow, wide, pyr, in that order
%                 Environ = ones(length(unitID),1) .* ee;
%                 Intensity = ones(length(unitID),1) .* sInt;
%                 data2add = data2add(all(~isnan(data2add),2),:); %excl NaN
%                 data.(fdNames{f}) = [data.(fdNames{f}); data2add];
%             end
%         end
%     end
% end


%% PLOT STIM EFFECTS ON FIRING RATES PER INTENSITY: ALL SESSIONS COMBINED FROM ALL ANIMALS
%identify session info
stim_famgoal = unique(allindex(allindex(:,6) == 1 & allindex(:,9) == 1, [1:2,7]),'rows'); %includes all mice
stim_novgoal = unique(allindex(allindex(:,6) ~= 1 & allindex(:,9) == 1, [1:2,7]),'rows');
stim_novsham = unique(allindex(allindex(:,6) ~= 1 & allindex(:,9) == 2, [1:2,7]),'rows');



%% behavior: Day x stim intensity
% create behavior data structure for plotting
xTime = mean(getBinEdges(binEdges_B),2);

clear data; fdNames = {'qSpeed_raw','qSpeed_norm','qLickR_raw','qLickR_norm'};

for f = 1:length(fdNames)
    data.(fdNames{f}) = [];
    for sInt = 1:length(stimNames)
        for ee = 1:length(params.environments) %for fam and nov environments
            env = params.environments{ee};
            data2add = AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}]);
            sessInfo = AllMap.(env).sessInfo.(['Volt_' stimNames{sInt}]);
            Environ = ones(length(data2add),1) .* ee;
            Intensity = ones(length(data2add),1) .* sInt;
            data2add = [data2add, sessInfo, Environ, Intensity]; %
            data.(fdNames{f}) = [data.(fdNames{f}); data2add];
            
        end
    end
    temp = data.(fdNames{f});
    data.(fdNames{f}) = temp(~isnan(temp(:,1)),:);
end

fdNames = {'SpeedxPos_raw','SpeedxPos_norm','LickRxPos_raw','LickRxPos_norm'};
for f = 1:length(fdNames)
    data.(fdNames{f}) = [];
    for sInt = 1:length(stimNames)
        for ee = 1:length(params.environments) %for fam and nov environments
            env = params.environments{ee};
            data2add = AllMap.(env).(fdNames{f}).(['Volt_' stimNames{sInt}]);
            sessInfo = AllMap.(env).sessInfo.(['Volt_' stimNames{sInt}]);
            data2add = data2add(all(~isnan(data2add),2),:); %excl NaN
            data.(fdNames{f}) = [data.(fdNames{f}); data2add];
            
        end
    end
end



behNames ={'Speed','LickR'};
fdNames = {'norm'};
% fdNames = {'raw','norm'};
for b = 1:length(behNames)
    clear behave
    for f = 1:length(fdNames)
        datTrace = data.([behNames{b}  'xPos_' fdNames{f}]);
        datQ = data.(['q' behNames{b} '_' fdNames{f}])(:,1);
        sessInfo = data.(['q' behNames{b} '_' fdNames{f}])(:,2:end);
        iNovGoalStim = ismember(sessInfo(:,1), params.goalshamMice) & ismember(sessInfo(:,2), stim_novgoal(:,2)) & sessInfo(:,4) == 2;
        iNovShamStim = ismember(sessInfo(:,1), params.goalshamMice) & ismember(sessInfo(:,2), stim_novsham(:,2)) & sessInfo(:,4) == 2;
        iFamGoalStim = ismember(sessInfo(:,1), params.goalshamMice) & ismember(sessInfo(:,2), stim_famgoal(:,2)) & sessInfo(:,4) == 1;

        fig = figure('units','inch','position',[0 0 6.5 4]);
        t = tiledlayout(2,4,'TileSpacing','compact','Units','inches','OuterPosition',[0 0 6.5 4]);
%         t = tiledlayout(2,4,'TileSpacing','compact');
        clear h
        
        for nD = 1:4
            h(nD) = nexttile; box off;
            if nD < 4
                for sInt = 1:length(stimInt)-1
                    %nov behavior
                    behave{nD}{sInt} = datQ(iNovGoalStim & sessInfo(:,3) == nD & sessInfo(:,5) == sInt, :);
                    
                    ops.ax = h(nD);
                    ops.error = 'sem';
                    ops.x_axis = xTime;
                    ops.alpha = 0.2;
                    ops.line_width = 2;
                    ops.color_area = params.colors_nov(sInt,:);
                    ops.color_line = params.colors_nov(sInt,:);
                    plot_areaerrorbar(datTrace(iNovGoalStim & sessInfo(:,3) == nD & sessInfo(:,5) == sInt, :), ops);
                    box off; hold on;
                    title(h(nD), ['Nov ' behNames{b} ' ' num2str(nD)])
                    axis square
                    xline(h(nD),0,'k:'); yline(h(nD),0,'k:');
                    
                    %save stats info
                    % fileID = fopen(fullfile(figdir, 'stats.txt'), 'a');
                    %fprintf(fileID,'%s\r\n', deblank([fdNames{f} 'Nov' behNames{b} '_Day' num2str(nD) '_stimIntensity_' stimNames{sInt} ' [min, 25th, 50th, 75th, max] = '...
                        % num2str( prctile(behave{nD}{sInt}, [0,25,50,75,100]),'% 5.4f') ]));
                    %fprintf(fileID,'%s\r\n', deblank([fdNames{f} 'Nov' behNames{b} '_Day' num2str(nD) '_stimIntensity_' stimNames{sInt} ' n = ' num2str(length(behave{nD}{sInt}))...
                        % ', ' num2str(mean(behave{nD}{sInt})) ' ± ' num2str(std(behave{nD}{sInt}) ./ sqrt(length(behave{nD}{sInt})))]));
                end
            else
                for sInt = 1:length(stimInt)-1
                    %fam behavior
                    behave{nD}{sInt} = datQ(iFamGoalStim & sessInfo(:,3) == 3 & sessInfo(:,5) == sInt, :);
                    
                    ops.ax = h(nD);
                    ops.error = 'sem';
                    ops.x_axis = xTime;
                    ops.alpha = 0.2;
                    ops.line_width = 2;
                    ops.color_area = params.colors_fam(sInt,:);
                    ops.color_line = params.colors_fam(sInt,:);
                    plot_areaerrorbar(datTrace(iFamGoalStim & sessInfo(:,3) == 3 & sessInfo(:,5) == sInt, :), ops);
                    box off; hold on;
                    title(h(nD), ['Fam '  behNames{b} ' 3'])
                    axis square
                    xline(h(nD),0,'k:'); yline(h(nD),0,'k:');
                    
                    %save stats info
                    % fileID = fopen(fullfile(figdir, 'stats.txt'), 'a');
                    %fprintf(fileID,'%s\r\n', deblank([fdNames{f} '_Fam' behNames{b} '_Day' num2str(nD) '_stimIntensity_' stimNames{sInt} ' [min, 25th, 50th, 75th, max] = '...
                        % num2str( prctile(behave{nD}{sInt}, [0,25,50,75,100]),'% 5.4f') ]));
                    %fprintf(fileID,'%s\r\n', deblank([fdNames{f} '_Fam' behNames{b} '_Day' num2str(nD) '_stimIntensity_' stimNames{sInt} ' n = ' num2str(length(behave{nD}{sInt}))...
                        % ', ' num2str(mean(behave{nD}{sInt})) ' ± ' num2str(std(behave{nD}{sInt}) ./ sqrt(length(behave{nD}{sInt})))]));
                    
                end
            end
        end
        
        for nD = 1:4
            h(nD+4) = nexttile; hold on;
            if nD < 4
                arrayfun( @(x) scatter(h(nD+4), ones(size(behave{nD}{x},1),1) .*x, behave{nD}{x}, 50, 'Jitter', 'on',...
                    'JitterAmount', 0.1, 'MarkerEdgeColor', params.colors_nov(x,:),...
                    'MarkerFaceColor', params.colors_nov(x,:),'MarkerFaceAlpha', 0.2), 1:3);
                errorbar(h(nD+4), 1:3, arrayfun( @(x) mean(behave{nD}{x}), 1:3),...
                    arrayfun( @(x) std(behave{nD}{x}), 1:3) ./ arrayfun( @(x) sqrt(length(behave{nD}{x})), 1:3),'k-')
                axis square
            else
                arrayfun( @(x) scatter(h(nD+4), ones(size(behave{nD}{x},1),1) .*x, behave{nD}{x}, 50, 'Jitter', 'on',...
                    'JitterAmount', 0.1, 'MarkerEdgeColor', params.colors_fam(x,:),...
                    'MarkerFaceColor', params.colors_fam(x,:),'MarkerFaceAlpha', 0.2), 1:3);
                errorbar(h(nD+4), 1:3, arrayfun( @(x) mean(behave{nD}{x}), 1:3),...
                    arrayfun( @(x) std(behave{nD}{x}), 1:3) ./ arrayfun( @(x) sqrt(length(behave{nD}{x})), 1:3),'k-')
                axis square
            end
            xlim(h(nD+4), [0.5 3.5]); xticks(h(nD+4), 1:3); xticklabels(h(nD+4), {'None','Low','High'});
        end
        
        xlim(h(1:4), [-1,2])
        linkaxes(h(1:4),'xy')
        linkaxes(h(5:8),'xy')
        xlabel(h(1),'Time to Stim/AZ onset (s)')
        ylabel(t,[fdNames{f} ' change in '  behNames{b} ' from baseline'])
        if b == 1
            figname = 'SuppFig06_B';
        elseif b == 2
            figname = 'SuppFig06_A';
        end
        makefigurepretty(gcf)
        savefigALP([figdir '/'], figname, 'filetype', 'pdf')
        %exportgraphics(t, fullfile(figdir, [fdNames{f} behNames{b} '_DayXIntensity.pdf']),...
            %'Resolution',300,'ContentType','vector','BackgroundColor','none');
%         set(gcf,'Position',get(0,'ScreenSize'),'PaperOrientation','landscape')
%         %print(gcf,figname,'-dpdf','-r300','-bestfit'); clf;
    end
end

