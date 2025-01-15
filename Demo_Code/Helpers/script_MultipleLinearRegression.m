%% Multiple linear regression
% This script takes binned firing rate maps (binned based on distance or
% time) and performs multiple linear regression using binned speed and lick
% rate behavioral maps as two predictors. The main purpose of the script is
% to regress out potential effects of position/time-related changes in
% behavior on the observed firing rates, and to create a residual firing
% rate map that cannot be attributed to simple changes in behavior. The
% outputs are 1) trial-by-trial residual firing rates of the same size as
% original binned rate map structure, and 2) averaged residual firing rates
% across trials per unit, using at least 5 correct trials only. 
% NJ created 05/19/2022
% NJ edited 07/12/2022 to list all trial labels instead of pre-separating
% them by trial type
% NJ edited 11/09/2022 to use cell_metrics structure from CellExplorer as
% celltype info and combine int/pyr into one output structure with celltype
% labels; saves unitIDs; saves population average data in the same folder


clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\test_logs\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir);

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));

%flags to create distance- or time-based residual firing rate maps 
regressLaps = 1;    %regress out firing rate maps for the entire track from 0-360 degrees (distance-based; in degrees)
regressDist2RZ = 1; %regress out firing rate maps for the areas around the reward zone, 60 degrees before and after the RZ (distance-based; in 5 deg increments)
regressTime2RZ = 1; %regress out firing rate maps for the time around the reward zone, 6 seconds before and after the RZ entry (time-based; in 100 ms increments)
plotScatter_laps = 0;
%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);
[sets, ~,~] = splitSessions2Set(allindex);
sets = sets(arrayfun( @(x) size(sets{x},1), 1:length(sets)) == 3); %only include files with 3-day fam/nov paradigm, NJ added 5/16/22
uniqSess = cell2mat(sets);


%% create an output structure of RZ/NRZ-centered ratemap across reward trials, separated by stim vs no stim
if regressDist2RZ
    %running structure across multiple sessions
    allsess_sessinfo = []; allsess_unitID = []; allsess_unitType = []; allsess_stats = [];
    allsess_mean_fam_raw_all = []; allsess_mean_fam_raw_hit = []; allsess_mean_fam_raw_miss = [];
    allsess_mean_fam_raw_1st = []; allsess_mean_fam_raw_2nd = [];
    allsess_mean_nov_raw_all = []; allsess_mean_nov_raw_hit = []; allsess_mean_nov_raw_miss = [];
    allsess_mean_nov_raw_1st = []; allsess_mean_nov_raw_2nd = [];
    allsess_mean_fam_residual_all = []; allsess_mean_fam_residual_hit = []; allsess_mean_fam_residual_miss = [];
    allsess_mean_fam_residual_1st = []; allsess_mean_fam_residual_2nd = [];
    allsess_mean_nov_residual_all = []; allsess_mean_nov_residual_hit = []; allsess_mean_nov_residual_miss = [];
    allsess_mean_nov_residual_1st = []; allsess_mean_nov_residual_2nd = [];
    
    for iSess = 1:length(uniqSess)
        currSess = uniqSess(iSess, :);
        animal = currSess(1);
        recDay = currSess(2);
        novelDay = currSess(3);
        sessStr = [params.iden num2str(animal) '_' num2str(recDay)];
        sessRegStr = [params.iden num2str(animal) '_' num2str(recDay) '_' params.brainReg{1}];
        
        %load binned firing rate maps with behavioral info
        fdir = fullfile(dirs.data2load, 'singlesess_ratemap_distance2RZ');        
        dist2RZfname = getlatestfile_with_string(fdir,...
            [sessStr '_createTrialByTrialDist2RewardRateMaps_Day ' num2str(novelDay)]);        
        distMap = load(fullfile(fdir,dist2RZfname));
        distMap = distMap.outmap;
        unitID = distMap.unitID;
        
        %% get cell type info
        iCellInMetrics = find(...
            arrayfun( @(x) strcmp(cell_metrics.sessionName{x}, sessRegStr),...
            1:length(cell_metrics.sessionName)));
        
        %check same unitIDs between outmap and cell_metrics structure, as
        %expected. 'idx_all' is 1:N within a single session, not the entire
        %cluster_metrics indices
        if sum(~ismember(unitID, cell_metrics.cluID(iCellInMetrics))) == 0
            idx_all = cell_metrics.cellID(iCellInMetrics(...
                ~ismember(iCellInMetrics, cell_metrics.tags.Bad))); %exclude units tagged as 'Bad'
            unitID = unitID(idx_all);
            unitTypeLabel = cell_metrics.putativeCellType(...
                iCellInMetrics(~ismember(iCellInMetrics, cell_metrics.tags.Bad)))';
        end
        
        
        %% collect behavioral data (smoothed speed & smoothed lick rate)
        position_binEdges = distMap.bin5.binEdges;
        
        famFiles    = distMap.sessionInfo(distMap.sessionInfo(:,6) == 1, 3);
        famCount    = distMap.bin5.fam.rawCount;
        famOccup    = distMap.bin5.fam.rawOccup;
        famSpeed    = distMap.bin5.fam.smoothSpeed;
        famLickrate = distMap.bin5.fam.smoothLickrate;
        
        novFiles    = distMap.sessionInfo(distMap.sessionInfo(:,6) ~= 1, 3);
        novCount    = distMap.bin5.nov.rawCount;
        novOccup    = distMap.bin5.nov.rawOccup;
        novSpeed    = distMap.bin5.nov.smoothSpeed;
        novLickrate = distMap.bin5.nov.smoothLickrate;
        
        %collect trial type and stim intensity info (in volts)
        famTrialLabel = distMap.bin5.fam.labels;
        famTrialVoltage = distMap.bin5.fam.stimVoltage;
        novTrialLabel = distMap.bin5.nov.labels;
        novTrialVoltage = distMap.bin5.nov.stimVoltage;
        trialLabelDict = distMap.bin5.fam.labels_dict;
        
        %organize X for fam and nov behavior that is consistent for all cells
        x1 = reshape([famSpeed; novSpeed],[],1); %predictor 1: speed
        x2 = reshape([famLickrate; novLickrate],[],1); %predictor 2: lick rate
        X = [ones(size(x1)) x1 x2 x1.*x2];
        
        all_stats_combined = [];
        all_fam_raw = [];
        all_nov_raw = [];
        all_fam_residual = [];
        all_nov_residual = [];
        
        %for each unit regardless of celltype
        for ii = 1:length(idx_all)
            f_rate = famCount(:,:,idx_all(ii)) ./ famOccup;
            n_rate = novCount(:,:,idx_all(ii)) ./ novOccup;
            f_size = size(f_rate); n_size = size(n_rate);
            combined_size = size([f_rate; n_rate]);
            
            y = reshape([f_rate;n_rate], [], 1); %combine familiar and novel data
            
            [b,bint,r,rint,stats] = regress(y,X);
            tempRes = reshape(r, combined_size);
            famRes = tempRes(1:f_size(1),:);
            novRes = tempRes(f_size(1)+1:end,:);
            
            %collect raw and residual firing data: trial-by-trial residuals per cell
            all_fam_raw(:,:,ii) = famCount(:,:,idx_all(ii)) ./ famOccup(:,:);
            all_fam_residual(:,:,ii) = famRes(:,:);
            
            all_nov_raw(:,:,ii) = novCount(:,:,idx_all(ii)) ./ novOccup(:,:);
            all_nov_residual(:,:,ii) = novRes(:,:);
            
            all_stats_combined = [all_stats_combined; stats];
        end
        
        %save output structure
        outputdir = fullfile(dirs.saveoutputstruct, 'trial-by-trial-residuals', 'Dist2RZ', datestr(now,'yymmdd'));
        if ~isfolder(outputdir)
            mkdir(outputdir)
        end
        outputdir = fullfile(outputdir,...
            [sessStr '_NovDay' num2str(novelDay) '_raw_vs_residuals_stim_singleModel.mat']);
        save(outputdir,'unitTypeLabel','unitID','allsess_unitType','position_binEdges',...
            'all_fam_raw', 'all_fam_residual',...
            'all_nov_raw','all_nov_residual','all_stats_combined',...
            'famSpeed','novSpeed',...
            'famLickrate','famLickrate','novLickrate','novLickrate',...
            'famTrialLabel','famTrialVoltage','novTrialLabel','novTrialVoltage',...
            'trialLabelDict');
        
        %% take an average raw and residuals per unit and combine all units.
        % the resulting values are size N x M with N = number of cells, M =
        % position bins. the string at the end of output denotes whther
        % all/success-only/misses-only trials were include for averaging
        % report NaN for averaged data if trial number < 5 for all units
        
        mean_fam_raw_hit = nan(length(unitID), length(position_binEdges) - 1);
        mean_fam_raw_miss = nan(length(unitID), length(position_binEdges) - 1);
        mean_nov_raw_hit = nan(length(unitID), length(position_binEdges) - 1);
        mean_nov_raw_miss = nan(length(unitID), length(position_binEdges) - 1);
        
        mean_fam_residual_hit = nan(length(unitID), length(position_binEdges) - 1);
        mean_fam_residual_miss = nan(length(unitID), length(position_binEdges) - 1);
        mean_nov_residual_hit = nan(length(unitID), length(position_binEdges) - 1);
        mean_nov_residual_miss = nan(length(unitID), length(position_binEdges) - 1);
        
        
        %fam data
        iTrial_rzCenter = find(famTrialLabel(:,2) == 1); %only use reward-centered trials
        iHitTrial = intersect(iTrial_rzCenter, find(famTrialLabel(:,3) == 1));
        iMissTrial = intersect(iTrial_rzCenter, find(famTrialLabel(:,3) == 0));
        i1stHalf = intersect(iTrial_rzCenter, find(famTrialLabel(:,1) == 1));
        i2ndHalf = intersect(iTrial_rzCenter, find(famTrialLabel(:,1) == 2));
        
        mean_fam_raw_all = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(iTrial_rzCenter,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_raw_1st = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(i1stHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_raw_2nd = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(i2ndHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_all = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(iTrial_rzCenter,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_1st = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(i1stHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_2nd = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(i2ndHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        if length(iHitTrial) >= 5
            mean_fam_raw_hit = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(iHitTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
            mean_fam_residual_hit = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(iHitTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
        end
        if length(iMissTrial) >= 5
            mean_fam_raw_miss = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(iMissTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
            mean_fam_residual_miss = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(iMissTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
        end
        
        %nov data
        iTrial_rzCenter = find(novTrialLabel(:,2) == 1); %only use reward-centered trials
        iHitTrial = intersect(iTrial_rzCenter, find(novTrialLabel(:,3) == 1));
        iMissTrial = intersect(iTrial_rzCenter, find(novTrialLabel(:,3) == 0));
        i1stHalf = intersect(iTrial_rzCenter, find(novTrialLabel(:,1) == 1));
        i2ndHalf = intersect(iTrial_rzCenter, find(novTrialLabel(:,1) == 2));
        
        mean_nov_raw_all = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(iTrial_rzCenter,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_raw_1st = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(i1stHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_raw_2nd = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(i2ndHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_all = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(iTrial_rzCenter,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_1st = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(i1stHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_2nd = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(i2ndHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        if length(iHitTrial) >= 5
            mean_nov_raw_hit = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(iHitTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
            mean_nov_residual_hit = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(iHitTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
        end
        if length(iMissTrial) >= 5
            mean_nov_raw_miss = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(iMissTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
            mean_nov_residual_miss = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(iMissTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
        end
        
        %combine all
        allsess_mean_fam_raw_all = [allsess_mean_fam_raw_all; mean_fam_raw_all];
        allsess_mean_fam_raw_1st = [allsess_mean_fam_raw_1st; mean_fam_raw_1st];
        allsess_mean_fam_raw_2nd = [allsess_mean_fam_raw_2nd; mean_fam_raw_2nd];
        allsess_mean_fam_raw_hit = [allsess_mean_fam_raw_hit; mean_fam_raw_hit];
        allsess_mean_fam_raw_miss = [allsess_mean_fam_raw_miss; mean_fam_raw_miss];
        allsess_mean_fam_residual_all = [allsess_mean_fam_residual_all; mean_fam_residual_all];
        allsess_mean_fam_residual_1st = [allsess_mean_fam_residual_1st; mean_fam_residual_1st];
        allsess_mean_fam_residual_2nd = [allsess_mean_fam_residual_2nd; mean_fam_residual_2nd];
        allsess_mean_fam_residual_hit = [allsess_mean_fam_residual_hit; mean_fam_residual_hit];
        allsess_mean_fam_residual_miss = [allsess_mean_fam_residual_miss; mean_fam_residual_miss];
        
        allsess_mean_nov_raw_all = [allsess_mean_nov_raw_all; mean_nov_raw_all];
        allsess_mean_nov_raw_1st = [allsess_mean_nov_raw_1st; mean_nov_raw_1st];
        allsess_mean_nov_raw_2nd = [allsess_mean_nov_raw_2nd; mean_nov_raw_2nd];
        allsess_mean_nov_raw_hit = [allsess_mean_nov_raw_hit; mean_nov_raw_hit];
        allsess_mean_nov_raw_miss = [allsess_mean_nov_raw_miss; mean_nov_raw_miss];
        allsess_mean_nov_residual_all = [allsess_mean_nov_residual_all; mean_nov_residual_all];
        allsess_mean_nov_residual_1st = [allsess_mean_nov_residual_1st; mean_nov_residual_1st];
        allsess_mean_nov_residual_2nd = [allsess_mean_nov_residual_2nd; mean_nov_residual_2nd];
        allsess_mean_nov_residual_hit = [allsess_mean_nov_residual_hit; mean_nov_residual_hit];
        allsess_mean_nov_residual_miss = [allsess_mean_nov_residual_miss; mean_nov_residual_miss];
        
        allsess_sessinfo = [allsess_sessinfo; repmat(currSess,length(unitID),1)];
        allsess_unitID = [allsess_unitID; unitID];
        allsess_unitType = [allsess_unitType; unitTypeLabel];
        allsess_stats = [allsess_stats; all_stats_combined];
    end
    outputdir = fullfile(dirs.saveoutputstruct, 'trial-by-trial-residuals', 'Dist2RZ', datestr(now,'yymmdd'));
    outputdir = fullfile(outputdir,'allsess_raw_vs_residuals_stim_singleModel.mat');
    save(outputdir, 'allsess_sessinfo','allsess_unitID','allsess_unitType','position_binEdges',...
        'allsess_mean_fam_raw_all','allsess_mean_fam_raw_hit','allsess_mean_fam_raw_miss',...
        'allsess_mean_fam_raw_1st','allsess_mean_fam_raw_2nd',...
        'allsess_mean_nov_raw_all','allsess_mean_nov_raw_hit','allsess_mean_nov_raw_miss',...
        'allsess_mean_nov_raw_1st','allsess_mean_nov_raw_2nd',...
        'allsess_mean_fam_residual_all','allsess_mean_fam_residual_hit','allsess_mean_fam_residual_miss',...
        'allsess_mean_fam_residual_1st','allsess_mean_fam_residual_2nd',...
        'allsess_mean_nov_residual_all','allsess_mean_nov_residual_hit','allsess_mean_nov_residual_miss',...
        'allsess_mean_nov_residual_1st','allsess_mean_nov_residual_2nd',...
        'allsess_stats');
end


%% create an output structure of lap-by-lap ratemap, separated by stim vs no stim
if regressLaps
    
    %running structure across multiple sessions
    allsess_sessinfo = []; allsess_unitID = []; allsess_unitType = []; allsess_stats = [];
    allsess_mean_fam_raw_all = [];
    allsess_mean_fam_raw_1st = []; %first session (first half of day)
    allsess_mean_fam_raw_2nd = []; %second session (seecond half of day)
    allsess_mean_nov_raw_all = [];
    allsess_mean_nov_raw_1st = [];
    allsess_mean_nov_raw_2nd = [];
    
    allsess_mean_fam_residual_all = [];
    allsess_mean_fam_residual_1st = [];
    allsess_mean_fam_residual_2nd = [];
    allsess_mean_nov_residual_all = [];
    allsess_mean_nov_residual_1st = [];
    allsess_mean_nov_residual_2nd = [];
    
    for iSess = 1:length(uniqSess)
        currSess = uniqSess(iSess, :);
        animal = currSess(1);
        if animal == 4
            params.iden = 'X';
        else
            params.ide = 'N';
        end
        recDay = currSess(2);
        novelDay = currSess(3);
        sessStr = [params.iden num2str(animal) '_' num2str(recDay)];
        sessRegStr = [params.iden num2str(animal) '_' num2str(recDay) '_' params.brainReg{1}];
        
        %load binned firing rate maps with behavioral info
        fdir = fullfile(dirs.data2load, 'singlesess_ratemap_laps');        
        lapfname = getlatestfile_with_string(fdir,...
            [sessStr '_createTrialByTrialRateMaps_Day ' num2str(novelDay)]);        
        lapMap = load(fullfile(fdir,lapfname));
        lapMap = lapMap.outmap;
        unitID = lapMap.unitID;        
        
        %% get cell type info
        iCellInMetrics = find(...
            arrayfun( @(x) strcmp(cell_metrics.sessionName{x}, sessRegStr),...
            1:length(cell_metrics.sessionName)));
        
        %check same unitIDs between outmap and cell_metrics structure, as
        %expected. 'idx_all' is 1:N within a single session, not the entire
        %cluster_metrics indices
        if sum(~ismember(unitID, cell_metrics.cluID(iCellInMetrics))) == 0
            idx_all = cell_metrics.cellID(iCellInMetrics(...
                ~ismember(iCellInMetrics, cell_metrics.tags.Bad))); %exclude units tagged as 'Bad'
            unitID = unitID(idx_all);
            unitTypeLabel = cell_metrics.putativeCellType(...
                iCellInMetrics(~ismember(iCellInMetrics, cell_metrics.tags.Bad)))';
        end
        
        %% collect behavioral data (smoothed speed & smoothed lick rate)
        position_binEdges = 0:5:360;
        
        famFiles    = lapMap.sessionInfo(lapMap.sessionInfo(:,6) == 1, 3);
        famCount    = lapMap.bin5.fam.rawCount;
        famOccup    = lapMap.bin5.fam.rawOccup;
        famSpeed    = lapMap.bin5.fam.smoothSpeed;
        famLickrate = lapMap.bin5.fam.smoothLickrate;
        
        novFiles    = lapMap.sessionInfo(lapMap.sessionInfo(:,6) ~= 1, 3);
        novCount    = lapMap.bin5.nov.rawCount;
        novOccup    = lapMap.bin5.nov.rawOccup;
        novSpeed    = lapMap.bin5.nov.smoothSpeed;
        novLickrate = lapMap.bin5.nov.smoothLickrate;
        
        %collect trial type and stim condition
        famTrialLabel = lapMap.bin5.fam.labels;
        famTrialVoltage = lapMap.bin5.fam.stimVoltage;
        novTrialLabel = lapMap.bin5.nov.labels;
        novTrialVoltage = lapMap.bin5.nov.stimVoltage;
        
        %organize X for fam and nov behavior that is consistent for all cells
        x1 = reshape([famSpeed; novSpeed],[],1); %predictor 1: speed
        x2 = reshape([famLickrate; novLickrate],[],1); %predictor 2: lick rate
        X = [ones(size(x1)) x1 x2 x1.*x2];
        
        all_stats_combined = [];
        all_fam_raw = [];
        all_nov_raw = [];
        all_fam_residual = [];
        all_nov_residual = [];
        
        %for each unit
        for ii = 1:length(idx_all)
            f_rate = famCount(:,:,idx_all(ii)) ./ famOccup(:,:,idx_all(ii));
            n_rate = novCount(:,:,idx_all(ii)) ./ novOccup(:,:,idx_all(ii));
            f_size = size(f_rate); n_size = size(n_rate);
            combined_size = size([f_rate; n_rate]);
            
            y = reshape([f_rate;n_rate], [], 1); %combine familiar and novel data
            
            [b,bint,r,rint,stats] = regress(y,X);
            tempRes = reshape(r, combined_size);
            famRes = tempRes(1:f_size(1),:);
            novRes = tempRes(f_size(1)+1:end,:);
            
            %collect raw and residual firing data: trial-by-trial residuals per cell
            all_fam_raw(:,:,ii) = famCount(:,:,idx_all(ii)) ./ famOccup(:,:,idx_all(ii));
            all_fam_residual(:,:,ii) = famRes(:,:);
            
            all_nov_raw(:,:,ii) = novCount(:,:,idx_all(ii)) ./ novOccup(:,:,idx_all(ii));
            all_nov_residual(:,:,ii) = novRes(:,:);
            
            all_stats_combined = [all_stats_combined; stats];
            
            % uncomment the following to plot data and model (related
            % with ExtendedDataFigure02_B example)
            if plotScatter_laps
                t = tiledlayout(1,1);
                scatter3(x1,x2,y,'blue')
                hold on
                x1fit = min(x1):max(x1);
                x2fit = min(x2):max(x2);
                [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
                YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
                mesh(X1FIT,X2FIT,YFIT)
                xlabel('Speed')
                ylabel('Lick rate (licks/s)')
                zlabel('Firing rate (spikes/s)')
                legend(['R-sq = ' num2str(stats(1))],'Location','best')
                legend boxoff
                hold off
                title(t, ['Familiar and Novel Combined - Interneuron ID' num2str(intidx(ii)) ' from Rec ' num2str(recInfo(ss))])

                figdir = fullfile(maindir, 'Demo_Figures', 'AcrossDays_MultipleLinearRegression_Laps_scatter_1Model');
                if ~isdir(figdir)
                    mkdir(figdir)
                end
                figname = fullfile(figdir, [params.iden num2str(animals(an), '%02d') '_' num2str(recInfo(ss)) '_Int_' num2str(ii, '%03d')]);
                print(gcf,figname,'-dpng')
                clf
            end
        end
        
        %% save output structure
        % individual sessions separately
        outputdir = fullfile(dirs.saveoutputstruct, 'trial-by-trial-residuals', 'Laps', datestr(now,'yymmdd'));
        if ~isfolder(outputdir)
            mkdir(outputdir)
        end
        outputdir = fullfile(outputdir,...
            [sessStr '_NovDay' num2str(novelDay) '_raw_vs_residuals_stim_singleModel.mat']);
        save(outputdir, 'unitTypeLabel','unitID','position_binEdges',...
            'all_fam_raw', 'all_fam_residual',...
            'all_nov_raw','all_nov_residual','all_stats_combined',...
            'famSpeed','novSpeed','famLickrate','novLickrate',...
            'famTrialLabel','famTrialVoltage','novTrialLabel','novTrialVoltage');
        
        
        %% take an average raw and residuals per unit and combine all units.
        % the resulting values are size N x M with N = number of cells, M =
        % position bins. for lap data, does not differentiate hits vs
        % misses and use all completed laps around the track as trials
        
        %fam data
        mean_fam_raw_all = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(:,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_raw_1st = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(famTrialLabel(:,1)==1,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_raw_2nd = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(famTrialLabel(:,1)==2,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_all = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(:,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_1st = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(famTrialLabel(:,1)==1,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_2nd = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(famTrialLabel(:,1)==2,:,x)),1:length(unitID), 'UniformOutput',false)');
        %nov data
        mean_nov_raw_all = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(:,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_raw_1st = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(novTrialLabel(:,1)==1,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_raw_2nd = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(novTrialLabel(:,1)==2,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_all = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(:,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_1st = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(novTrialLabel(:,1)==1,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_2nd = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(novTrialLabel(:,1)==2,:,x)),1:length(unitID), 'UniformOutput',false)');
        
        
        %combine all
        allsess_mean_fam_raw_all = [allsess_mean_fam_raw_all; mean_fam_raw_all];
        allsess_mean_fam_raw_1st = [allsess_mean_fam_raw_1st; mean_fam_raw_1st];
        allsess_mean_fam_raw_2nd = [allsess_mean_fam_raw_2nd; mean_fam_raw_2nd];
        allsess_mean_fam_residual_all = [allsess_mean_fam_residual_all; mean_fam_residual_all];
        allsess_mean_fam_residual_1st = [allsess_mean_fam_residual_1st; mean_fam_residual_1st];
        allsess_mean_fam_residual_2nd = [allsess_mean_fam_residual_2nd; mean_fam_residual_2nd];
        
        allsess_mean_nov_raw_all = [allsess_mean_nov_raw_all; mean_nov_raw_all];
        allsess_mean_nov_raw_1st = [allsess_mean_nov_raw_1st; mean_nov_raw_1st];
        allsess_mean_nov_raw_2nd = [allsess_mean_nov_raw_2nd; mean_nov_raw_2nd];
        allsess_mean_nov_residual_all = [allsess_mean_nov_residual_all; mean_nov_residual_all];
        allsess_mean_nov_residual_1st = [allsess_mean_nov_residual_1st; mean_nov_residual_1st];
        allsess_mean_nov_residual_2nd = [allsess_mean_nov_residual_2nd; mean_nov_residual_2nd];

        
        allsess_sessinfo = [allsess_sessinfo; repmat(currSess,length(unitID),1)];
        allsess_unitID = [allsess_unitID; unitID];
        allsess_unitType = [allsess_unitType; unitTypeLabel];
        allsess_stats = [allsess_stats; all_stats_combined];
    end
    outputdir = fullfile(dirs.saveoutputstruct, 'trial-by-trial-residuals', 'Laps', datestr(now,'yymmdd'));
    outputdir = fullfile(outputdir,'allsess_raw_vs_residuals_stim_singleModel.mat');
    save(outputdir, 'allsess_sessinfo','allsess_unitID','allsess_unitType','position_binEdges',...
        'allsess_mean_fam_raw_all','allsess_mean_nov_raw_all',...
        'allsess_mean_fam_residual_all','allsess_mean_nov_residual_all',...
        'allsess_mean_fam_raw_1st','allsess_mean_fam_raw_2nd',...
        'allsess_mean_nov_raw_1st','allsess_mean_nov_raw_2nd',...
        'allsess_mean_fam_residual_1st','allsess_mean_fam_residual_2nd',...
        'allsess_mean_nov_residual_1st','allsess_mean_nov_residual_2nd','allsess_stats');
end

%% create an output structure of Time-based residual firing rates
if regressTime2RZ
    %running structure across multiple sessions
    allsess_sessinfo = []; allsess_unitID = []; allsess_unitType = []; allsess_stats = [];
    allsess_mean_fam_raw_all = []; allsess_mean_fam_raw_hit = []; allsess_mean_fam_raw_miss = [];
    allsess_mean_fam_raw_1st = []; allsess_mean_fam_raw_2nd = [];
    allsess_mean_nov_raw_all = []; allsess_mean_nov_raw_hit = []; allsess_mean_nov_raw_miss = [];
    allsess_mean_nov_raw_1st = []; allsess_mean_nov_raw_2nd = [];
    allsess_mean_fam_residual_all = []; allsess_mean_fam_residual_hit = []; allsess_mean_fam_residual_miss = [];
    allsess_mean_fam_residual_1st = []; allsess_mean_fam_residual_2nd = [];
    allsess_mean_nov_residual_all = []; allsess_mean_nov_residual_hit = []; allsess_mean_nov_residual_miss = [];
    allsess_mean_nov_residual_1st = []; allsess_mean_nov_residual_2nd = [];
    
    for iSess = 1:length(uniqSess)
        currSess = uniqSess(iSess, :);
        animal = currSess(1);
        recDay = currSess(2);
        novelDay = currSess(3);
        sessStr = [params.iden num2str(animal) '_' num2str(recDay)];
        sessRegStr = [params.iden num2str(animal) '_' num2str(recDay) '_' params.brainReg{1}];
        
        %load binned firing rate maps with behavioral info
        fdir = fullfile(dirs.data2load, 'singlesess_ratemap_time2RZ');        
        time2RZfname = getlatestfile_with_string(fdir,...
            [sessStr '_createTrialByTrialTime2RewardRateMaps_Day ' num2str(novelDay)]);        
        timeMap = load(fullfile(fdir,time2RZfname));
        timeMap = timeMap.outmap;
        unitID = timeMap.unitID;
        
        
        %% get cell type info
        iCellInMetrics = find(...
            arrayfun( @(x) strcmp(cell_metrics.sessionName{x}, sessRegStr),...
            1:length(cell_metrics.sessionName)));
        
        %check same unitIDs between outmap and cell_metrics structure, as
        %expected. 'idx_all' is 1:N within a single session, not the entire
        %cluster_metrics indices
        if sum(~ismember(unitID, cell_metrics.cluID(iCellInMetrics))) == 0
            idx_all = cell_metrics.cellID(iCellInMetrics(...
                ~ismember(iCellInMetrics, cell_metrics.tags.Bad))); %exclude units tagged as 'Bad'
            unitID = unitID(idx_all);
            unitTypeLabel = cell_metrics.putativeCellType(...
                iCellInMetrics(~ismember(iCellInMetrics, cell_metrics.tags.Bad)))';
        end
                
        
        %% collect behavioral data (smoothed speed & smoothed lick rate)
        time_binEdges = timeMap.bin100.fam.binEdges;
        
        famFiles    = timeMap.sessionInfo(timeMap.sessionInfo(:,6) == 1, 3);
        famRatemap    = timeMap.bin100.fam.ratemap;
        famOccup    = timeMap.bin100.fam.rawOccup;
        famSpeed    = timeMap.bin100.fam.smoothSpeed;
        famLickrate = timeMap.bin100.fam.smoothLickrate;
        
        novFiles    = timeMap.sessionInfo(timeMap.sessionInfo(:,6) ~= 1, 3);
        novRatemap    = timeMap.bin100.nov.ratemap;
        novOccup    = timeMap.bin100.nov.rawOccup;
        novSpeed    = timeMap.bin100.nov.smoothSpeed;
        novLickrate = timeMap.bin100.nov.smoothLickrate;
        
        %collect trial type and stim intensity info (in volts)
        famTrialLabel = timeMap.bin100.fam.labels;
        famTrialVoltage = timeMap.bin100.fam.stimVoltage;
        novTrialLabel = timeMap.bin100.nov.labels;
        novTrialVoltage = timeMap.bin100.nov.stimVoltage;
        trialLabelDict = timeMap.bin100.fam.labels_dict;
        
        %organize X for fam and nov behavior that is consistent for all cells
        x1 = reshape([famSpeed; novSpeed],[],1); %predictor 1: speed
        x2 = reshape([famLickrate; novLickrate],[],1); %predictor 2: lick rate
        X = [ones(size(x1)) x1 x2 x1.*x2];
        
        all_stats_combined = [];
        all_fam_raw = [];
        all_nov_raw = [];
        all_fam_residual = [];
        all_nov_residual = [];
        
        %for each unit regardless of celltype
        for ii = 1:length(idx_all)
            f_rate = famRatemap(:,:,idx_all(ii));
            n_rate = novRatemap(:,:,idx_all(ii));
            f_size = size(f_rate); n_size = size(n_rate);
            combined_size = size([f_rate; n_rate]);
            
            y = reshape([f_rate;n_rate], [], 1); %combine familiar and novel data
            
            [b,bint,r,rint,stats] = regress(y,X);
            tempRes = reshape(r, combined_size);
            famRes = tempRes(1:f_size(1),:);
            novRes = tempRes(f_size(1)+1:end,:);
            
            %collect raw and residual firing data: trial-by-trial residuals per cell
            all_fam_raw(:,:,ii) = f_rate;
            all_fam_residual(:,:,ii) = famRes(:,:);
            
            all_nov_raw(:,:,ii) = n_rate;
            all_nov_residual(:,:,ii) = novRes(:,:);
            
            all_stats_combined = [all_stats_combined; stats];
        end
        
        %save output structure
        outputdir = fullfile(dirs.saveoutputstruct, 'trial-by-trial-residuals', 'Time2RZ', datestr(now,'yymmdd'));
        if ~isfolder(outputdir)
            mkdir(outputdir)
        end
        outputdir = fullfile(outputdir,...
            [sessStr '_NovDay' num2str(novelDay) '_raw_vs_residuals_stim_singleModel.mat']);
        save(outputdir,'unitTypeLabel','unitID','time_binEdges',...
            'all_fam_raw', 'all_fam_residual',...
            'all_nov_raw','all_nov_residual','all_stats_combined',...
            'famSpeed','novSpeed',...
            'famLickrate','famLickrate','novLickrate','novLickrate',...
            'famTrialLabel','famTrialVoltage','novTrialLabel','novTrialVoltage',...
            'trialLabelDict');
        
        %% take an average raw and residuals per unit and combine all units.
        % the resulting values are size N x M with N = number of cells, M =
        % position bins. the string at the end of output denotes whther
        % all/success-only/misses-only trials were include for averaging
        % report NaN for averaged data if trial number < 5 for all units
        
        mean_fam_raw_hit = nan(length(unitID), length(time_binEdges) - 1);
        mean_fam_raw_miss = nan(length(unitID), length(time_binEdges) - 1);
        mean_nov_raw_hit = nan(length(unitID), length(time_binEdges) - 1);
        mean_nov_raw_miss = nan(length(unitID), length(time_binEdges) - 1);
        
        mean_fam_residual_hit = nan(length(unitID), length(time_binEdges) - 1);
        mean_fam_residual_miss = nan(length(unitID), length(time_binEdges) - 1);
        mean_nov_residual_hit = nan(length(unitID), length(time_binEdges) - 1);
        mean_nov_residual_miss = nan(length(unitID), length(time_binEdges) - 1);
        
        
        %fam data
        iTrial_rzCenter = find(famTrialLabel(:,2) == 1); %only use reward-centered trials
        iHitTrial = intersect(iTrial_rzCenter, find(famTrialLabel(:,3) == 1));
        iMissTrial = intersect(iTrial_rzCenter, find(famTrialLabel(:,3) == 0));
        i1stHalf = intersect(iTrial_rzCenter, find(famTrialLabel(:,1) == 1));
        i2ndHalf = intersect(iTrial_rzCenter, find(famTrialLabel(:,1) == 2));
        
        mean_fam_raw_all = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(iTrial_rzCenter,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_raw_1st = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(i1stHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_raw_2nd = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(i2ndHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_all = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(iTrial_rzCenter,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_1st = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(i1stHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_fam_residual_2nd = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(i2ndHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        if length(iHitTrial) >= 5
            mean_fam_raw_hit = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(iHitTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
            mean_fam_residual_hit = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(iHitTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
        end
        if length(iMissTrial) >= 5
            mean_fam_raw_miss = cell2mat(arrayfun( @(x) nanmean(all_fam_raw(iMissTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
            mean_fam_residual_miss = cell2mat(arrayfun( @(x) nanmean(all_fam_residual(iMissTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
        end
        
        %nov data
        iTrial_rzCenter = find(novTrialLabel(:,2) == 1); %only use reward-centered trials
        iHitTrial = intersect(iTrial_rzCenter, find(novTrialLabel(:,3) == 1));
        iMissTrial = intersect(iTrial_rzCenter, find(novTrialLabel(:,3) == 0));
        i1stHalf = intersect(iTrial_rzCenter, find(novTrialLabel(:,1) == 1));
        i2ndHalf = intersect(iTrial_rzCenter, find(novTrialLabel(:,1) == 2));
        
        mean_nov_raw_all = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(iTrial_rzCenter,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_raw_1st = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(i1stHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_raw_2nd = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(i2ndHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_all = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(iTrial_rzCenter,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_1st = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(i1stHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        mean_nov_residual_2nd = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(i2ndHalf,:,x)),1:length(unitID), 'UniformOutput',false)');
        if length(iHitTrial) >= 5
            mean_nov_raw_hit = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(iHitTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
            mean_nov_residual_hit = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(iHitTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
        end
        if length(iMissTrial) >= 5
            mean_nov_raw_miss = cell2mat(arrayfun( @(x) nanmean(all_nov_raw(iMissTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
            mean_nov_residual_miss = cell2mat(arrayfun( @(x) nanmean(all_nov_residual(iMissTrial,:,x)),1:length(unitID), 'UniformOutput',false)');
        end
        
        %combine all
        allsess_mean_fam_raw_all = [allsess_mean_fam_raw_all; mean_fam_raw_all];
        allsess_mean_fam_raw_1st = [allsess_mean_fam_raw_1st; mean_fam_raw_1st];
        allsess_mean_fam_raw_2nd = [allsess_mean_fam_raw_2nd; mean_fam_raw_2nd];
        allsess_mean_fam_raw_hit = [allsess_mean_fam_raw_hit; mean_fam_raw_hit];
        allsess_mean_fam_raw_miss = [allsess_mean_fam_raw_miss; mean_fam_raw_miss];
        allsess_mean_fam_residual_all = [allsess_mean_fam_residual_all; mean_fam_residual_all];
        allsess_mean_fam_residual_1st = [allsess_mean_fam_residual_1st; mean_fam_residual_1st];
        allsess_mean_fam_residual_2nd = [allsess_mean_fam_residual_2nd; mean_fam_residual_2nd];
        allsess_mean_fam_residual_hit = [allsess_mean_fam_residual_hit; mean_fam_residual_hit];
        allsess_mean_fam_residual_miss = [allsess_mean_fam_residual_miss; mean_fam_residual_miss];
        
        allsess_mean_nov_raw_all = [allsess_mean_nov_raw_all; mean_nov_raw_all];
        allsess_mean_nov_raw_1st = [allsess_mean_nov_raw_1st; mean_nov_raw_1st];
        allsess_mean_nov_raw_2nd = [allsess_mean_nov_raw_2nd; mean_nov_raw_2nd];
        allsess_mean_nov_raw_hit = [allsess_mean_nov_raw_hit; mean_nov_raw_hit];
        allsess_mean_nov_raw_miss = [allsess_mean_nov_raw_miss; mean_nov_raw_miss];
        allsess_mean_nov_residual_all = [allsess_mean_nov_residual_all; mean_nov_residual_all];
        allsess_mean_nov_residual_1st = [allsess_mean_nov_residual_1st; mean_nov_residual_1st];
        allsess_mean_nov_residual_2nd = [allsess_mean_nov_residual_2nd; mean_nov_residual_2nd];
        allsess_mean_nov_residual_hit = [allsess_mean_nov_residual_hit; mean_nov_residual_hit];
        allsess_mean_nov_residual_miss = [allsess_mean_nov_residual_miss; mean_nov_residual_miss];
        
        allsess_sessinfo = [allsess_sessinfo; repmat(currSess,length(unitID),1)];
        allsess_unitID = [allsess_unitID; unitID];
        allsess_unitType = [allsess_unitType; unitTypeLabel];
        allsess_stats = [allsess_stats; all_stats_combined];
    end
    outputdir = fullfile(dirs.saveoutputstruct, 'trial-by-trial-residuals', 'Time2RZ', datestr(now,'yymmdd'));
    outputdir = fullfile(outputdir,'allsess_raw_vs_residuals_stim_singleModel.mat');
    save(outputdir, 'allsess_sessinfo','allsess_unitID','allsess_unitType','time_binEdges',...
        'allsess_mean_fam_raw_all','allsess_mean_fam_raw_hit','allsess_mean_fam_raw_miss',...
        'allsess_mean_fam_raw_1st','allsess_mean_fam_raw_2nd',...
        'allsess_mean_nov_raw_all','allsess_mean_nov_raw_hit','allsess_mean_nov_raw_miss',...
        'allsess_mean_nov_raw_1st','allsess_mean_nov_raw_2nd',...
        'allsess_mean_fam_residual_all','allsess_mean_fam_residual_hit','allsess_mean_fam_residual_miss',...
        'allsess_mean_fam_residual_1st','allsess_mean_fam_residual_2nd',...
        'allsess_mean_nov_residual_all','allsess_mean_nov_residual_hit','allsess_mean_nov_residual_miss',...
        'allsess_mean_nov_residual_1st','allsess_mean_nov_residual_2nd',...
        'allsess_stats');
end



