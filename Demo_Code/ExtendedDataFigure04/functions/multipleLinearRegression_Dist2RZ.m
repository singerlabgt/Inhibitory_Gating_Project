function [allsess_unitType, allsess_mean_fam_residual_hit, allsess_sessinfo, position_binEdges] = multipleLinearRegression_Dist2RZ(SavePath, cell_metrics)
%% Multiple linear regression

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023\';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir);
dirs.saveoutputstruct = SavePath;
% %load cell type info
% load(fullfile(dirs.data2load, 'cell_metrics.mat'));

%flags to create distance- or time-based residual firing rate maps 
regressLaps = 1;    %regress out firing rate maps for the entire track from 0-360 degrees (distance-based; in degrees)
regressDist2RZ = 1; %regress out firing rate maps for the areas around the reward zone, 60 degrees before and after the RZ (distance-based; in 5 deg increments)

%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);
uniqSess = unique(allindex(:, [1:2, 7]), 'rows');

%% create an output structure of RZ/NRZ-centered ratemap across reward trials, separated by stim vs no stim
if regressDist2RZ
    %running structure across multiple sessions
    allsess_sessinfo = []; allsess_unitID = []; allsess_unitType = []; allsess_stats = [];
    allsess_mean_fam_raw_all = []; allsess_mean_fam_raw_hit = []; allsess_mean_fam_raw_miss = [];
    allsess_mean_fam_raw_1st = []; allsess_mean_fam_raw_2nd = [];
    allsess_mean_fam_residual_all = []; allsess_mean_fam_residual_hit = []; allsess_mean_fam_residual_miss = [];
    allsess_mean_fam_residual_1st = []; allsess_mean_fam_residual_2nd = [];
    
    for iSess = 1:length(uniqSess)
        currSess = uniqSess(iSess, :);
        animal = currSess(1);
        if animal == 4
            params.iden = 'X';
        else
            params.iden = 'N';
        end
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
        
        %collect trial type and stim intensity info (in volts)
        famTrialLabel = distMap.bin5.fam.labels;
        famTrialVoltage = distMap.bin5.fam.stimVoltage;
        trialLabelDict = distMap.bin5.fam.labels_dict;
        
        %organize X for fam and nov behavior that is consistent for all cells
        x1 = reshape([famSpeed],[],1); %predictor 1: speed
        x2 = reshape([famLickrate],[],1); %predictor 2: lick rate
        X = [ones(size(x1)) x1 x2 x1.*x2];
        
        all_stats_combined = [];
        all_fam_raw = [];
        all_fam_residual = [];
        
        %for each unit regardless of celltype
        for ii = 1:length(idx_all)
            f_rate = famCount(:,:,idx_all(ii)) ./ famOccup;
            f_size = size(f_rate); 
            combined_size = size(f_rate);
            y = reshape(f_rate, [], 1);    
            [b,bint,r,rint,stats] = regress(y,X);
            tempRes = reshape(r, combined_size);
            famRes = tempRes(1:f_size(1),:);
            
            %collect raw and residual firing data: trial-by-trial residuals per cell
            all_fam_raw(:,:,ii) = famCount(:,:,idx_all(ii)) ./ famOccup(:,:);
            all_fam_residual(:,:,ii) = famRes(:,:);
            all_stats_combined = [all_stats_combined; stats];
        end
        
        %% take an average raw and residuals per unit and combine all units.
        % the resulting values are size N x M with N = number of cells, M =
        % position bins. the string at the end of output denotes whther
        % all/success-only/misses-only trials were include for averaging
        % report NaN for averaged data if trial number < 5 for all units
        
        mean_fam_raw_hit = nan(length(unitID), length(position_binEdges) - 1);
        mean_fam_raw_miss = nan(length(unitID), length(position_binEdges) - 1);
        
        mean_fam_residual_hit = nan(length(unitID), length(position_binEdges) - 1);
        mean_fam_residual_miss = nan(length(unitID), length(position_binEdges) - 1);       
        
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
                
        allsess_sessinfo = [allsess_sessinfo; repmat(currSess,length(unitID),1)];
        allsess_unitID = [allsess_unitID; unitID];
        allsess_unitType = [allsess_unitType; unitTypeLabel];
        allsess_stats = [allsess_stats; all_stats_combined];
    end
end
