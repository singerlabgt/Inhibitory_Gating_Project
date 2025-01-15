%% Jeong et al. 2023 MANUSCRIPT - FIGURE 04_F,G,H,I
% NJeong 03/23/2023

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir);

%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);
sessions = unique(allindex(:,1:2),'rows');

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Interneuron','Pyramidal Cell'}; %names used in our manuscript

%load NeuronFile structure
load(fullfile(dirs.data2load, 'NeuronFile.mat'));

%load RipData structure
load(fullfile(dirs.data2load, 'RipData_250ms.mat'));

%set paramters
params.pfBinsize = 5; %in deg
params.timeAroundMid = 0.125; %in s
params.decodingBin_s = params.timeAroundMid * 2; %use one timebin per ripple (i.e. do not break into multiple timebins)
params.decodingBin_samp = params.decodingBin_s .* params.samprate;
probThreshold = 0; %use time bins with spatial probability above median spatial probability
binEdges = -30:params.pfBinsize :70;
midBins = mean(getBinEdges(binEdges),2);

%% gather animal's current positions for each ripple event 
for iSess = 1:length(RipData)
    if ~isempty(RipData(iSess).Distance2RZ)
        
        for iRipple = 1:length(RipData(iSess).Distance2RZ)
%             RipData(iSess).Distance2RZ(iRipple).regressB = RipData(iSess).Distance2RZ(iRipple).regress.B(1); 
            currFile = RipData(iSess).Distance2RZ(iRipple).index(3); %find out which file that ripple came from
            env = RipData(iSess).Distance2RZ(iRipple).index(6); %identify current environment for that reipple 
            if env == 1
                dynFieldName = 'FamRewardTrials';
            else
                dynFieldName = 'NovRewardTrials';
            end
            
            %get index of ripple mid-point from data structure 
            iRipMid = RipData(iSess).Distance2RZ(iRipple).rippleTime(1) + params.timeAroundMid.*params.samprate; %add the time subtracted to get the ripple mid-point               
           
            %check if the current ripple was during VR session 
            iTrialInFile = find(currFile == vertcat(NeuronFile(iSess).(dynFieldName).currFile)); 
            iStartEnd = [NeuronFile(iSess).(dynFieldName)(iTrialInFile(1)).EphysIdx(1), NeuronFile(iSess).(dynFieldName)(iTrialInFile(end)).EphysIdx(end)];
            if iRipMid < iStartEnd(1) || iRipMid > iStartEnd(2)
                RipData(iSess).Distance2RZ(iRipple).animalPos_360 = nan; %mark as NaN for invalid ripples outside VR session 
                RipData(iSess).Distance2RZ(iRipple).animalPos_100 = nan;
            else            
                %find ephys index for the ripple event in current file
                ephysIdx = arrayfun( @(x) lookup2(iRipMid, NeuronFile(iSess).(dynFieldName)(x).EphysIdx), 1:length(NeuronFile(iSess).(dynFieldName)))';
                
                %find virmen index equivalent to ephys index
                
                virmenIdx = max(find(ephysIdx ~= 1 & currFile == vertcat(NeuronFile(iSess).(dynFieldName).currFile) ));
                
                %save animal position when ripple occurred, both in terms of
                %track position (0-360) or relative to RZ (-30 to 70)
                RipData(iSess).Distance2RZ(iRipple).animalPos_360 = NeuronFile(iSess).(dynFieldName)(virmenIdx).Pos360(ephysIdx(virmenIdx));
                RipData(iSess).Distance2RZ(iRipple).animalPos_100 = NeuronFile(iSess).(dynFieldName)(virmenIdx).Pos_30to70(ephysIdx(virmenIdx));
            end
        end
    end
end

animalPos360 = []; animalPos100 = []; ripDur_s = [];
for iSess = 1:length(sessions)
    if ~isempty(RipData(iSess).Distance2RZ)
        animalPos360 = [animalPos360; arrayfun( @(x) RipData(iSess).Distance2RZ(x).animalPos_360, 1:length(RipData(iSess).Distance2RZ))'];
        animalPos100 = [animalPos100; arrayfun( @(x) RipData(iSess).Distance2RZ(x).animalPos_100, 1:length(RipData(iSess).Distance2RZ))'];
        ripDur_s = [ripDur_s; arrayfun( @(x) RipData(iSess).Distance2RZ(x).rippleDuration,...
            1:length(RipData(iSess).Distance2RZ))'];
    end
end
ripDur_s(ripDur_s > 2) = nan;


%% collect all session information per ripple (one row of session info per ripple)
uniqMice = unique(sessions(:,1));
for iMouse = 1:length(uniqMice)
    index_rips = []; 
    index = find(ismember(sessions(:,1), uniqMice(iMouse)));
    for ii = 1:length(index)
        if ~isempty(RipData(index(ii)).Distance2RZ)
            index_rips = [index_rips; [cell2mat({RipData(index(ii)).Distance2RZ.index}'), [RipData(index(ii)).Distance2RZ.animalPos_100]'] ];
        end
    end
    pos(iMouse).rips = index_rips;
end
index_rips = vertcat(pos(:).rips);


%% calc distribution of decoded position during ripples - animal's current position 
%if decoded - actual postion > 0, forward replay (RZ in front of animal)
%and if < 0, backward replay (reactivation of RZ behind the animal), and if
%equal to zero, then currrent position reactivation 

mean_spatialProb_per_ripple = [];
for iSess = 1:length(sessions)
    if ~isempty(RipData(iSess).Distance2RZ)
        fam_decoded_minus_actual = []; nov_decoded_minus_actual = [];
        fam_decoded = []; nov_decoded = [];        
        allProb = cell2mat( arrayfun( @(x) ...
            nanmean(RipData(iSess).Distance2RZ(x).spatialProb,2), ...
            1:length(RipData(iSess).Distance2RZ), 'UniformOutput',false));
        sessThresh = prctile(max(allProb,[],1), probThreshold); %incldue only above threshold set per session
%         sessThresh = 0; %incldue all ripples for now
        
        %find decoded position bin of the highest spatial probability (i.e.
        %one value per ripple event) 
        for iRipple = 1:length(RipData(iSess).Distance2RZ)
            %collapse all time bins into one probability density plot 
            spaProb = RipData(iSess).Distance2RZ(iRipple).spatialProb; 
            mean_spatialProb_per_ripple = [mean_spatialProb_per_ripple; spaProb'];
            
            %find peak spatial probability per time bin
            [m, idx] = max(spaProb);
            
            %find position bins of peak spatial prob above session threshold
            rPos_incl = idx(m > sessThresh)';
            
            %collect data based on environment
            if RipData(iSess).Distance2RZ(iRipple).index(6) == 1
                fam_decoded = [fam_decoded; midBins(rPos_incl)];
                fam_decoded_minus_actual = [fam_decoded_minus_actual; midBins(rPos_incl) - RipData(iSess).Distance2RZ(iRipple).animalPos_100];
            else
                nov_decoded = [nov_decoded; midBins(rPos_incl)];
                nov_decoded_minus_actual = [nov_decoded_minus_actual; midBins(rPos_incl) - RipData(iSess).Distance2RZ(iRipple).animalPos_100];
            end
        end
        
        rPosData(iSess).fam_decoded = fam_decoded;
        rPosData(iSess).fam_decoded_minus_actual = fam_decoded_minus_actual;
        rPosData(iSess).nov_decoded = nov_decoded;
        rPosData(iSess).nov_decoded_minus_actual = nov_decoded_minus_actual;
    end
end

% plot spatial probability per condition
pvMice2incl = params.goalshamMice; 
maxProb = max(mean_spatialProb_per_ripple,[],2); 
Th = prctile(maxProb,probThreshold); 
ripDurTh = 0; %not applying ripple duration thereshold 

% define criteria to pick out the right ripples 
aboveTh = maxProb > Th;
famSess = index_rips(:,6)==1;
novSess = index_rips(:,6)~=1;
wtMice = ismember(index_rips(:,1),params.WTmice); 
pvMice = ismember(index_rips(:,1), pvMice2incl); 
goalRipple = abs(index_rips(:,10)) <= 10; %include ripples that occurred within 10 deg of RZ start (AZ+RZ only)
noStim = index_rips(:,8)==0;
goalStim = index_rips(:,9)==1;
shamStim = index_rips(:,9)==2; 
longRip = ripDur_s > ripDurTh; 

%find # goal ripples with peak spaital prob at either near goal (-10 to 20 deg),
%or far from goal (40 to 70 deg), followed by the total
%number of ripples for each of 4 conditions (WT fam, WT nov, goal stim, and
%sham stim)

%% Fig. 5F,G,H
fig5 = figure('units','inch','position',[0 0 5 3]);
t = tiledlayout(2,3,'TileSpacing','compact','Units','inches','OuterPosition',[0 0 5 3]);

colormap(gray);
h(1) = nexttile;
temp = mean_spatialProb_per_ripple(pvMice & novSess & aboveTh & longRip & goalStim & goalRipple,:);
tempIdx = index_rips(pvMice & novSess & aboveTh & longRip & goalStim & goalRipple,:); 
[~,b] = max(temp,[],2); [~,b2] = sort(b); 
uniqSess = unique(tempIdx(:,1:2),'rows'); 
for ua = 1:length(uniqSess)
    nRips.goalstimNov(ua,:) = [...
        sum(ismember(b(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2)), 5:10)),...
        sum(ismember(b(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2)), 15:20)),...
        sum(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2))];
end
imagesc(temp(b2,:)./max(temp(b2,:),[],2)); box off; 
xticks([1,7,20]); xticklabels(binEdges([1,7,21]))
title('GS in Nov')

h(2) = nexttile;
temp = mean_spatialProb_per_ripple(pvMice & novSess & aboveTh & longRip & shamStim & goalRipple,:);
tempIdx = index_rips(pvMice & novSess & aboveTh & longRip & shamStim & goalRipple,:); 
[~,b] = max(temp,[],2); [~,b2] = sort(b); 
uniqSess = unique(tempIdx(:,1:2),'rows'); 
for ua = 1:length(uniqSess)
    nRips.shamstimNov(ua,:) = [...
        sum(ismember(b(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2)), 5:10)),...
        sum(ismember(b(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2)), 15:20)),...
        sum(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2))];
end
imagesc(temp(b2,:)./max(temp(b2,:),[],2)); box off; 
xticks([1,7,20]); xticklabels(binEdges([1,7,21]))
title('SS in Nov')

h(3) = nexttile; hold on;
ripNames = {'goalstimNov','shamstimNov'};
meandata = nan(length(ripNames),2);
semdata = nan(length(ripNames),2);
colors = [params.colors_goalstim(3,:); params.colors_shamstim(3,:)];
clear invpoints;
for ii = 1:length(ripNames)
    temp = nRips.(ripNames{ii}); 
    temp = temp(:,1:2) ./ temp(:,3); 
    disp(ripNames(ii))
    disp(num2str(size(temp)))
    meandata(ii,:) = nanmean(temp); 
    semdata(ii,:) = nanstd(temp) ./ sqrt(size(temp,1));
    invpoints{ii,:} = temp;     
end
% arrayfun( @(x) scatter(h(3), ones(size(invpoints{1},1),1) .*x-0.1, invpoints{1}(:,x),20,...
%     'MarkerFaceColor', colors(1,:), 'MarkerFaceAlpha', 0.5,...
%     'MarkerEdgeColor', colors(1,:), 'MarkerEdgeAlpha', 0.5,...
%     'jitter', 'on', 'jitterAmount', 0.01),1:2);
% arrayfun( @(x) scatter(h(3), ones(size(invpoints{2},1),1) .*x+0.1, invpoints{2}(:,x),20,...
%     'MarkerFaceColor', colors(2,:), 'MarkerFaceAlpha', 0.5,...
%     'MarkerEdgeColor', colors(2,:), 'MarkerEdgeAlpha', 0.5,...
%     'jitter', 'on', 'jitterAmount', 0.01),1:2);
arrayfun( @(x)...
    errorbar(h(3),meandata(x,:),semdata(x,:),'LineWidth',2,...
    'Color',colors(x,:)),1:length(ripNames)); 
box off; hold off;
xlim([0.5,2.5]); xticks(1:2); xticklabels({'Near-goal','Far-goal'});
ylabel('Proportion of ripples')

h(4) = nexttile;
temp = mean_spatialProb_per_ripple(wtMice & famSess & aboveTh & longRip & goalRipple,:);
tempIdx = index_rips(wtMice & famSess & aboveTh & longRip & goalRipple,:); 
[~,b] = max(temp,[],2); [~,b2] = sort(b); 
uniqSess = unique(tempIdx(:,1:2),'rows'); 
for ua = 1:length(uniqSess)
    nRips.wtFam(ua,:) = [...
        sum(ismember(b(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2)), 5:10)),...
        sum(ismember(b(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2)), 15:20)),...
        sum(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2))];
end
imagesc(temp(b2,:)./max(temp(b2,:),[],2)); box off; 
xticks([1,7,20]); xticklabels(binEdges([1,7,21]))
title('WT/Fam')

h(5) = nexttile;
temp = mean_spatialProb_per_ripple(wtMice & novSess & aboveTh & longRip & goalRipple,:);
tempIdx = index_rips(wtMice & novSess & aboveTh & longRip & goalRipple,:); 
[~,b] = max(temp,[],2); [~,b2] = sort(b); 
uniqSess = unique(tempIdx(:,1:2),'rows'); 
for ua = 1:length(uniqSess)
    nRips.wtNov(ua,:) = [...
        sum(ismember(b(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2)), 5:10)),...
        sum(ismember(b(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2)), 15:20)),...
        sum(tempIdx(:,1)==uniqSess(ua,1) & tempIdx(:,2)==uniqSess(ua,2))];
end
imagesc(temp(b2,:)./max(temp(b2,:),[],2)); box off; 
xticks([1,7,20]); xticklabels(binEdges([1,7,21]))
title('WT/Nov')

h(6) = nexttile; hold on; 
ripNames = {'wtFam','wtNov'};
meandata = nan(length(ripNames),2);
semdata = nan(length(ripNames),2);
colors = [params.colors_fam(3,:); params.colors_nov(2,:)];
clear invpoints
for ii = 1:length(ripNames)
    temp = nRips.(ripNames{ii}); 
    temp = temp(:,1:2) ./ temp(:,3); 
    disp(ripNames(ii))
    disp(num2str(size(temp)))
    meandata(ii,:) = nanmean(temp); 
    semdata(ii,:) = nanstd(temp) ./ sqrt(size(temp,1));
    invpoints{ii,:} = temp;
end

arrayfun( @(x)...
    errorbar(h(6),meandata(x,:),semdata(x,:),'LineWidth',2,...
    'Color',colors(x,:), ...
    'DisplayName',ripNames{x}),1:length(ripNames)); 
arrayfun( @(x)...
    errorbar(h(6),meandata(x,:),semdata(x,:),'LineWidth',2,...
    'Color',colors(x,:)),1:length(ripNames)); 
box off; hold off;
xlim([0.5,2.5]); xticks(1:2); xticklabels({'Near-goal','Far-goal'});
ylabel('Proportion of ripples')

linkaxes(h(1:2),'xy')
linkaxes(h(4:5),'xy')
xlabel(t,'Decoded distance to RZ during ripple (deg)')
ylabel(t,'Ripple count')
cb = colorbar;
cb.Layout.Tile = 'east';

%% save figure
makefigurepretty(gcf)
figname = 'Figure04_FGH';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')