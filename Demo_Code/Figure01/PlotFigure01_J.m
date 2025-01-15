%% Jeong et al. 2023 MANUSCRIPT - FIGURE01_J
% NJeong 03/28/2023

clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));
celltypes = {'Narrow Interneuron','Pyramidal Cell'}; %CellExplorer's default naming conventions
celltypeNames = {'NS Int.','Pyr.'}; %names used in our manuscript
for ct = 1:length(celltypes) %get indices of units per specified cell type
    cellT{ct} = find(...
        strcmp(cell_metrics.putativeCellType, celltypes{ct})...
        & ~ismember(1:length(cell_metrics.cellID), cell_metrics.tags.Bad));
end

%get all indices
[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.animals);
allindex = allindex(animal_idx,:); %filter based on animals to include 
allindex = allindex(animal_idx,:); %filter based on animals to include 

[sessions, sessID] = unique(allindex(:,[1:2, 7]),'rows'); %session = [animalID, recording date, novelty day]

%set variables
q = 0.01; %alpha level set for False Discovery Rate (FDR)
bins_preRZ = 1:35; %-60 to 0 degrees, prior to RZ
position = mean(getBinEdges(-60:2:70),2);
position = position(bins_preRZ);

stats_observed_fam = []; %stats of observed data in the order of r-square, pval, coeefficient estimate
stats_observed_nov = [];
y_pred_fam = []; y_pred_nov = [];
y_resp_fam = []; y_resp_nov = [];
lbl_fam = []; lbl_nov = [];


for iSess = 1:length(sessions)
    animal = sessions(iSess, 1);
    recDay = sessions(iSess, 2);
    novelDay = sessions(iSess, 3);
    params.iden = alliden{sessID(iSess)};
    
    %% load ratemap 
    originalstructdir = fullfile(dirs.data2load, 'singlesess_ratemap_distance2RZ');
    
    origFile = getlatestfile_with_string(originalstructdir,...
        [params.iden num2str(animal) '_' num2str(recDay) '_createTrialByTrialDist2RewardRateMaps_Day ' num2str(novelDay)]);
    load(fullfile(originalstructdir,origFile)) %load 'outmap'
    
    %% fit generalized linear model to firing rate pre-RZ as a function of position
    numUnits = length(outmap.unitID);
    %fam data
    observed_ratemap = outmap.bin2.fam.rawCount ./ outmap.bin2.fam.rawOccup;
    observed_ratemap(outmap.bin2.fam.rawSpeed < params.speedTh) = nan;
    observed_ratemap = observed_ratemap(outmap.bin2.fam.labels(:,2) == 1, :, :);
    observed_ratemap = nanmean( observed_ratemap(:,:,:), 1);
    observed_ratemap = reshape(observed_ratemap,[],size(observed_ratemap,3));
    
    for iUnit = 1:numUnits
        y = observed_ratemap(:,iUnit);
        y = y(bins_preRZ);
        tbl = table(y, position);
        glm = fitglm(tbl,'y ~ 1 + position');
        
        stats_observed_fam = [stats_observed_fam; glm.Rsquared.Ordinary, glm.Coefficients.pValue(2), glm.Coefficients.Estimate(2)];
        y_pred_fam = [y_pred_fam, glm.predict];
        y_resp_fam = [y_resp_fam, y];
        lbl_fam = [ lbl_fam; sessions(iSess,:), outmap.unitID(iUnit)];
    end
    
    %novel data
    observed_ratemap = outmap.bin2.nov.rawCount ./ outmap.bin2.nov.rawOccup;
    if ~isempty(observed_ratemap)
        observed_ratemap(outmap.bin2.nov.rawSpeed < params.speedTh) = nan;
        observed_ratemap = observed_ratemap(outmap.bin2.nov.labels(:,2) == 1, :, :);
        observed_ratemap = nanmean( observed_ratemap(:,:,:), 1);
        observed_ratemap = reshape(observed_ratemap,[],size(observed_ratemap,3));
        
        for iUnit = 1:numUnits
            y = observed_ratemap(:,iUnit);
            y = y(bins_preRZ);
            tbl = table(y, position);
            glm = fitglm(tbl,'y ~ 1 + position');
            
            stats_observed_nov = [stats_observed_nov; glm.Rsquared.Ordinary, glm.Coefficients.pValue(2), glm.Coefficients.Estimate(2)];
            y_pred_nov = [y_pred_nov, glm.predict];
            y_resp_nov = [y_resp_nov, y];
            lbl_nov = [ lbl_nov; sessions(iSess,:), outmap.unitID(iUnit)];
        end
    end
end


%% create a label structure to identify correct session and cell info
lbl_fam(:,5) = nan; 
for ii = 1:length(lbl_fam)
    tempStr = [params.iden num2str(lbl_fam(ii,1)) '_' num2str(lbl_fam(ii,2)) '_CA3' ];
    idx = find(strcmp(cell_metrics.sessionName, tempStr) & lbl_fam(ii,4) == cell_metrics.cluID); 
    if ismember(idx, cellT{1})
        lbl_fam(ii,5) = 1; %narrow interneuron
    elseif ismember(idx, cellT{2})
        lbl_fam(ii,5) = 2; %pyramidal cell
    else
        lbl_fam(ii,5) = nan;
    end
end

lbl_nov(:,5) = nan; 
for ii = 1:length(lbl_nov)
    tempStr = [params.iden num2str(lbl_nov(ii,1)) '_' num2str(lbl_nov(ii,2)) '_CA3' ];
    idx = find(strcmp(cell_metrics.sessionName, tempStr) & lbl_nov(ii,4) == cell_metrics.cluID); 
    if ismember(idx, cellT{1})
        lbl_nov(ii,5) = 1; %narrow interneuron
    elseif ismember(idx, cellT{2})
        lbl_nov(ii,5) = 2; %pyramidal cell
    else
        lbl_nov(ii,5) = nan;
    end
end


%% determine neurons with significant ramp down/up activity 
[h, crit_p, adj_p]=fdr_bh(stats_observed_fam(:,2),q,'pdep','yes');
sigNeurons_fam = find(adj_p < q);

[h, crit_p, adj_p]=fdr_bh(stats_observed_nov(:,2),q,'pdep','yes');
sigNeurons_nov = find(adj_p < q);

sigNarrowInts_fam = intersect(sigNeurons_fam, find(lbl_fam(:,5) == 1));
sigNarrowInts_fam = intersect(sigNarrowInts_fam, find(stats_observed_fam(:,3) < 0)); %only include units with negative slope
sigNarrowInts_nov = intersect(sigNeurons_nov, find(lbl_nov(:,5) == 1));
sigNarrowInts_nov = intersect(sigNarrowInts_nov, find(stats_observed_nov(:,3) < 0)); %only include units with negative slope

sigPyrs_fam = intersect(sigNeurons_fam, find(lbl_fam(:,5) == 2));
sigPyrs_fam = intersect(sigPyrs_fam, find(stats_observed_fam(:,3) >= 0)); %only include units with positive slope
sigPyrs_nov = intersect(sigNeurons_nov, find(lbl_nov(:,5) == 2));
sigPyrs_nov = intersect(sigPyrs_nov, find(stats_observed_nov(:,3) >= 0)); %only include units with positive slope



%% plot proportion of narrow/wide interneurons with significant negative slope
%averaged across sessions, over days 1-3
%wild-type narrow interneurons fam + nov day 1-3, one data point per recDay

fig = figure('units','inch','position',[0 0 3.5 2]);
ax = arrayfun( @(x) subplot(1,length(celltypes),x,...
    'NextPlot','add','Box','off','XLim',[0.5 3.5],'XTick',1:3), 1:length(celltypes));

for ct = 1:length(celltypes)
    tempFam = intersect(sigNeurons_fam, find(lbl_fam(:,5) == ct));
    tempNov = intersect(sigNeurons_nov, find(lbl_nov(:,5) == ct));
    
    lmeFam = []; lmeFamlbl = [];
    lmeNov = []; lmeNovlbl = [];
    
    for ii = 1:3
        %fam 
        tempList = lbl_fam(...
            ismember(lbl_fam(:,1),params.WTmice)...
            & lbl_fam(:,3) == ii...
            & lbl_fam(:,5) == ct, 1:3);
        recDays = unique(tempList,'rows');
        nTotal = arrayfun( @(x) sum(tempList(:,1) == recDays(x,1) & tempList(:,2) == recDays(x,2)), 1:size(recDays,1));
        tempList = lbl_fam(...
            intersect(tempFam,...
            find(ismember(lbl_fam(:,1),params.WTmice)...
            & lbl_fam(:,3) == ii & lbl_fam(:,5) == ct)), 1:2);
        nSigUnits = arrayfun( @(x) sum(tempList(:,1) == recDays(x,1) & tempList(:,2) == recDays(x,2)), 1:size(recDays,1));
        propSig_fam{ii} = nSigUnits ./ nTotal;
        lmeFam = [lmeFam; [recDays,propSig_fam{ii}']]; lmeFamlbl = [lmeFamlbl; repmat({'Fam'},size(recDays,1),1)];
        %nov
        tempList = lbl_nov(...
            ismember(lbl_nov(:,1),params.WTmice)...
            & lbl_nov(:,3) == ii...
            & lbl_nov(:,5) == ct, 1:3);
        recDays = unique(tempList,'rows');
        nTotal = arrayfun( @(x) sum(tempList(:,1) == recDays(x,1) & tempList(:,2) == recDays(x,2)), 1:size(recDays,1));
        tempList = lbl_nov(...
            intersect(tempNov,...
            find(ismember(lbl_nov(:,1),params.WTmice)...
            & lbl_nov(:,3) == ii & lbl_nov(:,5) == ct)), 1:2);
        nSigUnits = arrayfun( @(x) sum(tempList(:,1) == recDays(x,1) & tempList(:,2) == recDays(x,2)), 1:size(recDays,1));
        propSig_nov{ii} = nSigUnits ./ nTotal;
    end
    
    %plotting errorbars  
    x = [1 2 3];
    vals = [arrayfun( @(x) mean(propSig_fam{x}), 1:3); arrayfun( @(x) mean(propSig_nov{x}), 1:3)];
    errs = [arrayfun( @(x) std(propSig_fam{x})./sqrt(length(propSig_fam{x})),1:3);...
        arrayfun( @(x) std(propSig_nov{x})./sqrt(length(propSig_nov{x})), 1:3)];
    errorbar(ax(ct), x, vals(1,:), errs(1,:), 'Color', params.colors_fam(end,:),'LineWidth',2)
    errorbar(ax(ct), x, vals(2,:), errs(2,:), 'Color', params.colors_nov(2,:),'LineWidth',2)    
    xlabel(ax(ct), 'Day')
    title(ax(ct), celltypeNames{ct})
    ylim([0, 0.6])
end
ylabel(ax(1), {'Proportion of units with'; 'significant activity change'})


%% save figure
makefigurepretty(gcf)
figname = 'Figure01_J';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')