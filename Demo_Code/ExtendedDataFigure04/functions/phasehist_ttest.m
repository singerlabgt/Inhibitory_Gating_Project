% The t-statistical approach:
% 1) For each potential subtypes (AAC, PVBC, CCK, bistratified), define the preferred phase as mentioned in the paper. 
% 2) Do one tail two-sample t-test, to compare the spike counts in the preferred/unpreferred phase. The length of each category corresponds to the preferred/unpreferred phase bin number (binwidth=18). Record p-values and t scores for each cell type t-test.
% 3) Do FDR estimation given the 4 cell types, and record adjusted p-values (q-values).
% 4) Thresholding: q-values<0.05. 
% If multiple cell types have q-values passing the threshold, then compare the first 2 cell types with largest t-scores:
% a. if t-score difference is larger than 0.5, then assign the cell type with the largest t-score to that neuron.
% b. if t-score difference <= 0.5, then assign both 2 cell types with the largest t-score to that neuron.
%% TODO: 1. SET SPIKE_THR (DONE); 2. SET TO SMALLER BINWIDTH (ABORTED)

% only compare the known PV+ HPC interneurons -- PVBC, AAC, Bistratified.
subtypes = ["pvbc", "aac", "bistratified", "cck"];
pvbc_mean = 271;
aac_mean = 185;
bistratified_mean = 79;
cck_mean = 155; 
pyr_mean = 20;
if universal_tol == 1
    pvbc_tol = shift_thr;
    aac_tol = shift_thr;
    bistratified_tol = shift_thr;
    cck_tol = 81*shift_factor;
else
    pvbc_tol = 68*scale_factor;
    aac_tol = 55*scale_factor;
    bistratified_tol = 92*scale_factor;
    cck_tol = 81*scale_factor;
end


for i = 1:length(NeuronFile)
    PYR_notbad = intersect(find(ismember(NeuronFile(i).TSIDtotal, NeuronFile(i).TSIDnotbad)), find(NeuronFile(i).CellType == 1));
    if pyr_curation == 1 & ~isempty(PYR_notbad)
        binnedphase_pyr = [];
        for iCell = 1:length(NeuronFile(i).CellType)
            if NeuronFile(i).SpikeCountInTheta_PYR{iCell} >= spike_thr
                binnedphase_pyr = [binnedphase_pyr; (NeuronFile(i).BinSize1ThetaPhase_PYR{iCell}) ./ max(NeuronFile(i).BinSize1ThetaPhase_PYR{iCell})];
            end
        end
        % determine the preferred phase of pyr cells and calculate the
        % CA1-CA3 shift (preferred phase for CA1 pyr: 20 deg)
        spike_prob_pyr = [];
        for iBin = 1:length(x_bins)
            if iBin < binwidth/2 
                spike_prob_pyr = [spike_prob_pyr, mean([binnedphase_pyr(1 : binwidth/2 + iBin), binnedphase_pyr(end - (binwidth/2 - iBin) : end)])];
            elseif iBin > length(x_bins) - binwidth/2 
                spike_prob_pyr = [spike_prob_pyr, mean([binnedphase_pyr(1 : iBin + binwidth/2 - length(x_bins)), binnedphase_pyr(iBin - binwidth/2: end)])];
            else
                spike_prob_pyr = [spike_prob_pyr, mean(iBin - binwidth/2 : iBin + binwidth/2)];
            end
        end
        [M, I] = max(spike_prob_pyr); % index of 20 bins
        % pyr_shift = I - 2; % minus 2: this is because pyramidal cells
        % fire at 20deg which corresponds to the 2nd x_bin. If positive,
        % then actual phase is later than ref, should set to earlier
        pyr_shift = I - pyr_mean;
    else
        pyr_shift = 0;
    end
    NS_notbad = intersect(find(ismember(NeuronFile(i).TSIDtotal, NeuronFile(i).TSIDnotbad)), find(NeuronFile(i).CellType == 2));
    for iCell = 1:length(NeuronFile(i).CellType)
        NeuronFile(i).ThetaBasedCellType{iCell} = [];
        NeuronFile(i).BinnedThetaPhase{iCell} = [];
        NeuronFile(i).SpikeCountInTheta{iCell} = [];
        if ismember(iCell, NS_notbad)
            data = NeuronFile(i).ThetaPhase{iCell};
            data = toDegrees('radians', data) + 180 ;
            data1 = zeros(size(data));
            for iphase = 1:length(data1)
                tophase = data(iphase) - pyr_shift;
                if tophase > 360
                    tophase = tophase - 360;
                end
                if tophase < 0
                    tophase = 360 - tophase + 1;
                end
                data1(iphase) = tophase;
            end
            data = data1;
            
            h = histogram(data, 'binWidth', binwidth, 'BinLimits', [0 360]); 
            binnedphase = h.BinCounts;
            % binnedphase0 = h.BinCounts;
            % binnedphase = zeros(size(binnedphase0));
            % for iBin = 1:length(x_bins)
            %     toBin = iBin - pyr_shift;
            %     if toBin > length(x_bins)
            %         toBin = toBin - length(x_bins);
            %     end
            %     if toBin < 0
            %         toBin = length(x_bins) - toBin + 1;
            %     end
            %     binnedphase(toBin) = binnedphase0(iBin);
            % end
            NeuronFile(i).BinnedThetaPhase{iCell} = binnedphase;
            NeuronFile(i).SpikeCountInTheta{iCell} = sum(binnedphase);
            if sum(binnedphase) <= spike_thr
                NeuronFile(i).ThetaBasedCellType{iCell} = ["ungrouped_ns_pv"];
                continue
            end
            % x_bins = h.BinEdges(1:end-1)' + h.BinWidth/2;
            pvals = [];
            tstats = [];
            for iType = 1:length(subtypes)
                subtype = char(subtypes(iType));
                subtype_mean = eval([subtype '_mean']);
                subtype_tol = eval([subtype '_tol']);
                lowerbound = subtype_mean - subtype_tol;
                upperbound = subtype_mean + subtype_tol;
                incircleperiod = [max(lowerbound, 0), min(upperbound, 360)];
                lowerboundexceedperiod = [];
                upperboundexceedperiod = [];
                if lowerbound < 0
                    lowerboundexceedperiod = [360 - lowerbound, 360];
                end
                if upperbound > 360
                    upperboundexceedperiod = [0, upperbound - 360];
                end
                includeperiods = [incircleperiod; lowerboundexceedperiod; upperboundexceedperiod];
                preferredangles = find(isExcluded(x_bins, includeperiods));
                preferredangles_spikecounts = binnedphase(preferredangles);
                baselineangles = find(~isExcluded(x_bins, includeperiods));
                baselineangles_spikecounts = binnedphase(baselineangles);
                [h, p, ci, stat] = ttest2(reshape(preferredangles_spikecounts, 1, []), reshape(baselineangles_spikecounts, 1, []), ...
                    "Alpha", 0.05, "Tail", "right", "Vartype", "unequal");
                % hs = [hs, h];
                pvals = [pvals, p];
                tstats = [tstats, stat.tstat];
            end
            % potentialsubtypes = subtypes(find(hs));
            [fdr, qval] = mafdr(pvals);
            hs = find(qval < qval_thr); % select subtypes passing qval threshold
            % select subtype(s) with largest t statistics
            % if largest 3 t stats differ less than 0.5, then assign 3
            % subtypes
            if apply_tstat_thr == 1
                if length(hs) > 1
                    [sortedtstats, sortedtstatsI] = sort(tstats, 'descend');
                    if sortedtstats(1) - sortedtstats(2) <= tstat_thr
                        if sortedtstats(1) - sortedtstats(3) <= tstat_thr
                            hs = intersect(hs, sortedtstatsI(1:3));
                        else
                            hs = intersect(hs, sortedtstatsI(1:2));
                        end
                    else
                        hs = sortedtstatsI(1);
                    end
                end
            end
            potentialsubtypes = subtypes(hs);
            potentialsubtypes_qval = qval(hs);
            potentialsubtypes_tstat = tstats(hs);
            
            bar(x_bins, binnedphase)
            ylim([0, max(binnedphase)*14/10])
            format short
            for ct = 1:length(hs)
                ctinfo = ['class=', char(potentialsubtypes(ct)), ', adj\_p=', num2str(potentialsubtypes_qval(ct)), ', tstat=' num2str(potentialsubtypes_tstat(ct))];
                text(0, max(binnedphase)*(14-ct)/10, ctinfo, "FontSize", 14);
            end
            tstat_dirname = ['qvalthr' num2str(qval_thr) 'tstatthr' num2str(tstat_thr)];
            if ~exist(fullfile(SavePath, 'phase2theta', 'tstat', tstat_dirname))
                mkdir(fullfile(SavePath, 'phase2theta', 'tstat', tstat_dirname))
            end
            saveas(gcf, fullfile(SavePath, 'phase2theta', 'tstat', tstat_dirname, ['NStstat_day' num2str(i) '_cell' num2str(iCell) '.png']))
            if ~isempty(potentialsubtypes)
                NeuronFile(i).ThetaBasedCellType{iCell} = potentialsubtypes;
                % disp(sum(NeuronFile(i).BinnedThetaPhase{iCell})) % total spike counts in theta period, to determine spike thr
            else 
                NeuronFile(i).ThetaBasedCellType{iCell} = ["ungrouped_ns_pv"];
            end
        end
    end
end

        
% stats
numtypes = [];
assignedtypes = [];
for i = 1:length(NeuronFile)
    for iCell = 1:length(NeuronFile(i).ThetaBasedCellType)
        if ~isempty(NeuronFile(i).ThetaBasedCellType{iCell})
            ctinfo = NeuronFile(i).ThetaBasedCellType{iCell};
            if length(ctinfo) == 1 & ctinfo == ["ungrouped_ns_pv"]
                numtypes = [numtypes, 0];
            else 
                numtypes = [numtypes, length(ctinfo)];
                assignedtypes_tmp = zeros(length(subtypes), 1);
                assignedtypes_tmp(ismember(subtypes, ctinfo)) = 1;
                assignedtypes = [assignedtypes, assignedtypes_tmp];
            end
        end
    end
end
numtypes_c = categorical(numtypes, [0 1 2], {'0' '1' '2'});
[n, c] = histcounts(numtypes_c);

pie(n)
title('Num of subtypes assigned')
lgd = legend(c);
saveas(gcf, fullfile(SavePath, 'phase2theta', 'tstat', 'subtype_num_of_assignment_pie.png'))

uniq_assignedtypes = unique(assignedtypes', 'rows');
type_categories = [];
type_counts = [];
for itype = 1:size(uniq_assignedtypes, 1)
    type_category = strjoin(subtypes(find(uniq_assignedtypes(itype, :))));
    type_categories = [type_categories, type_category];
    type_count = sum(ismember(assignedtypes', uniq_assignedtypes(itype, :), "rows"));
    type_counts = [type_counts, type_count];
    disp([type_category, num2str(type_count)])
end

b = bar(categorical(type_categories), type_counts);
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);
text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
saveas(gcf, fullfile(SavePath, 'phase2theta', 'tstat', 'subtype_counts.png'))
