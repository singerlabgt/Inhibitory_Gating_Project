% The t-statistical approach:
% 1) For each potential subtypes (AAC, PVBC, CCK, bistratified), define the preferred phase as mentioned in the paper. 
% 2) Do one tail two-sample t-test, to compare the spike counts in the preferred/unpreferred phase. The length of each category corresponds to the preferred/unpreferred phase bin number (binwidth=18). Record p-values and t scores for each cell type t-test.
% 3) Do FDR estimation given the 4 cell types, and record adjusted p-values (q-values).
% 4) Thresholding: q-values<0.05. 
% If multiple cell types have q-values passing the threshold, then compare the first 2 cell types with largest t-scores:
% a. if t-score difference is larger than 0.5, then assign the cell type with the largest t-score to that neuron.
% b. if t-score difference <= 0.5, then assign both 2 cell types with the largest t-score to that neuron.

for i = 1:length(NeuronFile)
    PYR_notbad = intersect(find(ismember(NeuronFile(i).TSIDtotal, NeuronFile(i).TSIDnotbad)), find(NeuronFile(i).CellType == 1));
    for iCell = 1:length(NeuronFile(i).CellType)
        NeuronFile(i).BinnedThetaPhase_PYR{iCell} = [];
        NeuronFile(i).SpikeCountInTheta_PYR{iCell} = [];
        if ismember(iCell, PYR_notbad)
            data = NeuronFile(i).ThetaPhase{iCell};
            data = toDegrees('radians', data) + 180;
            
            h = histogram(data, 'binWidth', 1, 'BinLimits', [0 360]); 
            binnedphase = h.BinCounts;
            NeuronFile(i).BinSize1ThetaPhase_PYR{iCell} = binnedphase;

            h = histogram(data, 'binWidth', binwidth, 'BinLimits', [0 360]); 
            binnedphase = h.BinCounts;
            NeuronFile(i).BinnedThetaPhase_PYR{iCell} = binnedphase;
            NeuronFile(i).SpikeCountInTheta_PYR{iCell} = sum(binnedphase);
            x_bins = h.BinEdges(1:end-1)' + h.BinWidth/2;
        end
    end
end
%%
t0 = tiledlayout(round(length(NeuronFile)^(1/2)), ceil(length(NeuronFile)^(1/2)));
for i = 1:length(NeuronFile)
    theta_phases = [];
    for iCell = 1:length(NeuronFile(i).BinnedThetaPhase_PYR)
        if NeuronFile(i).SpikeCountInTheta_PYR{iCell} > spike_thr
            theta_phases = [theta_phases; (NeuronFile(i).BinnedThetaPhase_PYR{iCell}) ./ max(NeuronFile(i).BinnedThetaPhase_PYR{iCell})];
        end
    end
    if ~isempty(theta_phases)
        ax = nexttile;
        bar(x_bins, theta_phases)
        title(['day' num2str(i)])
    end
end
%%
close all
t1 = tiledlayout(round(length(NeuronFile)^(1/2)), ceil(length(NeuronFile)^(1/2)));
colors = params.colors_pvbc;
for i = 1:length(NeuronFile)
    theta_phases = [];
    for iCell = 1:length(NeuronFile(i).BinnedThetaPhase_PYR)
        if NeuronFile(i).SpikeCountInTheta_PYR{iCell} > spike_thr
            theta_phases = [theta_phases; (NeuronFile(i).BinnedThetaPhase_PYR{iCell}) ./ max(NeuronFile(i).BinnedThetaPhase_PYR{iCell})];
        end
    end
    if ~isempty(theta_phases)
        ax = nexttile;
        ops.ax     = ax;
        ops.x_axis = x_bins;
        ops.color_area = colors(end,:);
        ops.color_line = colors(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 2;
        ops.error      = 'sem';
        plot_areaerrorbar(theta_phases, ops); hold on; box off;
    end
end
%%
close all

t2 = tiledlayout(round(length(NeuronFile)^(1/2)), ceil(length(NeuronFile)^(1/2)));
for i = 1:length(NeuronFile)
    theta_phases = [];
    for iCell = 1:length(NeuronFile(i).BinnedThetaPhase_PYR)
        if NeuronFile(i).SpikeCountInTheta_PYR{iCell} > spike_thr
            theta_phases = [theta_phases; (NeuronFile(i).BinnedThetaPhase_PYR{iCell}) ./ max(NeuronFile(i).BinnedThetaPhase_PYR{iCell})];
        end
    end
    theta_phases = nanmean(theta_phases, 1);
    if ~isempty(theta_phases)
        ax = nexttile;
        bar(x_bins, theta_phases)
        title(['day' num2str(i)])
    end
end
%%
close all