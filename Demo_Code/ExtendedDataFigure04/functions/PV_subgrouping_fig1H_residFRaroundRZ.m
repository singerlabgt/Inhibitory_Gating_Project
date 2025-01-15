function t = PV_subgrouping_fig1H_residFRaroundRZ(mode, params, celltypes, celltypeNames, populationAvg, baseBins, position_binEdges)
%% Fig. 1H: residual firing curves around RZ
% figure('units','inch','position', [50, 150, 200, 800])
figure('Position', [50 150 210,800])
t = tiledlayout(5, 1,'TileSpacing','Compact','Padding','Compact');
for ct = 1:length(celltypes)
    if strcmp(celltypes{ct}, 'aac')
        colors = params.colors_aac;
    elseif strcmp(celltypes{ct}, 'cck')
        colors = params.colors_cck;
    elseif strcmp(celltypes{ct}, 'bistratified')
        colors = params.colors_bistratified;
    elseif strcmp(celltypes{ct}, 'pvbc')
        colors = params.colors_pvbc;
    elseif strcmp(celltypes{ct}, 'ungrouped_ns_pv')
        colors = params.colors_ungrouped_ns_pv;
    end

    temp = populationAvg{ct};
    temp = (temp - nanmean(temp(:,baseBins),2)) .* 100;
    ax = nexttile;
    ops.ax     = ax;
    ops.x_axis = mean(getBinEdges(position_binEdges),2);
    ops.color_area = colors(end,:);
    ops.color_line = colors(end,:);
    ops.alpha      = 0.2;
    ops.line_width = 2;
    ops.error      = 'sem';
    plot_areaerrorbar(temp, ops); hold on; box off;

    %t-test with Bonferroni correction for multiple comparisons
    h = nan(size(temp,2),1);
    p = nan(size(temp,2),1);
    stats = cell(size(temp,2),1);
    for iT = 1:size(temp,2)
        [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT));
    end
    sig = find(p < 0.05 / numel(p));
    arrayfun( @(ii) scatter(ax, position_binEdges(sig(ii)), 10+ct*2,...
        'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
    % originally put out of the for loop -- modified by xiao 10/8/2023
    if mode == 'time'
        xlim([-6, 6])
        ylim([-30, 40])
        set(ax, 'XTick',[-6, 0, 6],...
            'XTickLabel',[-6, 0, 6]);
        xlabel(ax, 'Time to familiar RZ (sec)');
    elseif mode == 'dist'
        xlim(ax,[-40 40])
        ylim([-50, 50])
        set(ax, 'XTick',position_binEdges(find(ismember(position_binEdges,[-40,0,10,40]))),...
            'XTickLabel',{'-40','0','10','40'});
        xlabel(ax, 'Distance to familiar RZ (deg)');
    else
        disp('mode did not match dist or time.')
    end
    xline(ax, 0, 'k:'); xline(ax, 10, 'k:'); yline(ax, 0, 'k:');
    ylabel(ax, {'Change in residual firing','rate from baseline (%)'})
    
    title([celltypeNames{ct} ' n=' num2str(size(temp, 1))])
    %     %one sample/paired samples permutation t-test with correction for multiple comparisons
    %     p = arrayfun( @(x) mult_comp_perm_t1(temp(:,x),5000,0,0.05,0,0), 1:size(temp,2));
    %     sig = find(p < 0.05);
    %     arrayfun( @(x) scatter(ax, position_binEdges(sig(x)), 10+ct,'MarkerEdgeColor',colors(3,:)),1:length(sig))
end
