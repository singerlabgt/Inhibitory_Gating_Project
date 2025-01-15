function t = PV_subgrouping_fig1G_heatmap(params, celltypes, celltypeNames, cellT, allsess_sessinfo, allsess_mean_fam_residual_hit, position_binEdges)
%% Fig. 1G
% f = figure('units','inch','position',[50 150 200,800]);
figure('Position', [50 150 250,800])
t = tiledlayout(5, 1,'TileSpacing','Compact','Padding','Compact');
% t = tiledlayout(5, 1, 'TileSpacing','compact','Units','inches','OuterPosition',[0 0 6.5 4]);
colormap(gray);
for ct = 1:length(celltypes)
    %min-max normalization per unit
    vMax = max(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
    vMin = min(allsess_mean_fam_residual_hit(cellT{ct},:),[],2);
    normMap_fam = (allsess_mean_fam_residual_hit(cellT{ct},:) - vMin) ./ (vMax - vMin);
    % nexttile([3 2]);
    nexttile(t)
    temp = normMap_fam(ismember(allsess_sessinfo(cellT{ct},1), params.WTmice),:); %select WT units only
    temp = temp(all(~isnan(temp),2),:);
    % populationAvg{ct} = temp(all(~isnan(temp),2),:);
    if strcmp(celltypes{ct}, 'Pyramidal Cell')
        [~,maxLoc] = max(temp,[],2,'omitnan'); [~,b2] = sort(maxLoc); %sort pyramidal cells by peak location
    else
        [~,minLoc] = min(temp,[],2,'omitnan'); [~,b2] = sort(minLoc); %sort interneurons by min location
        % disp(length(b2))
    end
    imagesc(temp(b2,:));box off;
    if ct == length(celltypes); cb = colorbar('Location','westoutside'); cb.Label.String = 'Normalized residual firing rate'; end
    if ct == 3; ylabel('Unit index'); end
    if ct == 5; xlabel('Distance to familiar RZ (deg)'); end
    % 
    tick0 = position_binEdges(1);
    tick_mid = position_binEdges(ceil(length(position_binEdges)/2));
    tick_end = -1 * position_binEdges(1);
    set(gca, 'XTick',[1,ceil(length(position_binEdges)/2),length(position_binEdges)-1],...
        'XTickLabel',[tick0, tick_mid, tick_end]);
    
    title(celltypeNames{ct});
end
