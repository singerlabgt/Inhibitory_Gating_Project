clear; close all;

% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 
load(fullfile(dirs.data2load, 'cell_metrics.mat'));

%get all indices
allindex = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %exclude VR manipulation sessions
allindex = allindex(ismember(allindex(:,1), params.animals),:);

mice2incl = params.goalshamMice;
% only include nov sessions
goalstimsess = unique(allindex(ismember(allindex(:,1),mice2incl) & allindex(:,6) ~= 1 & allindex(:,9) == 1, [1:2,7]),'rows'); 
shamstimsess = unique(allindex(ismember(allindex(:,1),mice2incl) & allindex(:,6) ~= 1 & allindex(:,9) == 2, [1:2,7]),'rows');

conditions = {'goal', 'sham'};
% 3 x 2 plots: NS int. firing rates at last N stim trials vs. N no-stim
% trials. plot effects of all light intensities.
for iCond = 1%:2 % only look at goal-stim
    post_highstim_num = []; post_lowstim_num = []; post_nostim_num = [];
    cond = conditions{iCond};
    if iCond == 1
        sess2incl = goalstimsess;
    else
        sess2incl = shamstimsess;
    end
    highstim_ratemap = [];
    lowstim_ratemap = [];
    nostim_ratemap = [];
    poststim_ratemap = [];
    post_highstim_ratemap = [];
    post_nostim_ratemap = [];
    post_lowstim_ratemap = [];
    base_ratemap = [];
    for iSess = 1:size(sess2incl, 1)
        if sess2incl(iSess, 1) == 4
            iden = 'X';
        else
            iden = 'N';
        end
        sessname = [iden num2str(sess2incl(iSess, 1)) '_' num2str(sess2incl(iSess, 2)) '_CA3'];
        celltype = cell_metrics.putativeCellType(find(ismember(cell_metrics.sessionName, sessname)));
        notbadidx = ~ismember(find(ismember(cell_metrics.sessionName, sessname)), cell_metrics.tags.Bad);
        unit2incl = intersect(find(notbadidx), find(ismember(celltype, 'Narrow Interneuron')));
        if isempty(unit2incl)
            continue
        end
        
        fname = dir(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ", [iden num2str(sess2incl(iSess, 1)) '_' num2str(sess2incl(iSess, 2)) '*mat']));
        load(fullfile(fname.folder, fname.name));
        % residual fr calculation
        x1 = reshape([outmap.bin5.fam.smoothSpeed; outmap.bin5.nov.smoothSpeed],[],1); %predictor 1: speed
        x2 = reshape([outmap.bin5.fam.smoothLickrate; outmap.bin5.nov.smoothLickrate],[],1); %predictor 2: lick rate
        X = [ones(size(x1)) x1 x2 x1.*x2];
        ratemap = nan(size(outmap.bin5.nov.ratemap));
        %for each unit regardless of celltype
        for ii = 1:length(outmap.unitID)
            f_rate = outmap.bin5.fam.ratemap(:,:,ii);
            n_rate = outmap.bin5.nov.ratemap(:,:,ii);
            f_size = size(f_rate); n_size = size(n_rate);
            combined_size = size([f_rate; n_rate]);
            
            y = reshape([f_rate;n_rate], [], 1); %combine familiar and novel data
            
            [b,bint,r,rint,stats] = regress(y,X);
            tempRes = reshape(r, combined_size);
            % famRes = tempRes(1:f_size(1),:);
            novRes = tempRes(f_size(1)+1:end,:);
            ratemap(:,:,ii) = novRes(:,:);
        end
        % separate stim/poststim trials
        RZcenteredtrials = find(outmap.bin5.nov.labels(:, 2) == 1);  % RZ-centered. TODO: filter correct trials later
        hit = outmap.bin5.nov.labels(RZcenteredtrials, 3) == 1;
        if isempty(RZcenteredtrials)
            continue
        end
        ratemap = ratemap(RZcenteredtrials, :, :);
        rcfileidx = outmap.bin5.nov.labels(RZcenteredtrials, 1); % separate post-stim trials per file
        [uniqfiles, uniqstartidx] = unique(rcfileidx);
        v = cellfun(@(x) max(x), outmap.bin5.nov.stimVoltage);
        v = v(RZcenteredtrials);

        highstim_trials = find(v > 1); 
        lowstim_trials = find(v > 0 & v < 1);
        
        post_nostim_trials = []; post_lowstim_trials = []; post_highstim_trials = [];
        for iFile = 1:length(uniqfiles)
            v_iFile = v(rcfileidx == uniqfiles(iFile));
            poststim_trialidx_start = find(v_iFile>0, 1, 'last' ) + 1;
            if mod(poststim_trialidx_start, 3) == 0
                poststim_trialidx_start = poststim_trialidx_start + 1; 
            end
            % separate post-RZ trials            
            vpattern = unique(v_iFile);
            for iPattern = 1:length(vpattern)
                % RZidx_iFile_iPattern = RZidx_iFile(find(v_iFile == vpattern(iPattern), 1));
                poststim_iPattern_idx = setdiff(find(v_iFile == vpattern(iPattern), 1) : 3 : length(v_iFile), 1 : poststim_trialidx_start - 1 );% what should be the pattern after stim session?
                if vpattern(iPattern) > 1
                    post_highstim_trials = [post_highstim_trials, poststim_iPattern_idx +  uniqstartidx(iFile) - 1];
                elseif vpattern(iPattern) < 1 && vpattern(iPattern) > 0
                    post_lowstim_trials = [post_lowstim_trials, poststim_iPattern_idx +  uniqstartidx(iFile) - 1];
                elseif vpattern(iPattern) == 0
                    post_nostim_trials = [post_nostim_trials, poststim_iPattern_idx +  uniqstartidx(iFile) - 1];
                end
            end
            
        end
        % group all post-RZ trials
        poststim_trials = [post_highstim_trials, post_lowstim_trials, post_nostim_trials];
        nostim_trials = setdiff(find(v == 0), poststim_trials);
        % norm and scaling
        x_axis = outmap.bin5.binEdges(1:end-1);
        % RZtrialstart = find(x_axis == -40); RZtrialend = find(x_axis == 40);
        RZtrialstart = 1; RZtrialend = length(x_axis);
        x_axis_new = x_axis(RZtrialstart:RZtrialend);

        highstimfr = nanmean(ratemap(intersect(find(hit), highstim_trials), RZtrialstart:RZtrialend, unit2incl), 1);
        lowstimfr = nanmean(ratemap(intersect(find(hit), lowstim_trials), RZtrialstart:RZtrialend, unit2incl), 1);
        nostimfr = nanmean(ratemap(intersect(find(hit), nostim_trials), RZtrialstart:RZtrialend, unit2incl), 1);
        post_highstimfr = nanmean(ratemap(intersect(find(hit), post_highstim_trials), RZtrialstart:RZtrialend, unit2incl), 1);
        post_lowstimfr = nanmean(ratemap(intersect(find(hit), post_lowstim_trials), RZtrialstart:RZtrialend, unit2incl), 1);
        post_nostimfr = nanmean(ratemap(intersect(find(hit), post_nostim_trials), RZtrialstart:RZtrialend, unit2incl), 1);
        post_stimfr = nanmean(ratemap(intersect(find(hit), poststim_trials), RZtrialstart:RZtrialend, unit2incl), 1);
        basefr = nanmean(ratemap(find(hit), :, unit2incl), 1);

        post_highstim_num = [post_highstim_num, length(intersect(find(hit), highstim_trials))];
        post_lowstim_num = [post_lowstim_num, length(intersect(find(hit), lowstim_trials))];
        post_nostim_num = [post_nostim_num, length(intersect(find(hit), nostim_trials))];

        highstim_ratemap = cat(3, highstim_ratemap, (highstimfr ) );
        lowstim_ratemap = cat(3, lowstim_ratemap, (lowstimfr ) );
        nostim_ratemap = cat(3, nostim_ratemap, (nostimfr ) );
        poststim_ratemap = cat(3, poststim_ratemap, (post_stimfr ) ); % grouping all post-stim
        post_lowstim_ratemap = cat(3, post_lowstim_ratemap, (post_lowstimfr ) );
        post_highstim_ratemap = cat(3, post_highstim_ratemap, (post_highstimfr ) );
        post_nostim_ratemap = cat(3, post_nostim_ratemap, (post_nostimfr ) );
        base_ratemap = cat(3, base_ratemap, basefr);
        
    end
    frmat = {}; 
    frmat{1} = highstim_ratemap; frmat{2} = lowstim_ratemap; frmat{3} = nostim_ratemap; 
    frmat{4} = post_highstim_ratemap; frmat{5} = post_lowstim_ratemap; frmat{6} = post_nostim_ratemap; 
    frmat{7} = poststim_ratemap; frmat{8} = base_ratemap;
    frnames = {'highstim', 'lowstim', 'nostim', 'post-high', 'post-low', ...
        'post-no', 'post-all', 'basefr'};
    % noramlization and scaling
    
    for iIntensity = 1:8
        frmat{iIntensity} = (squeeze(frmat{iIntensity}))';
    end
    for iIntensity = 1:8 % norm
        % frmat{iIntensity} = (frmat{iIntensity} - min(frmat{iIntensity}, [], 2)) ./ (max(frmat{iIntensity}, [], 2) - min(frmat{iIntensity}, [], 2));
        frmat{iIntensity} = (frmat{iIntensity} - min(frmat{end}, [], 2)) ./ (max(frmat{end}, [], 2) - min(frmat{end}, [], 2));
    end
    for iIntensity = 1:8 % scale
        frmat{iIntensity} = (frmat{iIntensity} - nanmean(frmat{end}(:,1:2), 2)) .* 100;
    end
    % plot
    figure; t = tiledlayout(3,3); 
    colors = params.colors_narrowInt;
    for iIntensity = 1:7
        ax = nexttile;
        ops.ax     = ax;
        
        ops.x_axis = x_axis_new;
        ops.color_area = colors(end,:);
        ops.color_line = colors(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 2;
        ops.error      = 'sem';
        temp = frmat{iIntensity};
        plot_areaerrorbar(temp, ops); hold on; box off;
        
        %t-test with Bonferroni correction for multiple comparisons
        h = nan(size(temp,2),1);
        p = nan(size(temp,2),1);
        stats = cell(size(temp,2),1);
        for iT = 1:size(temp,2)
            [h(iT),p(iT),~,stats{iT}] = ttest(temp(:,iT));
        end
        sig = find(p < 0.05 / numel(p));
        arrayfun( @(ii) scatter(ax, outmap.bin5.binEdges(sig(ii)), 1,...
            'Marker','|','MarkerEdgeColor',colors(end,:)),1:length(sig));
        title([frnames{iIntensity} cond])
        ylabel('resid FR (%)'); xlabel('Dist2RZ')
        xline(0, '--'); yline(0, '--');% ylim([-20 40]) %xlim([-40 40])
    end
    %%% save fig
    makefigurepretty(gcf)
    figname = 'SuppFigure06';
    savefigALP([figdir '/'], figname, 'filetype', 'pdf')

    figure; t = tiledlayout(3,3); 
    colors = params.colors_narrowInt;
    for iIntensity = 1:7
        ax = nexttile;
        temp = frmat{iIntensity};
        plot(ax, x_axis_new, temp);
        ylabel('resid FR (%)'); xlabel('Dist2RZ')
        xline(0, '--'); yline(0, '--')
        title(frnames{iIntensity})
    end
end


