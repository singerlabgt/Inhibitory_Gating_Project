%% 08062024: plot raw behavioral data for 2 animal examples
close all; clear
%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end

[dirs, params] = getDefaultParameters(maindir); 

[allindex, alliden] = getSessionIndex(dirs.spreadsheet, params, 'behavioronly');
allindex = table2array(allindex);
vronly_idx = allindex(:,6) > 3;
allindex(vronly_idx, :) = []; %exclude VR manipulation sessions
alliden(vronly_idx, :) = []; %exclude VR manipulation sessions
animal_idx = ismember(allindex(:,1), params.animals);
allindex = allindex(animal_idx,:); %filter based on animals to include 
[sessions, sessID] = unique(allindex(:,[1:2, 7]),'rows'); %session = [animalID, recording date, novelty day]

%% plot 2 FAM examples, lickrate
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
animals = [11 21];
identifier = {'N', 'N'};
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", "*.mat")).name};
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        acrossdaymat = [];
        for iDay = 1:length(date_iNovday)
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", allfiles{fileidx}))
            trialidx = outmap.bin5.fam.labels(:, 2) == 1; % not necessarily got a reward. just make sure it's NRZ-centered
            acrossdaymat = [acrossdaymat; outmap.bin5.fam.smoothLickrate(trialidx, :)];
        end
        ops.ax     = ax;
        ops.x_axis = mean(getBinEdges(outmap.bin5.binEdges),2);
        ops.color_area = params.colors_fam(end,:);
        ops.color_line = params.colors_fam(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 2;
        ops.error      = 'sem';
        plot_areaerrorbar(acrossdaymat, ops); hold on; box off;
        xline(ax, 0, 'k:');
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {['Lickrate (Hz)']});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to familiar RZ (deg)');
        end
        xticks(ax, [-50 0 50])
        ylim([0 4])
        yticks([0 2 4])
        if iAnim == 1
            title(ax, ['Day ' num2str(iNovday)])
        end
    end
end

makefigurepretty(gcf,1)
figname = 'SuppFigure04_WT_Fam_LickRate';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')


% 08212024 by xz: lap by lap trial data
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_laps/", "*.mat")).name};
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        acrossdaymat = [];
        for iDay = 1%:length(date_iNovday)
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_laps/", allfiles{fileidx}))
            acrossdaymat = [acrossdaymat; outmap.bin5.fam.smoothLickrate];
        end
        lickCountsSmooth = acrossdaymat;
        vMax = nanmax(lickCountsSmooth(:));
        y = [0 0 30 30];
        x1 = [outmap.bin5.fam.thetaReward(1), outmap.bin5.fam.thetaReward(1)+10, outmap.bin5.fam.thetaReward(1)+10, outmap.bin5.fam.thetaReward(1)];
        x2 = [outmap.bin5.fam.thetaReward(2), outmap.bin5.fam.thetaReward(2)+10, outmap.bin5.fam.thetaReward(2)+10, outmap.bin5.fam.thetaReward(2)];
        x3 = [outmap.bin5.fam.thetaReward(3), outmap.bin5.fam.thetaReward(3)+10, outmap.bin5.fam.thetaReward(3)+10, outmap.bin5.fam.thetaReward(3)];
        hold on
        binsize=1; binEdges = 1:5:360;
        fill(ax, x1/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x2/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x3/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        plot(ax, binEdges, lickCountsSmooth', 'color',[.5 .5 .5]);
        plot(ax, binEdges, nanmean(lickCountsSmooth), 'color', params.colors_fam(end,:),'linewidth',2);
        xlim(ax, [0 360])
        hold off 
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Lickrate (Hz)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Track position (bin)');
        end

        if iAnim == 1
        ylim(ax, [0 8])
        else
            ylim([0 6])
        end
                xticks([0 180 360])

        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end

        disp(['Familiar Lickrate, Mouse ' iden num2str(anim)])
    end


%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure04_WT_Fam_PerTrial_Lickrate';
savefigALP([figdir '/'], figname, 'filetype', 'pdf') 
%% plot 2 NOV examples, lickrate

% RZ-centered data
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", "*.mat")).name};
iTile = 1;
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile(iTile);
        acrossdaymat = [];
        for iDay = 1%:length(date_iNovday) % ONLY PICK 1 NOV TRACK
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", allfiles{fileidx}))
            trialidx = outmap.bin5.nov.labels(:, 2) == 1; % not necessarily got a reward. just make sure it's RZ-centered
            acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothLickrate(trialidx, :)];
        end
        ops.ax     = ax;
        ops.x_axis = mean(getBinEdges(outmap.bin5.binEdges),2);
        ops.color_area = params.colors_nov(end,:);
        ops.color_line = params.colors_nov(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 2;
        ops.error      = 'sem';
        plot_areaerrorbar(acrossdaymat, ops); hold on; box off; 
        xline(ax, 0, 'k:');
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {['Lickrate (Hz)']});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to novel RZ (deg)');
        end
        if iAnim == 1
        ylim(ax,[0 4])
        yticks([0 2 4])
        else
            ylim([0 2])
            yticks([0 1 2])
        end
                xticks(ax, [-50 0 50])
        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
        iTile = iTile + 1;
    end
end

makefigurepretty(gcf,1)
figname = 'SuppFigure04_WT_Nov_LickRate';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')


% 08212024 by xz: lap by lap trial data
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_laps/", "*.mat")).name};
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        acrossdaymat = [];
        for iDay = 1%:length(date_iNovday)
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_laps/", allfiles{fileidx}))
            acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothLickrate];
        end
        lickCountsSmooth = acrossdaymat;
        vMax = nanmax(lickCountsSmooth(:));
        y = [0 0 30 30];
        x1 = [outmap.bin5.nov.thetaReward(1), outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)];
        x2 = [outmap.bin5.nov.thetaReward(2), outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)];
        x3 = [outmap.bin5.nov.thetaReward(3), outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)];
        hold on
        binsize=1; binEdges = 1:5:360;
        fill(ax, x1/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x2/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x3/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        plot(ax, binEdges, lickCountsSmooth', 'color',[.5 .5 .5]);
        plot(ax, binEdges, nanmean(lickCountsSmooth), 'color', params.colors_nov(end,:),'linewidth',2);
        xlim(ax, [0 360])
        hold off 
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Lickrate (Hz)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Track position (bin)');
        end
        if iAnim == 1
            ylim(ax, [0 8])
        else
            ylim([0 6])
            yticks([0 2 4 6])
        end
                xticks([0 180 360])

        % xlim(ax, [x2(1) - 30, x2(2) + 30])
        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end

        disp(['WT Novel Lickrate, Mouse ' iden num2str(anim)])
    end


%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure04_WT_Novel_PerTrial_Lickrate';
savefigALP([figdir '/'], figname, 'filetype', 'pdf') 
%% plot 2 NOV goal-stim examples, lickrate

% RZ-centered data
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
animals = [57 65];
identifier = {'N', 'N'};
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", "*.mat")).name};
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        acrossdaymat = [];
        for iDay = 1%:length(date_iNovday) % ONLY PICK 1 NOV TRACK
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", allfiles{fileidx}))
            trialidx = outmap.bin5.nov.labels(:, 2) == 1; % not necessarily got a reward. just make sure it's RZ-centered
            acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothLickrate(trialidx, :)];
        end
        ops.ax     = ax;
        ops.x_axis = mean(getBinEdges(outmap.bin5.binEdges),2);
        ops.color_area = params.colors_goalstim(end,:);
        ops.color_line = params.colors_goalstim(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 2;
        ops.error      = 'sem';
        plot_areaerrorbar(acrossdaymat, ops); hold on; box off; 
        xline(ax, 0, 'k:');
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {['Lickrate (Hz)']});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to novel RZ (deg)');
        end
        ylim(ax,[0 5])
        yticks([0 2 4])
                xticks(ax, [-50 0 50])
        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end
end

makefigurepretty(gcf,1)
figname = 'SuppFigure04_PV_Goal_Lickrate';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

% 08212024 by xz: lap by lap trial data
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_laps/", "*.mat")).name};
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        acrossdaymat = [];
        for iDay = 1%:length(date_iNovday)
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_laps/", allfiles{fileidx}))
            acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothLickrate];
        end
        lickCountsSmooth = acrossdaymat;
        vMax = nanmax(lickCountsSmooth(:));
        y = [0 0 30 30];
        x1 = [outmap.bin5.nov.thetaReward(1), outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)];
        x2 = [outmap.bin5.nov.thetaReward(2), outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)];
        x3 = [outmap.bin5.nov.thetaReward(3), outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)];
        hold on
        binsize=1; binEdges = 1:5:360;
        fill(ax, x1/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x2/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x3/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        plot(ax, binEdges, lickCountsSmooth', 'color',[.5 .5 .5]);
        plot(ax, binEdges, nanmean(lickCountsSmooth), 'color', params.colors_goalstim(end,:),'linewidth',2);
        xlim(ax, [0 360])
        hold off 
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Lickrate (Hz)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Track position (bin)');
        end
        ylim(ax, [0 6])
        xticks([0 180 360])

        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end

        disp(['PV Goal Lickrate, Mouse ' iden num2str(anim)])
    end


%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure04_PV_Goal_PerTrial_Lickrate';
savefigALP([figdir '/'], figname, 'filetype', 'pdf') 

%% sham-stim examples, lickrate

% RZ-centered data
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
animals = [57 65];
identifier = {'N', 'N'};
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", "*.mat")).name};
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        acrossdaymat = [];
        iDay = 2;
        fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
        load(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", allfiles{fileidx}))
        trialidx = outmap.bin5.nov.labels(:, 2) == 1; % not necessarily got a reward. just make sure it's RZ-centered
        acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothLickrate(trialidx, :)];
        
        % keep raw data. Don't do normalizatoin or scaling. 
        ops.ax     = ax;
        ops.x_axis = mean(getBinEdges(outmap.bin5.binEdges),2);
        ops.color_area = params.colors_shamstim(end,:);
        ops.color_line = params.colors_shamstim(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 2;
        ops.error      = 'sem';
        plot_areaerrorbar(acrossdaymat, ops); hold on; box off; 
        xline(ax, 0, 'k:');
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {['Lickrate (Hz)']});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to novel RZ (deg)');
        end
        xticks(ax, [-50 0 50])
        ylim(ax,[0 5])
        yticks([0 2 4])

        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end
end

makefigurepretty(gcf,1)
figname = 'SuppFigure04_PV_Sham_Lickrate';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

% 08212024 by xz: lap by lap trial data
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_laps/", "*.mat")).name};
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        acrossdaymat = [];
        iDay = 2;
        fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
        load(fullfile(dirs.data2load, "singlesess_ratemap_laps/", allfiles{fileidx}))
        acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothLickrate];

        lickCountsSmooth = acrossdaymat;
        vMax = nanmax(lickCountsSmooth(:));
        y = [0 0 30 30];
        x1 = [outmap.bin5.nov.thetaReward(1), outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)];
        x2 = [outmap.bin5.nov.thetaReward(2), outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)];
        x3 = [outmap.bin5.nov.thetaReward(3), outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)];
        hold on
        binsize=1; binEdges = 1:5:360;
        fill(ax, x1/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x2/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x3/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        plot(ax, binEdges, lickCountsSmooth', 'color',[.5 .5 .5]);
        plot(ax, binEdges, nanmean(lickCountsSmooth), 'color', params.colors_shamstim(end,:),'linewidth',2);
        xlim(ax, [0 360])
        hold off 
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Lickrate (Hz)'});
        elseif iNovday == 2 && iAnim ==2
            xlabel(ax, 'Track position (bin)');
        end
        ylim(ax, [0 6])
        xticks([0 180 360])
        if iAnim == 1
            title(ax, ['Day ' num2str(iNovday)])
        end
    end
    disp(['PV Sham Lickrate, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure04_PV_Sham_PerTrial_Lickrate';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')