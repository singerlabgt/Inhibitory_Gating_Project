%% 08062024: plot raw behavioral data for 2 animal examples
%% 08212024: xz added lap by lap data. show this in MS supplement
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

%% plot 2 FAM examples, velocity

% RZ-centered avg data
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
        for iDay = 1%:length(date_iNovday)
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", allfiles{fileidx}))
            trialidx = outmap.bin5.fam.labels(:, 2) == 1; % not necessarily got a reward. just make sure it's NRZ-centered
            acrossdaymat = [acrossdaymat; outmap.bin5.fam.smoothSpeed(trialidx, :)];
        end
        ops.ax     = ax;
        ops.x_axis = mean(getBinEdges(outmap.bin5.binEdges),2);
        ops.color_area = params.colors_fam(end,:);
        ops.color_line = params.colors_fam(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 1;
        ops.error      = 'sem';
        plot_areaerrorbar(acrossdaymat, ops); hold on; box off; 
        xline(ax, 0, 'k:');
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Velocity (deg/s)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to familiar RZ (deg)');
        end
        ylim(ax, [5 25])
        xticks(ax, [-50 0 50])
        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end

    disp(['Familiar Velocity, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure03_WT_Fam_Velocity';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')


% 08212024 by xz: lap by lap trial data
disp('WT FAM')
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
            acrossdaymat = [acrossdaymat; outmap.bin5.fam.smoothSpeed];
        end
        % adapted from plotPerformance_xz.m
        velocCountsSmooth = acrossdaymat;
        vMax = nanmax(velocCountsSmooth(:));
        y = [0 0 30 30];
        x1 = [outmap.bin5.fam.thetaReward(1), outmap.bin5.fam.thetaReward(1)+10, outmap.bin5.fam.thetaReward(1)+10, outmap.bin5.fam.thetaReward(1)];
        x2 = [outmap.bin5.fam.thetaReward(2), outmap.bin5.fam.thetaReward(2)+10, outmap.bin5.fam.thetaReward(2)+10, outmap.bin5.fam.thetaReward(2)];
        x3 = [outmap.bin5.fam.thetaReward(3), outmap.bin5.fam.thetaReward(3)+10, outmap.bin5.fam.thetaReward(3)+10, outmap.bin5.fam.thetaReward(3)];
        hold on
        binsize=1; binEdges = 1:5:360;
        fill(ax, x1/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x2/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x3/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        plot(ax, binEdges, velocCountsSmooth', 'color',[.5 .5 .5]);
        disp(['animal=' iden num2str(anim) ', day=' num2str(iNovday) ', trialnum=' num2str(size(velocCountsSmooth, 1))])
        plot(ax, binEdges, nanmean(velocCountsSmooth), 'color', params.colors_fam(end,:),'linewidth',1);
        xlim(ax, [0 360])
        hold off 
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Velocity (deg/s)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to familiar RZ (deg)');
        end
        ylim(ax, [0 30])
        xticks(ax, [0 180 360])
        if iAnim == 1
            title(ax, ['Day ' num2str(iNovday)])
        end
    end

    disp(['Familiar Velocity, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure03_WT_Fam_PerTrial_Velocity';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% plot 2 NOV examples, velocity

% RZ-centered data
fig = figure('units','inch','position',[0 0 3 2]);
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", "*.mat")).name};

for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2);
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        hold on
        acrossdaymat = [];
        for iDay = 1%:length(date_iNovday) % ONLY PICK 1 NOV TRACK
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", allfiles{fileidx}))
            trialidx = outmap.bin5.nov.labels(:, 2) == 1; % not necessarily got a reward. just make sure it's RZ-centered
            acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothSpeed(trialidx, :)];
        end
        ops.ax     = ax;
        ops.x_axis = mean(getBinEdges(outmap.bin5.binEdges),2);
        ops.color_area = params.colors_nov(end,:);
        ops.color_line = params.colors_nov(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 1;
        ops.error      = 'sem';
        plot_areaerrorbar(acrossdaymat, ops); hold on; box off; 
        xline(ax, 0, 'k:');
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Velocity (deg/s)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to novel RZ (deg)');
        end
        ylim(ax, [10 20])
        xticks(ax, [-50 0 50])
        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end

    disp(['Novel Velocity, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure03_WT_Nov_Velocity';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

% 08212024 by xz: lap by lap trial data
disp('WT NOV')
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
            acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothSpeed];
        end
        velocCountsSmooth = acrossdaymat;
        vMax = nanmax(velocCountsSmooth(:));
        y = [0 0 30 30];
        x1 = [outmap.bin5.nov.thetaReward(1), outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)];
        x2 = [outmap.bin5.nov.thetaReward(2), outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)];
        x3 = [outmap.bin5.nov.thetaReward(3), outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)];
        hold on
        binsize=1; binEdges = 1:5:360;
        fill(ax, x1/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x2/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x3/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        plot(ax, binEdges, velocCountsSmooth', 'color',[.5 .5 .5]);
        disp(['animal=' iden num2str(anim) ', day=' num2str(iNovday) ', trialnum=' num2str(size(velocCountsSmooth, 1))])
        plot(ax, binEdges, nanmean(velocCountsSmooth), 'color', params.colors_nov(end,:),'linewidth',1);
        xlim(ax, [0 360])
        hold off 
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Velocity (deg/s)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to novel RZ (deg)');
        end
        xticks(ax, [0 180 360])
        if iAnim == 1
            title(ax, ['Day ' num2str(iNovday)])
        end
    end

    disp(['Novel Velocity, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure03_WT_Nov_PerTrial_Velocity';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')


%% plot 2 NOV goal-stim examples, velocity

% RZ-centered
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
animals = [57 65];
identifier = {'N', 'N'};
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", "*.mat")).name};
for iAnim = 1:length(animals)
    anim = animals(iAnim);
    iden = identifier{iAnim};
    date = sessions(sessions(:, 1) == anim, 2); disp(num2str(date))
    novelexposureday = sessions(sessions(:, 1) == anim, 3);
    for iNovday = 1:3
        date_iNovday = date(novelexposureday == iNovday);
        ax = nexttile;
        acrossdaymat = [];
        for iDay = 1%:length(date_iNovday) % ONLY PICK 1 NOV TRACK
            fileidx = find(cellfun(@(f) startsWith(f, [iden num2str(anim) '_' num2str(date_iNovday(iDay))]), allfiles));
            load(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", allfiles{fileidx}))
            trialidx = outmap.bin5.nov.labels(:, 2) == 1; % not necessarily got a reward. just make sure it's RZ-centered
            acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothSpeed(trialidx, :)];
        end
        % keep raw data. Don't do normalizatoin or scaling. 
        ops.ax     = ax;
        ops.x_axis = mean(getBinEdges(outmap.bin5.binEdges),2);
        ops.color_area = params.colors_goalstim(end,:);
        ops.color_line = params.colors_goalstim(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 1;
        ops.error      = 'sem';
        plot_areaerrorbar(acrossdaymat, ops); hold on; box off; 
        xline(ax, 0, 'k:');
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Velocity (deg/s)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to novel RZ (deg)');
        end
        ylim(ax, [2 8])
        xticks(ax, [-50 0 50])
        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end

    disp(['PV Goal Velocity, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure03_PV_Goal_Velocity';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

% 08212024 by xz: lap by lap trial data
disp('NOV goal-stim')
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
            acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothSpeed];
        end
        velocCountsSmooth = acrossdaymat;
        vMax = nanmax(velocCountsSmooth(:));
        y = [0 0 15 15];
        x1 = [outmap.bin5.nov.thetaReward(1), outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)];
        x2 = [outmap.bin5.nov.thetaReward(2), outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)];
        x3 = [outmap.bin5.nov.thetaReward(3), outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)];
        hold on
        binsize=1; binEdges = 1:5:360;
        fill(ax, x1/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x2/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x3/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        plot(ax, binEdges, velocCountsSmooth', 'color',[.5 .5 .5]);
        disp(['animal=' iden num2str(anim) ', day=' num2str(iNovday) ', trialnum=' num2str(size(velocCountsSmooth, 1))])
        plot(ax, binEdges, nanmean(velocCountsSmooth), 'color', params.colors_goalstim(end,:),'linewidth',1);
        xlim(ax, [0 360])
        hold off 
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Velocity (deg/s)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to RZ (deg)');
        end
        ylim(ax, [0 15])
        xticks(ax, [0 180 360])
        if iAnim == 1
            title(ax, ['Day ' num2str(iNovday)])
        end
    end

    disp(['PV Goal Velocity, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure03_PV_Goal_PerTrial_Velocity';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% plot 2 NOV sham-stim examples, velocity
% RZ-centered
fig = figure('units','inch','position',[0 0 3 2]);
t = tiledlayout(2, 3,'TileSpacing','compact');
animals = [57 65];
identifier = {'N', 'N'};
allfiles = {dir(fullfile(dirs.data2load, "singlesess_ratemap_distance2RZ/", "*.mat")).name};
iTile = 1;
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
        acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothSpeed(trialidx, :)];
  
        % keep raw data. Don't do normalizatoin or scaling. 
        ops.ax     = ax;
        ops.x_axis = mean(getBinEdges(outmap.bin5.binEdges),2);
        ops.color_area = params.colors_shamstim(end,:);
        ops.color_line = params.colors_shamstim(end,:);
        ops.alpha      = 0.2;
        ops.line_width = 1;
        ops.error      = 'sem';
        plot_areaerrorbar(acrossdaymat, ops); hold on; box off; 
        xline(ax, 0, 'k:');
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Velocity (deg/s)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to novel RZ (deg)');
        end
        ylim(ax, [2 8])
        xticks(ax, [-50 0 50])
        if iAnim == 1
        title(ax, ['Day ' num2str(iNovday)])
        end
    end
    disp(['PV Sham Velocity, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure03_PV_sham_Velocity';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

% 08212024 by xz: lap by lap trial data
disp('NOV sham-stim')
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
        acrossdaymat = [acrossdaymat; outmap.bin5.nov.smoothSpeed];
        velocCountsSmooth = acrossdaymat;
        vMax = nanmax(velocCountsSmooth(:));
        y = [0 0 15 15];
        x1 = [outmap.bin5.nov.thetaReward(1), outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)+10, outmap.bin5.nov.thetaReward(1)];
        x2 = [outmap.bin5.nov.thetaReward(2), outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)+10, outmap.bin5.nov.thetaReward(2)];
        x3 = [outmap.bin5.nov.thetaReward(3), outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)+10, outmap.bin5.nov.thetaReward(3)];
        hold on
        binsize=1; binEdges = 1:5:360;
        fill(ax, x1/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x2/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        fill(ax, x3/binsize, y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        plot(ax, binEdges, velocCountsSmooth', 'color',[.5 .5 .5]);
        disp(['animal=' iden num2str(anim) ', day=' num2str(iNovday) ', trialnum=' num2str(size(velocCountsSmooth, 1))])
        plot(ax, binEdges, nanmean(velocCountsSmooth), 'color', params.colors_shamstim(end,:),'linewidth',1);
        xlim(ax, [0 360])
        hold off 
        if iNovday == 1 && iAnim == 1
            ylabel(ax, {'Velocity (deg/s)'});
        elseif iNovday == 2 && iAnim == 2
            xlabel(ax, 'Distance to RZ (deg)');
        end
        yticks(ax, [0 5 10 15])
        xticks(ax, [0 180 360])
        if iAnim == 1
            title(ax, ['Day ' num2str(iNovday)])
        end
    end
    disp(['PV Sham Velocity, Mouse ' iden num2str(anim)])
end

%%% save figure
makefigurepretty(gcf,1)
figname = 'SuppFigure04_PV_Sham_PerTrial_Velocity';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

