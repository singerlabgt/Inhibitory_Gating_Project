function swr_centered_rate = out2binnedSpikeCounts(is_shuffled, out, cell_metrics, allindex, dirs, params, identifier, ident_peranimal, animals, sessions, downsamplefactor, ephys_samprate, ...
    secBefore, secAfter, minRips, pre2start_binnum, start2mid_binnum, mid2end_binnum, end2post_binnum, iNarrow)
% input ripple out.mat and allindex
% output binned spiking rate centered at midind of ripple during the day.
% The time is normalized to m bins pre-start/end-post, n bins
% start-mid/mid-end.
% TODO: spikeThr
if is_shuffled
    repeat_time = length(out.sessindex) / length(out.rip_startind);
else
    repeat_time = 1;
end
for iSess = 1:height(allindex)
    currSess = allindex(iSess, :);
    animal = currSess.Animal; 
    recDay = currSess.Date;
    currFile = currSess.Recording;
    
    %% binned spikes centered at SWR midpoint
    %load cluster file
    animal_fname = [identifier{iSess}, num2str(animal) '_' num2str(recDay)];
    fdir = fullfile(dirs.processeddata_nj, animal_fname, params.brainReg{1});
    rawclusters = load(fullfile(fdir, dirs.clusfolder, ['rawclusters' num2str(currFile) '.mat']));
    rawclusters = rawclusters.rawclusters;
    ripindex = find(out.sessindex(:,1) == animal ...
        & out.sessindex(:,2) == recDay & out.sessindex(:,3) == currFile);
    if ~is_shuffled
        midripindex = out.rip_midind(ripindex) .* downsamplefactor;
        startripindex = out.rip_startind(ripindex) .* downsamplefactor;
        endripindex = out.rip_endind(ripindex) .* downsamplefactor;
    else
        midripindex = out.shuffled_midind(ripindex) .* downsamplefactor;
        startripindex = out.shuffled_startind(ripindex) .* downsamplefactor;
        endripindex = out.shuffled_endind(ripindex) .* downsamplefactor;
    end
    pre_rip = midripindex - secBefore *ephys_samprate; %X sec before rip midpoint
    post_rip = midripindex + secAfter *ephys_samprate;
    rip_start = startripindex - 0.035 *ephys_samprate; %35ms before rip start
    rip_end = endripindex - 0.035 *ephys_samprate; 
    % calculate 4 time invtervals
    pre2start_t = startripindex - pre_rip;
    start2mid_t = midripindex - startripindex;
    mid2end_t = endripindex - midripindex;
    end2post_t = post_rip - endripindex;
    ripstats{iSess'} = [pre2start_t(1:repeat_time:end), start2mid_t(1:repeat_time:end), mid2end_t(1:repeat_time:end), end2post_t(1:repeat_time:end)];

    %spike indices nUnits x nRipples (X sec before and after midpoint)
    clear temp_pre2start temp_start2mid temp_mid2end temp_end2post
    for iUnit = 1:length(rawclusters) % CAN SET THE ZERO POINT AS THE STARTGING POINT
        spikeInds = rawclusters(iUnit).spikeInds;
        temp_pre2start(iUnit,:) = arrayfun(@ (x) spikeInds(find(isExcluded(spikeInds,...
            [pre_rip(x), startripindex(x)]))) - midripindex(x), 1:length(pre_rip), 'UniformOutput', false);
        temp_start2mid(iUnit,:) = arrayfun(@ (x) spikeInds(find(isExcluded(spikeInds,...
            [startripindex(x), midripindex(x)]))) - midripindex(x), 1:length(pre_rip), 'UniformOutput', false);
        temp_mid2end(iUnit,:) = arrayfun(@ (x) spikeInds(find(isExcluded(spikeInds,...
            [midripindex(x), endripindex(x)]))) - midripindex(x), 1:length(pre_rip), 'UniformOutput', false);
        temp_end2post(iUnit,:) = arrayfun(@ (x) spikeInds(find(isExcluded(spikeInds,...
            [endripindex(x), post_rip(x)]))) - midripindex(x), 1:length(pre_rip), 'UniformOutput', false);
    end
    spikeseries_pre2start{iSess'} = temp_pre2start;
    spikeseries_start2mid{iSess'} = temp_start2mid;
    spikeseries_mid2end{iSess'} = temp_mid2end;
    spikeseries_end2post{iSess'} = temp_end2post;
%     spikesduringSWR{iSess'} = temp_spikecount;
    sessidx = find(startsWith(cell_metrics.sessionName, [ident_peranimal{animals == currSess.Animal} num2str(currSess.Animal) '_' num2str(currSess.Date)]));
    spikeAbsID{iSess'} = intersect(sessidx, find(ismember(cell_metrics.cluID, [rawclusters.ID])));% absolute ID according to cell_metrics
end

%include sessions with at least X ripple periods only
nRips = arrayfun( @(x) size(spikeseries_pre2start{x},2), 1:length(spikeseries_pre2start))' / repeat_time;
sess2incl = find(nRips >= minRips); % minRips=0, include all sessions


%% convert each cell from per-session to per-day
clear swr_centered_rate swr_centered_std spikeseries_pre2start_byday spikeseries_start2mid_byday spikeseries_mid2end_byday spikeseries_end2post_byday ripstats_byday spikeAbsID_byday
for iDay = 1:length(sessions)
    animal = sessions(iDay, 1);
    recDay = sessions(iDay, 2);
    all_sess_idx = intersect(sess2incl, find(allindex.Animal == animal & allindex.Date == recDay));
    spikeseries_pre2start_byday{iDay} = [spikeseries_pre2start{all_sess_idx}];
    spikeseries_start2mid_byday{iDay} = [spikeseries_start2mid{all_sess_idx}];
    spikeseries_mid2end_byday{iDay} = [spikeseries_mid2end{all_sess_idx}];
    spikeseries_end2post_byday{iDay} = [spikeseries_end2post{all_sess_idx}];
    ripstats_byday{iDay} = vertcat(ripstats{all_sess_idx});
    spikeAbsID_byday{iDay} = unique([spikeAbsID{all_sess_idx}]);
end
clear spikeseries_pre2start spikeseries_start2mid spikeseries_mid2end spikeseries_end2post spikesduringSWR spikeAbsID ripstats
%% for each ripple calculate the histgram for each unit
%defined edges for SWR-midpoint-centered firing rates
for iDay = 1:length(sessions)
    NSunits = intersect(spikeAbsID_byday{iDay}, iNarrow); % if only narrow PV units are used
    temp = cell(length(NSunits), 4);
    for iRip = 1:size(ripstats_byday{iDay}, 1)
        ripstat = ripstats_byday{iDay}(iRip, :);
        start2mid_binsize = ripstat(2) / start2mid_binnum;
        mid2end_binsize = ripstat(3) / mid2end_binnum;
        edges_pre2start = -1 * (pre2start_binnum + start2mid_binnum) * start2mid_binsize : start2mid_binsize : -1 * start2mid_binnum * start2mid_binsize; % 10 bins, same bin size as star2mid
        edges_start2mid = -1 * start2mid_binnum * start2mid_binsize : start2mid_binsize : 0;
        edges_mid2end = 0 : mid2end_binsize : mid2end_binnum * mid2end_binsize;
        edges_end2post = mid2end_binnum * mid2end_binsize : mid2end_binsize : (mid2end_binnum + end2post_binnum) * mid2end_binsize;

        if ~isnan(ripstat(1)) 
            start2mid_binsize = ripstat(2) / start2mid_binnum;
            mid2end_binsize = ripstat(3) / mid2end_binnum;
            edges_pre2start = -1 * (pre2start_binnum + start2mid_binnum) * start2mid_binsize : start2mid_binsize : -1 * start2mid_binnum * start2mid_binsize; % 10 bins, same bin size as star2mid
            edges_start2mid = -1 * start2mid_binnum * start2mid_binsize : start2mid_binsize : 0;
            edges_mid2end = 0 : mid2end_binsize : mid2end_binnum * mid2end_binsize;
            edges_end2post = mid2end_binnum * mid2end_binsize : mid2end_binsize : (mid2end_binnum + end2post_binnum) * mid2end_binsize;

            for iUnit = 1:length(NSunits)
                temp{iUnit, 1} = [temp{iUnit, 1};histcounts(vertcat(spikeseries_pre2start_byday{iDay}{iUnit,repeat_time*(iRip-1) + 1 : 1 : repeat_time*iRip}), edges_pre2start) ./ repeat_time];
                temp{iUnit, 2} = [temp{iUnit, 2}; histcounts(vertcat(spikeseries_start2mid_byday{iDay}{iUnit,repeat_time*(iRip-1) + 1 : 1 : repeat_time*iRip}), edges_start2mid) ./ repeat_time];
                temp{iUnit, 3} = [temp{iUnit, 3}; histcounts(vertcat(spikeseries_mid2end_byday{iDay}{iUnit,repeat_time*(iRip-1) + 1 : 1 : repeat_time*iRip}), edges_mid2end) ./ repeat_time];
                temp{iUnit, 4} = [temp{iUnit, 4}; histcounts(vertcat(spikeseries_end2post_byday{iDay}{iUnit,repeat_time*(iRip-1) + 1 : 1 : repeat_time*iRip}), edges_end2post) ./ repeat_time];
            end
        end
    end
    rate = cellfun(@sum, temp, 'UniformOutput', false);
    rate_normalized = cellfun(@mean, temp, 'UniformOutput', false); % for backup usage. not mentioned in the paper
    std = cellfun(@std, temp, 'UniformOutput', false);
    
    swr_centered_rate{iDay, 1} = temp;
    swr_centered_rate{iDay, 2} = NSunits;

end