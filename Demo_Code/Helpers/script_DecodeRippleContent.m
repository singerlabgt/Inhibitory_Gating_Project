%script_decode_ripple_content.m [EXAMPLE ONLY]
%
% ALP 12/09/21
% adapted from script_bayesiandecoding_rippples ALP 5/26/22
% NJ 08/19/22 adapted from Abby's code for novelty project - uses spike
% ephys indices for ripple timing
%
% This script uses a single time window per sharp-wave ripple event to
% decode the most likely information the SWR has about the distance to the
% closest reward zone. Outputs a giant 'RipData' output structure with all
% ripple events and the bins with the highest spatial probability.
% This script requires access to the Singer Lab ProcessedData folder on
% server for the extracted 'ripples' and 'bestRippleChan' matfiles. Also
% requires 'placecodingoutput' matfile from script_PlaceCodingProperties.m


clear; close all;

%% update your main directory to manuscript folder 
maindir = '\\ad.gatech.edu\bme\labs\singer\Nuri\Code\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%

%load default parameters
[dirs, params] = getDefaultParameters(maindir); 

%load cell type info
load(fullfile(dirs.data2load, 'cell_metrics.mat'));

%load NeuronFile data structure 
load(fullfile(dirs.data2load, getlatestfile_with_string(dirs.data2load, 'NeuronFile') ))

%get all session index
allindex = getSessionIndex(dirs.spreadsheet, params);
allindex = table2array(allindex);
allindex(allindex(:,6) > 3, :) = []; %remove control VR sessions
sessions = unique(allindex(:,[1:2 7]),'rows');

%set parameters 
params.nDeg = 100; %in deg
if params.nDeg > 300
    trialT = 'Laps';
else
    trialT = 'Distance2RZ';
end
params.nonThetaRipples = 0; %1 for nontheta only ripples, 0 for all ripples
params.pfBinsize = 5; %in deg
params.timeAroundMid = 0.125; %in s
params.decodingBin_s = params.timeAroundMid * 2; %use one timebin per ripple (i.e. do not break into multiple timebins)
params.decodingBin_samp = params.decodingBin_s .* params.samprate;


%% run decoding
for iSess = 1:size(sessions,1)
    tempindex = allindex(allindex(:,1) == sessions(iSess,1) & allindex(:,2) == sessions(iSess,2),:);
    tempData = []; 
    for ii = 1:size(tempindex,1)
        index = tempindex(ii,:); 
        switch index(6)
            case 1
                env = 'fam';
            case 2
                env = 'nov';
            case 3
                env = 'nov';
            otherwise
                disp('environment unidentified')
        end
        anprocesseddatadir = fullfile(dirs.processeddata,...
            [params.iden num2str(index(1)) '_' num2str(index(2))], params.brainReg{1});
        
        %%%%% ---------- load files -------------        
        %%% load spike indicies and ratemap
        spikeData = NeuronFile(iSess).SpikeIdx(ii).data;
        ratemapdir = fullfile('Y:\singer\Nuri\OutputStructs\placecoding', trialT);
        fname = getlatestfile_with_string(ratemapdir,...
            [params.iden num2str(index(1)) '_' num2str(index(2))]);
        load(fullfile(ratemapdir, fname))
        
        
        %%% load best ripple channel
        ripplechanfile = fullfile(anprocesseddatadir,...
            ['bestRippleChan' num2str(index(1)) '_' num2str(index(2)) '.mat']);
        if isfile(ripplechanfile)
            load(ripplechanfile)
        end
        rippleChan = bestRippleChan.channel;
        
        %%% load lfp
        anperiodsdatadir = fullfile(anprocesseddatadir, num2str(rippleChan));
        load(fullfile(anperiodsdatadir, ['ripples' num2str(index(3)) '.mat']))
        % give structures easy names for interfacing
        SWRData = ripples{index(1)}{index(2)}{index(3)};
        
        
        if ~isempty(SWRData.startind) && isfield(outmap.(['bin' num2str(params.pfBinsize)]).(env), 'ratemap_original')
            load(fullfile(anperiodsdatadir, ['nonthetas' num2str(index(3)) '.mat']))
            load(fullfile(anperiodsdatadir, ['eeg', num2str(index(3)) '.mat']))
            NTData = nonthetas{index(1)}{index(2)}{index(3)};
            LFPData = eeg{index(1)}{index(2)}{index(3)};
            downsamp = params.samprate/SWRData.samprate;
            
            %%% what ripples should we use?
            if params.nonThetaRipples == 1
                nonthetaper = [NTData.startind NTData.endind];
                isNTRip = isExcluded(SWRData.midind, nonthetaper);
                isNTRip = logical(isNTRip);
            else
                isNTRip = logical(ones(1,length(SWRData.midind)))';
            end
            rippleDurS = (SWRData.endind - SWRData.startind)./SWRData.samprate; %in s
            rippleInds = [SWRData.midind(isNTRip).*downsamp - params.timeAroundMid.*params.samprate,...
                SWRData.midind(isNTRip).*downsamp + params.timeAroundMid.*params.samprate]; %look @250ms around the middle of the ripple
            
            %%% get position for each ripple estimate position using bayes
            tempData = [tempData, getdecodedrippletrajectory_bayes_NJ(...
                spikeData, outmap, rippleInds, rippleDurS, params, index, cell_metrics)];
            
        end
    end
    RipData(iSess).(trialT) = tempData;        
end
save(fullfile(dirs.data2load,['RipData_' num2str(params.decodingBin_s *1000) 'ms_' datestr(now,'yymmdd') '.mat']),'RipData','params','sessions','-v7.3')

