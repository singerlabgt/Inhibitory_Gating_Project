function [dirs, params] = getDefaultParameters(userPath)
%This function spits out default parameters set for all analyses in the
%paper. Input requires user to specify the path to "Code" and "Data"
%folders. Outputs useful information such as sampling rate, binsize, color
%schemes for figure making. 

%% default directories 
dirs.data2load = fullfile(userPath, 'Demo_Data');
dirs.spreadsheet = fullfile(dirs.data2load, 'VR_NoveltySpreadsheet.xlsx'); 
dirs.saveoutputstruct = fullfile(userPath, 'Demo_Data' ,'IntermedOutput');
if ~isfolder(dirs.saveoutputstruct); mkdir(dirs.saveoutputstruct); end

%% general
params.iden = 'N';                                                          %typical animal ID starts with the letter N followed by a number
params.animals = [4, 11,18,21,24,45,46,47,48,50,52,53,54,57,61,63,62,65];      %animal IDs (numbers only)
params.WTmice = [11,18,21,24,45,46,47];                                     %animal IDs (numbers only) for all wild-type animals
params.goalshamMice = [4, 57,61,63,62,65];                                       %animal IDs (numbers only) for PVxAi32 animals exposed to both goal and sham stim
params.brainReg = {'CA3'};                                                  %brain region(s)
params.probeChannels = {1:64};                                              %64-channel NeuroNexus probe
params.binsize_ms = 1;                                                      %in ms for monoconnex calculations
params.binsize_deg = 2;                                                     %in degrees for behavioral analyses
params.binsize_s = 60;                                                      %in sec for firing rate stability across session
params.samprate = 30000;                                                    %in Hz for ephys acquisition system

%% color specification for grouped figures 
n_colors = cbrewer('seq','Greens',12);
params.colors_nov = n_colors(6:2:end, :);                                   %shades of green for novel day 1-3
f_colors = cbrewer('seq','Greys',12);
params.colors_fam = f_colors(7:2:end, :);                                   %shades of grey for familiar day 1-3
params.colors_shamstim = hex2rgb({'#ffba61','#ffa42e','#fa8d00'});          %shades of orange for sham stim 
params.colors_goalstim = hex2rgb({'#61a5ff','#2e88ff','#006cfa'});          %shades of blue for goal stim
params.colors_narrowInt = hex2rgb({'#d2dbec','#7894c5','#36489b'});         %shades of dark blue for NS interneuron
params.colors_wideInt = hex2rgb({'#e3f4f7','#aadde8','#72c7d9'});           %shades of lighter blue for WS interneuron
params.colors_pyr = hex2rgb({'#f5d6d6','#e08585','#cb3433'});               %shades of red for Pyramidal cell

%% place coding
params.environments = {'fam','nov'};
params.speedTh = 2;                                                         %in deg/s for occupancy normalized rate maps
params.spatialinfoTh = 95;                                                  %in n'th percentile for sptial info
params.cueSize = 10;                                                        %in degrees, size of wall cue (zone)

%% set whether to overwrite existing files or not
params.rewrite.ratemaps_lap = 0; 
params.rewrite.ratemaps_dist2rz = 0; 
params.rewrite.ratemaps_time2rz = 0; 

