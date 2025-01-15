%% getSessionIndex
% This function is adapted from getallindex_SpikeGadgetsNJ to simplify and
% implement array2table function
% NJ created 05.13.2020
% NJ updated 05.29.20 to make baseline files also accessible
%
%   INPUTS
%       spreadsheetfile - directory of spreadsheet containing session info
%       params - parameters containing animal info; output of getDirectoriesAndParams.m
%       OPTIONAL 3RD INPUT: 'baselineonly' - to output only baseline files in
%       between VR behavioral sessions
%
%   OUTPUT
%      allindex - table of all relevant animals and sessions; columns
%      associated with variable names as indicated in table

function [allindex, alliden] = getSessionIndex(spreadsheetfile, params, varargin)


sheetnumber = params.brainReg{1};   %eg. 'CA3'
animals2incl = params.animals;      %animal IDs to include
[data, text, ~] = xlsread(spreadsheetfile,sheetnumber);
allindex = array2table(data(:,1:9), 'VariableNames', text(1,2:10));


opts = detectImportOptions(spreadsheetfile);
opts.Sheet = sheetnumber;
ephysT = readtable(spreadsheetfile, opts);
alliden = ephysT.ID(2:end);


if nargin < 3
    sessionType = 'behavioronly'; %default, need no other specification
else
    sessionType = varargin{1}; %if want baseline files, 3rd input should be 'baselineonly'
end

switch sessionType
    case 'behavioronly'
        
        %only include behavioral files (by default) - good rec quality, no baseline behavior
        %files, and the animals of interest
        filter_idx = allindex.Include & allindex.VR ~= 10 & ...
            allindex.VR ~= 0 & ismember(allindex.Animal, animals2incl);
        allindex = allindex( filter_idx, :);
        alliden = alliden( filter_idx, :);
    case 'baselineonly'
        
        %only include baseline files with no visual cues (VR on but
        %projector covered with felt)
        filter_idx = allindex.Include & allindex.VR ~= 10 & ...
            allindex.VR == 0 & ismember(allindex.Animal, animals2incl);
        allindex = allindex(filter_idx , :);
        alliden = alliden( filter_idx, :);
    case 'all'
        
        %include both baseline and VR files 
        filter_idx = allindex.Include & allindex.VR ~= 10 & ...
            ismember(allindex.Animal, animals2incl);
        allindex = allindex( filter_idx, :);
        alliden = alliden( filter_idx, :);
    otherwise
        disp('Type of session not specified. Indicate baselineonly or behavioronly')
end
end

