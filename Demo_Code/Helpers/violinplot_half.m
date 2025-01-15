function [violins,plotinfo] = violinplot_half(data, cats, varargin)
%Violinplots plots violin plots of some data and categories
%   VIOLINPLOT(DATA) plots a violin of a double vector DATA
%
%   VIOLINPLOT(DATAMATRIX) plots violins for each column in
%   DATAMATRIX.
%
%   VIOLINPLOT(TABLE), VIOLINPLOT(STRUCT), VIOLINPLOT(DATASET)
%   plots violins for each column in TABLE, each field in STRUCT, and
%   each variable in DATASET. The violins are labeled according to
%   the table/dataset variable name or the struct field name.
%
%   VIOLINPLOT(DATAMATRIX, CATEGORYNAMES) plots violins for each
%   column in DATAMATRIX and labels them according to the names in the
%   cell-of-strings CATEGORYNAMES.
%
%   VIOLINPLOT(DATA, CATEGORIES) where double vector DATA and vector
%   CATEGORIES are of equal length; plots violins for each category in
%   DATA.
%
%   violins = VIOLINPLOT(...) returns an object array of
%   <a href="matlab:help('Violin')">Violin</a> objects.
%
%   VIOLINPLOT(..., 'PARAM1', val1, 'PARAM2', val2, ...)
%   specifies optional name/value pairs for all violins:
%     'Width'        Width of the violin in axis space.
%                    Defaults to 0.3
%     'Bandwidth'    Bandwidth of the kernel density estimate.
%                    Should be between 10% and 40% of the data range.
%     'ViolinColor'  Fill color of the violin area and data points.
%                    Defaults to the next default color cycle.
%     'ViolinAlpha'  Transparency of the violin area and data points.
%                    Defaults to 0.3.
%     'EdgeColor'    Color of the violin area outline.
%                    Defaults to [0.5 0.5 0.5]
%     'BoxColor'     Color of the box, whiskers, and the outlines of
%                    the median point and the notch indicators.
%                    Defaults to [0.5 0.5 0.5]
%     'MedianColor'  Fill color of the median and notch indicators.
%                    Defaults to [1 1 1]
%     'ShowData'     Whether to show data points.
%                    Defaults to true
%     'ShowNotches'  Whether to show notch indicators.
%                    Defaults to false
%     'ShowMean'     Whether to show mean indicator
%                    Defaults to false
%     'GroupOrder'   Cell of category names in order to be plotted.
%                    Defaults to alphabetical ordering
%
% Copyright (c) 2016, Bastian Bechtold
% This code is released under the terms of the BSD 3-clause license
%
% -----------
% HALF VIOLIN
% edited from original function 'violinplot', to support half violin plots.
% run comparison to original function to verify changes. 
% currently ONLY tested using STRUCTURE data input. half violins are dependent
% on the order of the numeric fields in the structure. non-numeric fields
% will be ignored.
%
%       data.(field1) data.(field2)     data.(field3) data.(field4)   etc
% 
%   add field 'xlabels' to your input structure to specify the xlabels to
%   display below your half violins! you can either have one entry per
%   group of half violins, or one entry per half violin.
%       data.xlabels = {'timepoint1' 'timepoint2'}
%       data.xlabels = {'time1group1', 'time1group2', 'time2group1', ...
%       'time1group2'};
%
%   add field 'legendlabels' to your input structure to specify labels for
%   a legend. only 2 legend fields are supported - one for each side of the
%   half violin plots. if labels names are not specified, a legend wil not
%   be plotted. 
%       data.legendlabels = {'group1', 'group2'}
%
%   added additional name/pair options as detailed below
%     'BandwidthScale'  Set the scale of the bandwidth for the kernel
%                       density estimate. Should be between 0.1 and 0.4 
%                       (10%-40%). The datarange will be automatically
%                       calculated for each data field. 
%                       Bandwidth = (bandwidthscale)*range(data.(field))
%
%     'ViolinColorMat'  Pass in a matrix of RGB color values. Each row will
%                       correspond to the color of the violin. 
%                       # rows should equal # data fields
%
% ALP 11/9/21

    hascategories = exist('cats','var') && not(isempty(cats));
    
    %parse the optional grouporder argument 
    %if it exists parse the categories order 
    % but also delete it from the arguments passed to Violin
    grouporder = {};
    idx=find(strcmp(varargin, 'GroupOrder'));
    if ~isempty(idx) && numel(varargin)>idx
        if iscell(varargin{idx+1})
            grouporder = varargin{idx+1};
            varargin(idx:idx+1)=[];
        else
            error('Second argument of ''GroupOrder'' optional arg must be a cell of category names')
        end
    end
    
    %if bandwidth scale
    idx = find(strcmp(varargin, 'BandWidthScale')); 
    if ~isempty(idx) && numel(varargin) > idx
        scaleval = varargin{idx+1}; 
        varargin(idx:idx+1) = [];
    end
    
    %if input color
    idx = find(strcmp(varargin, 'ViolinColorMat')); 
    if ~isempty(idx) && numel(varargin) > idx
        colors = varargin{idx+1}; 
        varargin(idx:idx+1) = [];
    end

    % tabular data
    if isa(data, 'dataset') || isstruct(data) || istable(data)
        if isa(data, 'dataset')
            colnames = data.Properties.VarNames;
        elseif istable(data)
            colnames = data.Properties.VariableNames;
        elseif isstruct(data)
            colnames = fieldnames(data);
        end
        catnames = {};
        for n=1:length(colnames)
            if isnumeric(data.(colnames{n}))
                catnames = [catnames colnames{n}];
            end
        end
        for n=1:length(catnames)
            thisData = data.(catnames{n});
            if exist('scaleval', 'var') %ALP 2/3/22 calculate the bandwidth if desired
                bandwidth = scaleval*range(thisData); 
                varargin = [varargin {'Bandwidth'}, {bandwidth}];
            end
            violins(n) = Violin_half(thisData, n, varargin{:});
            plotinfo{n} = catnames{n}; %sanity check for you to make sure your data plotted in the correct order
        end
        
        %%% alp 11/9/21 to deal with adding xlabels. you can either specify
        %%% one label for each group of 2 violins, or have a label for each
        %%% half violin
        if isfield(data, 'xlabels')
            xlabels = data.xlabels; 
        else
            xlabels = catnames; 
        end
        if length(xlabels) == length(catnames)/2
            xtickvect = 1:1:length(catnames)/2; 
        else
            xticks = [0.8; 1.2]; %seems like a good distance to prevent text overlap
            xticks = repmat(xticks, [1 size(thisData,2)/2]); 
            addvect = [0:1:size(thisData,2)/2-1; 0:1:size(thisData,2)/2-1];
            xtickvect = xticks+addvect; 
            xtickvect = reshape(xtickvect, [1, size(thisData,2)]);
            
%             xticks = [0.8; 1.2]; %seems like a good distance to prevent text overlap
%             xticks = repmat(xticks, [1 length(catnames)/2]); 
%             addvect = [0:1:length(catnames)/2-1; 0:1:length(catnames)/2-1];
%             xtickvect = xticks+addvect; 
%             xtickvect = reshape(xtickvect, [1, length(catnames)]);
        end
    
        set(gca, 'XTick', xtickvect, 'XTickLabels', xlabels);
        
        % alp 11/9/21 add support for a legend if desired
        if isfield(data, 'legendlabels')
            legend([violins(end-1).ViolinPlot violins(end).ViolinPlot], data.legendlabels{:})
            %legend('Location', 'best')
            legend('boxoff')
        end
        
        %%% alp 2/3/22 add support for changing the color
        if exist('colors', 'var')
            for n = 1:length(violins)
                violins(n).ViolinColor = colors(n,:);
            end
        end

    % 1D data, one category for each data point
    elseif hascategories && numel(data) == numel(cats)

        if isempty(grouporder)
            cats = categorical(cats);
        else
            cats = categorical(cats, grouporder);
        end

        catnames = (unique(cats)); % this ignores categories without any data
        catnames_labels = {};
        for n = 1:length(catnames)
            thisCat = catnames(n);
            catnames_labels{n} = char(thisCat);
            thisData = data(cats == thisCat);
            violins(n) = Violin_half(thisData, n, varargin{:});
        end
        set(gca, 'XTick', 1:1:length(catnames), 'XTickLabels', catnames_labels);

    % 1D data, no categories
    elseif not(hascategories) && isvector(data)
        violins = Violin_half(data, 1, varargin{:});
        set(gca, 'XTick', 1);

    % 2D data with or without categories
    elseif ismatrix(data)
        for n=1:size(data, 2)
            thisData = data(:, n);
            violins(n) = Violin_half(thisData, n, varargin{:});
        end
        set(gca, 'XTick', 1:size(data, 2));
        if hascategories && length(cats) == size(data, 2)
            set(gca, 'XTickLabels', cats);
        end

    end
    
    %alp 11/9/21 cuz prettier i think
    set(gca, 'TickDir', 'out')
end
