function allindex = selectindextable(spreadsheetfile, varargin)

%% get selection info
sheetnumber = 1;                %default is to use first sheet

filters.include = 1;            %default is to include files marked as 1
filters.stimulation = nan;      %default is to select none
filters.behavior = nan;         %default is to select none
filters.experimentday = nan;    %default is to select none
filters.animal = nan;           %default is to select none
filters.datesincluded = nan;    %default is to include all
filters.datesexcluded = nan;    %default is to exclude none

%get input arguments
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'rewritefileinfo'
            rewritefileinfo = varargin{option+1};
        case 'sheetnumber'
            sheetnumber = varargin{option+1};
        case 'include'
            filters.include = varargin{option+1};
        case 'stimulation'
            filters.stimulation = varargin{option+1};
        case 'behavior'
            filters.behavior = varargin{option+1};
        case 'experimentday'
            filters.experimentday = varargin{option+1};
        case 'animal'
            filters.animal = varargin{option+1};
        case 'datesincluded'
            filters.datesincluded = varargin{option+1};
        case 'datesexcluded'
            filters.datesexcluded = varargin{option+1};
        otherwise
            disp(['Warning: ' varargin{option} ' is an invalid input argument'])
    end
end

%% read in the spreadsheet
opts = detectImportOptions(spreadsheetfile);
opts.Sheet = sheetnumber;
ephysT = readtable(spreadsheetfile, opts);
headerNames = ephysT.Properties.VariableNames;

%% apply selected filters
fnames = fieldnames(filters);
filtersApplied = fnames(cellfun(@(x) any(~isnan(filters.(x))), fnames)); % TODO - extract string based filters

% get included rows for each filter
allFilters = zeros(size(ephysT,1),numel(filtersApplied));
for filtIdx = 1:numel(filtersApplied)
    filtName = filtersApplied{filtIdx};        
    if any(strcmp(lower(headerNames), filtName)) %if the filter name applies to a specific column of the table
        filtCol = headerNames{strcmp(lower(headerNames), filtName)};
        allFilters(:,filtIdx) = any(ephysT.(filtCol) == filters.(filtName),2); %gets any matching criteria if multiple filter items
    elseif contains(filtName, 'date') %special case to deal with dates
        if isempty(filters.datesexcluded)
            filters.datesexcluded = nan;
        end
        if isempty(filters.datesincluded)
            filters.datesincluded = nan;
        end
        if isnan(filters.datesincluded) %default is to include all dates
            allFilters(:,filtIdx) = all((ephysT.Date ~= filters.datesexcluded),2);
        else
            allFilters(:,filtIdx) = (any(ephysT.Date == filters.datesincluded, 2) & (ephysT.Date ~= filters.datesexcluded));
        end
    else
        disp(['Warning: no corresponding column found for the ' filtName...
            ' filter. Check the spreadsheet header names'])
    end
end

% apply all filters to the table
allindex = ephysT(all(allFilters,2),:);

end