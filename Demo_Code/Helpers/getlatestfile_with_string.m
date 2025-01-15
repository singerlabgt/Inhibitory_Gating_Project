function latestfile = getlatestfile_with_string(directory,str2locate,varargin)
%This function returns the latest file from the directory passsed as input
%argument
%NJ edited 04.13.20 to find latest file that contains a specific string

%Get the directory contents
dirc = dir(directory);

%Filter out all the folders.
dirc = dirc(find(~cellfun(@isdir,{dirc(:).name})));


if nargin > 2 

    dirc = dirc(contains({dirc.name}, str2locate) & contains({dirc.name}, varargin{1}));
else 
    dirc = dirc(contains({dirc.name}, str2locate) );
end

%I contains the index to the biggest number which is the latest file
[~,I] = max([dirc(:).datenum]);

if ~isempty(I)
    latestfile = dirc(I).name;
else
    latestfile = [];
end

end