function savefigALP(datadir, figname, varargin)
%makes figure directory and saves a png and matlab fig
% inputs: datadir, figname
%ALP 4.27.18
%edited ALP 6/30/2022 to save a text file with the name of the functions
%that generated these figures

date = 0; info = 0; ftype = 0;
for option = 1:length(varargin)
    switch varargin{option}
        case 'date'
            date = 1;
        case 'info'
            info = 1;
        case 'filetype'
            ftype = 1;
            filetype = varargin{option+1};
    end
end

if date
    datadir = [datadir,  datestr(now, 'yymmdd'), '\'];
end

if ~exist(datadir)
    mkdir(datadir);
end

set(gcf, 'Visible', 'on')

infoname = ['info_', num2str(datestr(now, 'yymmdd')), '.txt'];
temp = dbstack;
temp2 = arrayfun(@(x) [temp(x).file '/'], 2:size(temp,1), 'UniformOutput', false);
metadata.name = strcat(cell2mat(temp2));
metadata.date = datestr(now);
f = gcf;
f.UserData = metadata;

infostr = ['fig name: ', figname, ' --- generated by: ', strcat(metadata.name), ' on ', num2str(datestr(now))];
fid = fopen([datadir, infoname], 'at');
fprintf(fid, '%s\r\n', strcat(infostr));
fclose(fid);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

filename = [datadir,figname];
set(gcf,'renderer','painters')
if ftype == 0
    saveas(gcf,filename,'png');
    saveas(gcf,filename,'fig');
    print(filename, '-dpdf');
elseif strcmp(filetype, 'pdf')
    print(filename, '-dpdf');
elseif strcmp(filetype, 'png')
    saveas(gcf,filename,'png');
elseif strcmp(filetype, 'fig')
    saveas(gcf,filename,'fig');
end
% print(filename, '-dsvg');

end