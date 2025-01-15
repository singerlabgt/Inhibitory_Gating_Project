function makefigurepretty(fh, varargin)
%makefigurepretty
%
%   inputs:     fh - figure handle
%ALP 2/8/22

figure(fh)
set(fh, 'color', 'w')

ms = 0;
if ~isempty(varargin)
    ms = 1;
end

%%% ----- axes ----
g = findall(fh, 'Type', 'axes');
for i = 1:numel(g)
    % ticks
    g(i).Units = 'inches';
    ticklength = 0.04; 
    pos = g(i).Position;
    g(i).TickDir = 'out';
    normlength = ticklength/max(pos(3:4));
    %     g(i).TickLength = [0.015 0.015];
    g(i).TickLength = [normlength normlength];
    g(i).FontName = 'Arial';
    
    if ms == 0
        g(i).FontSize = 8;
        g(i).LabelFontSizeMultiplier = 1.25;  %10pt
        g(i).Title.FontSize = 12;
        
    elseif ms == 1
        g(i).FontSize = 6;
        g(i).LabelFontSizeMultiplier = 1.25;  %7.5 pt
        g(i).Title.FontSize = 8.5;
    end
    
    axes(g(i))
    box off
    
    
    % move axes away - matlab 2019b
    
end

%%% ---- title ----
t = sgtitle; 
t.FontSize = 14; 

end

