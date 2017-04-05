% h = plotMultERP(epochs, channel [, names, colors, delaycorrection, yticksize, vscale, xscalepos, legendpos, smoothing])
%
%       Plots ERPs from any number of given epoched EEG sets, for a single
%       channel.
%
% In:
%       epochs - 1-by-n cell containing 1-by-m cells of epoched EEG set(s)
%                to average and plot. multiple sets are first averaged within,
%                then across sets. i.e., the length of epochs determines the
%                amount of ERPs plotted, the length of epochs{1} the amount of
%                sets averaged together for the first ERP.
%       channel - index of channel to plot
%
% Optional:
%       names - cell of legend entries for the ERPs (default: no legend)
%       colors - array of colors for the ERPs (default hsv(length(epochs)))
%       delaycorrection - delay to visually correct for in seconds, i.e.
%                         where the zero point should be on the x axis
%       yticksize - scale of the y-axis, in microvolts (default 0 attempts
%                   to find some value automatically)
%       vscale - scaling factor for the vertical spacing of elements.
%                vertical spacing is relative to yticksize; adjusting the
%                vscale can be useful when the default settings don't work
%                out too well. (default 1)
%       xscalepos - position of the x-axis scale indicator (default 7):
%                   1   2 | 3   4
%                   ------+------
%                   5   6 | 7   8
%       legendpos - position of the legend (0 disables, default 4):
%                         | 
%                   ------+------
%                         |
%                   1    234    5
%       smoothing - whether or not to antialias the curves. switch 'off'
%                   if you wish to save the figure in a vector format.
%                   (default 'on')
%
% Out:  
%       h - handle of the generated figure
%
% Usage example:
%       >> plotMultERP({{EEGa1, EEGa2}, {EEGb1}, {EEGc1}}, 14, {'case a',
%          'case b', 'case c'}, winter(3), 0.25)
% 
%                       Laurens R Krol, 2015
%                       Team PhyPA, Department of Biological Psychology and
%                       Neuroergonomics, Berlin Institute of Technology

% 2016-01-21 lrk
%   - Improved automatic y scaling for values below 1
%   - Changed the y scale label format to %+1.1d instead of %+d
%   - Added a brief description of the file
% 2015-11-06 First version

function h = plotMultERP(epochs, channel, names, colors, delaycorrection, yticksize, vscale, xscalepos, legendpos, smoothing)

% setting defaults
if (~exist('colors', 'var')); colors = hsv(length(epochs)); end
if (~exist('delaycorrection', 'var')); delaycorrection = 0; end
if (~exist('yticksize', 'var')); yticksize = 0; end
if (~exist('vscale', 'var')); vscale = 1; end
if (~exist('xscalepos', 'var')); xscalepos = 7; end
if (~exist('legendpos', 'var')); legendpos = 4; end
if (~exist('names', 'var')); legendpos = 0; end
if (~exist('smoothing', 'var')); smoothing = 'on'; end

% getting mean ERPs
erps = [];
for n = 1:length(epochs)
    conditionerp = [];
    for m = 1:length(epochs{n})
        conditionerp = [conditionerp; mean(epochs{n}{m}.data(channel,:,:), 3)];
    end
    erps = [erps; mean(conditionerp, 1)];
end

% applying delay correction
xmin = epochs{1}{1}.xmin - delaycorrection;
xmax = epochs{1}{1}.xmax - delaycorrection;

% getting x values
x = xmin:1/epochs{1}{1}.srate:xmax;

% drawing figure
h = figure('Color', [1 1 1]);
hold on;
for i = 1:size(erps, 1)
    curves(i) = plot(x, erps(i,:), 'Color', colors(i,:), 'LineSmoothing', smoothing);
end

% getting data ranges
ymax = max(max(erps));
ymin = min(min(erps));
yrange = ymax - ymin;
xrange = xmax - xmin;

% setting plot limits
ylim([-max(abs([ymax ymin])), max(abs([ymax ymin]))])
xlim([xmin xmax]);

% getting xticks, then increasing limits slightly to give more drawing room
xticks = get(gca, 'XTick');
xmargin = xrange / 10;
xlim([xmin - xmargin, xmax + xmargin]);

% setting width of the tick lines
xtickwidth = yrange / 50;
ytickwidth = (xlim / ylim) * xtickwidth;

% drawing the y-axis
if yticksize == 0
    if yrange < 2
        yticksize = 10^floor(log10(yrange));
    else
        yticksize = ceil(yrange / 4);
    end
end
yaxisx = [-ytickwidth ytickwidth 0         0          -ytickwidth ytickwidth];
yaxisy = [yticksize   yticksize  yticksize -yticksize -yticksize  -yticksize];
yaxis = line(yaxisx, yaxisy, 'Color', [0 0 0]);

% drawing y-axis labels
labelyp = text(0, double(yticksize * 1.25 * vscale), [num2str(yticksize, '%+1.1d') '{\mu}V'], 'VerticalAlignment', 'cap', 'HorizontalAlignment', 'center');
labelyn = text(0, double(yticksize * -1.25 * vscale), [num2str(-yticksize, '%+1.1d') '{\mu}V'], 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');

% drawing the x-axis
xaxisx = [];
xaxisy = [];
for t = 1:length(xticks);
    if t == length(xticks)
        xaxisx = [xaxisx, xticks(t) xticks(t)   xticks(t)];
        xaxisy = [xaxisy, 0         xtickwidth  -xtickwidth];
    else
        xaxisx = [xaxisx, xticks(t) xticks(t)   xticks(t)    xticks(t) xticks(t+1)];
        xaxisy = [xaxisy, 0         xtickwidth  -xtickwidth  0         0];
    end
end
xaxis = line(xaxisx, xaxisy, 'Color', [0 0 0]);

% drawing x-axis scale and label
xaxislabely = [-yticksize+xtickwidth, -yticksize-xtickwidth, -yticksize, -yticksize, -yticksize+xtickwidth, -yticksize-xtickwidth];
textposy = yticksize * -1.25 * vscale;
if xscalepos == 1 || xscalepos == 5
    xpos1 = xticks(1);
    xpos2 = xticks(2);
    if xscalepos == 1
        xaxislabely = xaxislabely * -1;
        textposy = textposy * -1;
    end
elseif xscalepos == 2 || xscalepos == 6
    xposticks = xticks(xticks < 0);
    xpos1 = xposticks(end-1);
    xpos2 = xposticks(end);
    if xscalepos == 2
        xaxislabely = xaxislabely * -1;
        textposy = textposy * -1;
    end
elseif xscalepos == 3 || xscalepos == 7
    xposticks = xticks(xticks > 0);
    xpos1 = xposticks(1);
    xpos2 = xposticks(2);
    if xscalepos == 3
        xaxislabely = xaxislabely * -1;
        textposy = textposy * -1;
    end
elseif xscalepos == 4 || xscalepos == 8
    xpos1 = xticks(end-1);
    xpos2 = xticks(end);
    if xscalepos == 4
        xaxislabely = xaxislabely * -1;
        textposy = textposy * -1;
    end
end
xsize = xpos2 - xpos1;
textposx = xpos1 + xsize / 2;
xaxislabelx = [xpos1 xpos1 xpos1 xpos2 xpos2 xpos2];
xaxislabel = line(xaxislabelx, xaxislabely, 'Color', [0 0 0]);
labelx = text(double(textposx), double(textposy), [num2str(xsize) 's'], 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');

% drawing channel name
chanlabel = text(0, double(yticksize * 1.4 * vscale^1.35), ['\bf' epochs{1}{1}.chanlocs(channel).labels], 'HorizontalAlignment', 'center');

% drawing legend
if legendpos == 0
    legend = text(0, 0, '');
else
    switch legendpos
        case 1
            x = xmin;
            align = 'left';
        case 2
            x = 0;
            align = 'right';
        case 3
            x = 0;
            align = 'center';
        case 4
            x = 0;
            align = 'left';
        case 5
            x = xmax;
            align = 'right';
    end
    
    legendtext = [];
    for i = 1:size(erps, 1)
        legendtext = [legendtext, '\color[rgb]{' num2str(colors(i,:)) '}' names{i}];
        if i < size(erps, 1)
            legendtext = [legendtext char(10)];
        end
    end
    legend = text(double(x), double(yticksize * -1.35 * vscale^1.35), legendtext, 'HorizontalAlignment', align, 'VerticalAlignment', 'top');
end

% setting figure drawing order, removing original axes
set(gca, 'Children', [curves, xaxis, yaxis, labelyp, labelyn, xaxislabel, labelx, chanlabel, legend]);
set(gca, 'Visible', 'off');

end