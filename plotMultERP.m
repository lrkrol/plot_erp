% h = plotMultERP(epochs, channel [, varargin])
%
%       Plots ERPs from any number of given epoched EEG sets, for a single
%       channel.
%
% In:
%       epochs - 1-by-n cell containing 1-by-m cells of epoched EEG set(s)
%                to average and plot. the length of epochs determines the
%                amount of ERPs plotted, the length of epochs{1} the amount of
%                sets merged together and averaged for the first ERP curve.
%       channel - label of channel to plot: note that EEG sets must have
%                 labels for the channel positions.
%
% Optional (name-value pairs):
%       avgmode - 'within', 'across', or 'auto'. whether to first average
%                 within the given datasets before averaging them together,
%                 or to first append all epochs across datasets and then take
%                 their mean. default 'auto' selects 'within' when when
%                 more than one dataset is provided for epochs{1}, 'across'
%                 otherwise.
%       labels - cell of legend entries for the ERPs. (default: no legend)
%       colors - array of colors for the ERPs. (default hsv(length(epochs)))
%       delaycorrection - delay to visually correct for in seconds, i.e.
%                         where the zero point should be on the x axis
%       yticksize - scale of the y-axis, in microvolts. (default 0 attempts
%                   to find some value automatically)
%       xticksize - scale of the x-axis, in seconds. (default 0 uses
%                   MATLAB's default plot ticks)
%       vscale - scaling factor for the vertical spacing of elements.
%                vertical spacing is relative to yticksize; adjusting the
%                vscale can be useful when the default settings don't work
%                out too well. (default 1)
%       xscalepos - position of the x-axis scale indicator (0 disables,
%                   default 7):
%                   1   2 | 3   4
%                   ------+------
%                   5   6 | 7   8
%       legendpos - position of the legend (0 disables, default 4):
%                         | 
%                   ------+------
%                         |
%                   1    234    5
%       figpos - [x y w h] position and size of the figure on the screen,
%                in pixels. (default [300 400 600 425])
%       fontsize - font size of all text in the figure. (default 10)
%       linewidth - line width for all graphs in the figure. (default 1)
%       smoothing - whether or not to antialias the curves. switch off
%                   if you wish to save the figure in a vector format.
%                   (0|1, default 1)
%       newfig - enable to plot figure in new window, switch off to plot
%                in existing (sub)figure. (0|1, default 1)
%
% Out:  
%       h - handle of the generated figure
%
% Usage example:
%       >> plotMultERP({{EEGa1, EEGa2}, {EEGb1}, {EEGc1}}, 'Cz', 'labels',
%          {'case a', 'case b', 'case c'}, 'colors', winter(3))
% 
%                       Laurens R Krol, 2015
%                       Team PhyPA, Biological Psychology and Neuroergonomics,
%                       Berlin Institute of Technology

% 2017-01-13 lrk
%   - Switched to inputParser to parse arguments
%   - Channel is now given as label, not as index, for consistency when
%     channel is missing in some datasets
%   - X-axis scale indicator can now be disabled
%   - Fixed a placement bug when xscalepos was 2 or 6
%   - Changed automatic y scaling to round to nearest decimal when below 1
%   - Added avgmode, figpos, fontsize, linewidth, xticksize, newfig functions/
%     parameters
% 2016-11-07 lrk
%   - Added number of averaged epochs (n=...) to legend
%   - Changed figure scaling to better fit window
%   - Fixed xscalepos positioning bug due to rounding error
% 2016-01-21 lrk
%   - Improved automatic y scaling for values below 1
%   - Changed the y scale label format to %+1.1d instead of %+d
%   - Added a brief description of the file
% 2015-11-06 First version

function h = plotMultERP(epochs, channel, varargin)

% parsing input
p = inputParser;

addRequired(p, 'epochs', @iscell);
addRequired(p, 'channel', @ischar);

addParamValue(p, 'avgmode', 'auto', @(x) any(validatestring(x,{'auto', 'within', 'across'})));
addParamValue(p, 'labels', {}, @(x) (length(x) == length(epochs)));
addParamValue(p, 'colors', hsv(length(epochs)), @isnumeric);
addParamValue(p, 'delaycorrection', 0, @isnumeric);
addParamValue(p, 'yticksize', 0, @isnumeric);
addParamValue(p, 'xticksize', 0, @isnumeric);
addParamValue(p, 'vscale', 1, @isnumeric);
addParamValue(p, 'xscalepos', 7, @isnumeric);
addParamValue(p, 'legendpos', 4, @isnumeric);
addParamValue(p, 'figpos', [300 400 600 425], @isnumeric);
addParamValue(p, 'fontsize', 10, @isnumeric);
addParamValue(p, 'linewidth', 1, @isnumeric);
addParamValue(p, 'smoothing', 1, @isnumeric);
addParamValue(p, 'newfig', 1, @isnumeric);

parse(p, epochs, channel, varargin{:})

avgmode = p.Results.avgmode;
labels = p.Results.labels;
colors = p.Results.colors;
delaycorrection = p.Results.delaycorrection;
yticksize = p.Results.yticksize;
xticksize = p.Results.xticksize;
vscale = p.Results.vscale;
xscalepos = p.Results.xscalepos;
legendpos = p.Results.legendpos;
smoothing = p.Results.smoothing;
figpos = p.Results.figpos;
fontsize = p.Results.fontsize;
linewidth = p.Results.linewidth;
newfig = p.Results.newfig;

if strcmp(avgmode, 'auto')
    if length(epochs{1}) > 1, avgmode = 'within';
    else avgmode = 'across'; end
end

if isempty(labels)
    for n = 1:length(epochs)
        labels = [labels, {sprintf('condition %d', n)}];
    end
end

% getting ERP data
erps = [];
numepochs = [];
for n = 1:length(epochs)
    conditionerp = [];
    conditionnumepochs = 0;
    for m = 1:length(epochs{n})
        channelidx = find(ismember({epochs{n}{m}.chanlocs.labels}, channel));
        if isempty(channelidx)
            warning('Channel %s not found in epochs{%d}{%d} (%s)', channel, n, m, epochs{n}{m}.setname);
        else 
            if strcmp(avgmode, 'within')
                conditionerp = [conditionerp; mean(epochs{n}{m}.data(channelidx,:,:), 3)];
                conditionnumepochs = conditionnumepochs + 1;
            else
                conditionerp = [conditionerp; squeeze(epochs{n}{m}.data(channelidx,:,:))'];
                conditionnumepochs = conditionnumepochs + size(epochs{n}{m}.data, 3);
            end
        end
    end
    erps = [erps; mean(conditionerp, 1)];
    numepochs = [numepochs, conditionnumepochs];
end

if any(all(isnan(erps), 2)), error('Could not generate all ERPs'); end

% getting x axis limits, applying delay correction
xmin = epochs{1}{1}.xmin - delaycorrection;
xmax = epochs{1}{1}.xmax - delaycorrection;

% getting x axis indices
x = xmin:1/epochs{1}{1}.srate:xmax;

% drawing figure
if newfig
    h = figure('units', 'pixels', 'Position', figpos, 'Color', 'w');
else
    h = NaN;
end
hold on;
if smoothing, lsmoothing = 'on';
else lsmoothing = 'off'; end
for i = 1:size(erps, 1)
    curves(i) = plot(x, erps(i,:), 'Color', colors(i,:), 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
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
if xticksize == 0
    xticks = get(gca, 'XTick');
else
    % spacing xticks around 0
    xticks = [0];
    
    i = 1;
    while ~(0 - i * xticksize < xmin)
        xticks = [xticks, -i * xticksize];
        i=i+1;
    end
    
    i = 1;
    while ~(0 + i * xticksize > xmax)
        xticks = [xticks, i * xticksize];
        i=i+1;
    end
    
    xticks = sort(xticks);
end
xmargin = xrange / 10;
xlim([xmin - xmargin, xmax + xmargin]);

% setting width of the tick lines
xtickwidth = yrange / 50;
ytickwidth = (xlim / ylim) * xtickwidth;

% drawing the y-axis
if yticksize == 0
    yticksize = yrange/4;
    if yticksize > 1, yticksize = round(yticksize);
    else yticksize = roundn(yticksize, -1); end
end
yaxisx = [-ytickwidth ytickwidth 0         0          -ytickwidth ytickwidth];
yaxisy = [yticksize   yticksize  yticksize -yticksize -yticksize  -yticksize];
yaxis = line(yaxisx, yaxisy, 'Color', [0 0 0]);

% drawing y-axis labels
if yticksize < 1, ylabelformat = '%+1.1f';
else ylabelformat = '%+d'; end
ylabelp = text(0, double(yticksize * 1.25 * vscale), [num2str(yticksize, ylabelformat) '{\mu}V'], 'VerticalAlignment', 'cap', 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
ylabeln = text(0, double(yticksize * -1.25 * vscale), [num2str(-yticksize, ylabelformat) '{\mu}V'], 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center', 'Color', [0 0 0]);

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
if xscalepos > 0
    xscaley = [-yticksize+xtickwidth, -yticksize-xtickwidth, -yticksize, -yticksize, -yticksize+xtickwidth, -yticksize-xtickwidth];
    xlabely = yticksize * -1.25 * vscale;
    if xscalepos == 1 || xscalepos == 5
        xpos1 = xticks(1);
        xpos2 = xticks(2);
        if xscalepos == 1
            xscaley = xscaley * -1;
            xlabely = xlabely * -1;
        end
    elseif xscalepos == 2 || xscalepos == 6
        xtempticks = xticks(xticks < 0.0001);    % taking rounding errors into account
        xpos1 = xtempticks(end-2);
        xpos2 = xtempticks(end-1);
        if xscalepos == 2
            xscaley = xscaley * -1;
            xlabely = xlabely * -1;
        end
    elseif xscalepos == 3 || xscalepos == 7
        xtempticks = xticks(xticks > 0.0001);
        xpos1 = xtempticks(1);
        xpos2 = xtempticks(2);
        if xscalepos == 3
            xscaley = xscaley * -1;
            xlabely = xlabely * -1;
        end
    elseif xscalepos == 4 || xscalepos == 8
        xpos1 = xticks(end-1);
        xpos2 = xticks(end);
        if xscalepos == 4
            xscaley = xscaley * -1;
            xlabely = xlabely * -1;
        end
    end
    xsize = xpos2 - xpos1;
    xlabelx = xpos1 + xsize / 2;
    xscalex = [xpos1 xpos1 xpos1 xpos2 xpos2 xpos2];
    xscale = line(xscalex, xscaley, 'Color', [0 0 0]);
    xlabel = text(double(xlabelx), double(xlabely), [num2str(xsize) 's'], 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
else
    xscale = line([0 0], [0 0]);
    xlabel = text();
end

% drawing channel name
chanlabel = text(0, double(yticksize * 1.4 * vscale^1.35), ['\bf' channel], 'HorizontalAlignment', 'center', 'Color', [0 0 0]);

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
        legendtext = [legendtext, '\color[rgb]{' num2str(colors(i,:)) '}' labels{i} ' (n=' num2str(numepochs(i)) ')'];
        if i < size(erps, 1)
            legendtext = [legendtext char(10)];
        end
    end
    legend = text(double(x), double(yticksize * -1.35 * vscale^1.35), legendtext, 'HorizontalAlignment', align, 'VerticalAlignment', 'top');
end

% setting figure drawing order, removing original axes, scaling figure to
% fill window, setting font size
set(gca, 'Children', [curves, xaxis, yaxis, ylabelp, ylabeln, xscale, xlabel, chanlabel, legend]);
set(gca, 'Visible', 'off');
set(findall(gcf,'type','text'), 'FontSize', fontsize);
if newfig, set(gca, 'Position', [0 .05 1 .90]); end

end