% [h, pvals] = plotDiffERP(epochs1, epochs2, channel [, varargin])
%
%       Plots two ERPs from epoched EEG data, and their difference wave,
%       for a single channel.
%
% In:
%       epochs1 - 1-by-n cell or struct of epoched EEG set(s) to average and plot.
%       epochs2 - 1-by-m cell or struct of epoched EEG set(s) to average, plot, and
%                 subtract from epochs1; the difference will also be plotted.
%       channel - label of channel to plot: note that EEG sets must have
%                 labels for the channel positions.
%
% Optional (name-value pairs):
%       avgmode - 'within', 'across', or 'auto'. whether to first average
%                 within the given datasets before averaging them together,
%                 or to first append all epochs across datasets and then take
%                 their mean. default 'auto' selects 'within' when when
%                 more than one dataset is provided for epochs1, 'across'
%                 otherwise.
%       plotstd - whether or not to plot standard errors of the mean for the 
%                 two average curves. 
%                 'fill' - plots a filled, semi-transparent area. note:
%                          this interferes with MATLAB's ability to plot
%                          smooth curves. when smoothing is off, these will
%                          not be semi-transparent to maintain vector
%                          compatibility.
%                 'lines' - plots upper and lower boundary curves. this
%                           option does not use transparency, thus does not
%                           interfere with MATLAB's graphical abilities.
%                 'none' - does not plot standard deviations (default).
%       labels - 1-by-2 cell of chars to be used as legend entries for the
%                two curves. (default {'condition 1', 'condition 2'})
%       colors - array of colors for the ERPs. (default predetermined colors)
%       delaycorrection - delay to visually correct for in seconds, i.e.
%                         where the zero point should be on the x axis
%       yticksize - scale of the y-axis, in microvolts. (default 0 attempts
%                   to find some value automatically)
%       xticksize - scale of the x-axis, in seconds. (default 0 uses
%                   MATLAB's default plot ticks)
%       vscale - scaling factor for the vertical spacing of elements.
%                vertical spacing is relative to yticksize; adjusting the
%                vscale can be useful when using different figure sizes.
%                (default 1)
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
%       permute - number of permutations per sample for statistics. if > 0,
%                 permutation tests will be calculated for every sample and
%                 p-values will be plotted as grey bars behind the plot.
%                 (default 0)
%       smoothing - whether or not to antialias the curves and use gradients
%                   for the p-values. switch off if you wish to save the
%                   figure in a vector format. (0|1, default 1)
%       newfig - enable to plot figure in new window, switch off to plot
%                in existing (sub)figure. (0|1, default 1)
%
% Out:  
%       h - handle of the generated figure
%       pvals - per-sample p-values as calculated using permutation tests,
%               if enabled.
%
% Usage example:
%       >> plotDiffERP(ALLEEG(1:2), {EEG3, EEG4}, 'Cz')
% 
%                       Laurens R Krol, 2015
%                       Team PhyPA, Biological Psychology and Neuroergonomics,
%                       Berlin Institute of Technology

% 2017-03-13 lrk
%   - Added function to plot standard errors of the mean (plotstd)
%   - Added waitbar for statistics
%   - Script now accepts structs as well as cells
% 2017-01-13 lrk
%   - Switched to inputParser to parse arguments
%   - Channel is now given as label, not as index, for consistency when
%     channel is missing in some datasets
%   - Labels are now given as single cell, not as two strings
%   - X-axis scale indicator can now be disabled
%   - Fixed a placement bug when xscalepos was 2 or 6
%   - Changed automatic y scaling to round to nearest decimal when below 1
%   - Added statistical testing and visualisation
%   - Added avgmode, color, figpos, fontsize, linewidth, xticksize, newfig
%     functions/parameters
% 2016-11-07 lrk
%   - Added number of averaged epochs (n=...) to legend
%   - Changed figure scaling to better fit window
%   - Fixed xscalepos positioning bug due to rounding error
% 2016-01-21 lrk
%   - Improved automatic y scaling for values below 1
%   - Changed the y scale label format to %+1.1d instead of %+d
%   - Added a brief description of the file
% 2015-11-06 lrk
%   - Instead of single EEG sets, epoch1 and epoch2 are now cells that may
%     contain more than one set each (e.g. one for each subject)
%   - Fixed a bug where sometimes, text() coordinates were not doubles
% 2015-06-24 lrk
%   - Fixed bug where the channel index for the chanlabel was hardcoded to 5
%   - Fixed vertical alignment of labels that vscale caused to lose alignment
%   - Forced signs on y scale labels
% 2015-05-29 First version

function [h, pvals] = plotDiffERP(epochs1, epochs2, channel, varargin)

% parsing input
p = inputParser;

addRequired(p, 'epochs1', @(x) any([iscell(x), isstruct(x)]));
addRequired(p, 'epochs2', @(x) any([iscell(x), isstruct(x)]));
addRequired(p, 'channel', @ischar);

addParamValue(p, 'avgmode', 'auto', @(x) any(validatestring(x,{'auto', 'within', 'across'})));
addParamValue(p, 'plotstd', 'none', @(x) any(validatestring(x,{'fill', 'lines', 'none'})));
addParamValue(p, 'labels', {'condition 1', 'condition 2'}, @(x) (all(size(x) == [1,2]) && ischar([x{:}])));
addParamValue(p, 'colors', [0 51 153; 51 153 102; 255 102 51] ./ 255, @isnumeric);
addParamValue(p, 'delaycorrection', 0, @isnumeric);
addParamValue(p, 'yticksize', 0, @isnumeric);
addParamValue(p, 'xticksize', 0, @isnumeric);
addParamValue(p, 'vscale', 1, @isnumeric);
addParamValue(p, 'xscalepos', 7, @isnumeric);
addParamValue(p, 'legendpos', 4, @isnumeric);
addParamValue(p, 'figpos', [300 400 600 425], @isnumeric);
addParamValue(p, 'fontsize', 10, @isnumeric);
addParamValue(p, 'linewidth', 1, @isnumeric);
addParamValue(p, 'permute', 0, @isnumeric);
addParamValue(p, 'smoothing', 1, @isnumeric);
addParamValue(p, 'newfig', 1, @isnumeric);

parse(p, epochs1, epochs2, channel, varargin{:})

avgmode = p.Results.avgmode;
plotstd = p.Results.plotstd;
labels = p.Results.labels;
colors = p.Results.colors;
delaycorrection = p.Results.delaycorrection;
yticksize = p.Results.yticksize;
xticksize = p.Results.xticksize;
vscale = p.Results.vscale;
xscalepos = p.Results.xscalepos;
legendpos = p.Results.legendpos;
figpos = p.Results.figpos;
fontsize = p.Results.fontsize;
linewidth = p.Results.linewidth;
permute = p.Results.permute;
smoothing = p.Results.smoothing;
newfig = p.Results.newfig;

% struct to cell
if isstruct(epochs1)
    newepochs1 = cell(1, length(epochs1));
    for i = 1:length(epochs1)
        newepochs1{i} = epochs1(i);
    end
    epochs1 = newepochs1;
end

if isstruct(epochs2)
    newepochs2 = cell(1, length(epochs2));
    for i = 1:length(epochs2)
        newepochs2{i} = epochs2(i);
    end
    epochs2 = newepochs2;
end

if strcmp(avgmode, 'auto')
    if length(epochs1) > 1, avgmode = 'within';
    else avgmode = 'across'; end
end

% getting ERP data
erp1 = [];
erp1n = 0;
for n = 1:length(epochs1)
    channelidx = find(ismember({epochs1{n}.chanlocs.labels}, channel));
    if isempty(channelidx)
        warning('Channel %s not found in epochs1{%d} (%s)', channel, n, epochs1{n}.setname);
    else
        if strcmp(avgmode, 'within')
            erp1 = [erp1; mean(epochs1{n}.data(channelidx,:,:), 3)];
            erp1n = erp1n + 1;
        else
            erp1 = [erp1; squeeze(epochs1{n}.data(channelidx,:,:))'];
            erp1n = erp1n + size(epochs1{n}.data, 3);
        end
    end
end

erp2 = [];
erp2n = 0;
for n = 1:length(epochs2)
    channelidx = find(ismember({epochs2{n}.chanlocs.labels}, channel));
    if isempty(channelidx)
        warning('Channel %s not found in epochs2{%d} (%s)', channel, n, epochs2{n}.setname);
    else
        if strcmp(avgmode, 'within')
            erp2 = [erp2; mean(epochs2{n}.data(channelidx,:,:), 3)];
            erp2n = erp2n + 1;
        else
            erp2 = [erp2; squeeze(epochs2{n}.data(channelidx,:,:))'];
            erp2n = erp2n + size(epochs2{n}.data, 3);
        end
    end
end

% calculating sample-by-sample statistics
pvals = [];
if permute > 0
    w = waitbar(0, 'Calculating statistics');
    for i = 1:size(erp1,2)
        w = waitbar(i / length(erp1));
        pvals = [pvals, permutationTest(erp1(:,i), erp2(:,i), permute)];
    end
    delete(w);
end

% getting standard errors
if size(erp1, 1) == 1
    warning('Cannot calculate standard error on single sample');
    std1 = 0;
else
    std1 = std(erp1, 1) / sqrt(erp1n);
end
if size(erp2, 1) == 1
    warning('Cannot calculate standard error on single sample');
    std2 = 0;
else
    std2 = std(erp2, 1) / sqrt(erp2n);
end

% getting final curve data
erp1 = mean(erp1, 1);
erp2 = mean(erp2, 1);
diff = erp1 - erp2;

if (all(isnan(diff))), error('Could not generate difference wave'); end

% getting x axis limits, applying delay correction
xmin = epochs1{1}.xmin - delaycorrection;
xmax = epochs1{1}.xmax - delaycorrection;

% getting x axis values
xvals = xmin:1/epochs1{1}.srate:xmax;

% drawing figure
if newfig
    h = figure('units', 'pixels', 'Position', figpos, 'Color', 'w');
else
    h = NaN;
end
hold on;
if smoothing, lsmoothing = 'on';
else lsmoothing = 'off'; end

% plotting ERP curves
curve1 = plot(xvals, erp1, 'Color', colors(1,:), 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
curve2 = plot(xvals, erp2, 'Color', colors(2,:), 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
curved = plot(xvals, diff, 'Color', colors(3,:), 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);

% plotting standard error curves
std1fill = [];
std2fill = [];
curve1std1 = [];
curve1std2 = [];
curve2std1 = [];
curve2std2 = [];
if strcmp(plotstd, 'fill')
    fillx = [xvals, fliplr(xvals)];
    filly = [erp1 + std1, fliplr(erp1 - std1)];
    std1fill = fill(fillx, filly, colors(1,:));
    if smoothing, set(std1fill, 'facealpha', .25); end
    set(std1fill, 'edgecolor', 'none');

    filly = [erp2 + std2, fliplr(erp2 - std2)];
    std2fill = fill(fillx, filly, colors(2,:));  
    if smoothing, set(std2fill, 'facealpha', .25); end
    set(std2fill, 'edgecolor', 'none');
elseif strcmp(plotstd, 'lines')
    alpha = .75;
    c1 = colors(1,:);
    c1 = [c1(1) + alpha * (1 - c1(1)), c1(2) + alpha * (1 - c1(2)), c1(3) + alpha * (1 - c1(3))];
    c2 = colors(2,:);
    c2 = [c2(1) + alpha * (1 - c2(1)), c2(2) + alpha * (1 - c2(2)), c2(3) + alpha * (1 - c2(3))];
    curve1std1 = plot(xvals, erp1 + std1, 'Color', c1, 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
    curve1std2 = plot(xvals, erp1 - std1, 'Color', c1, 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
    curve2std1 = plot(xvals, erp2 + std1, 'Color', c2, 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
    curve2std2 = plot(xvals, erp2 - std1, 'Color', c2, 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
end

% getting data ranges
if strcmp(plotstd, 'none')
    ymax = max([erp1, erp2, diff]);
    ymin = min([erp1, erp2, diff]);
else
    ymax = max([erp1 + std1, erp2 + std2, diff]);
    ymin = min([erp1 - std1, erp2 - std2, diff]);
end
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
        i = i + 1;
    end
    
    i = 1;
    while ~(0 + i * xticksize > xmax)
        xticks = [xticks, i * xticksize];
        i = i + 1;
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
if mod(yticksize,1) == 0, ylabelformat = '%+d';
else ylabelformat = '%+1.1f'; end
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

    legend = text(double(x), double(yticksize * -1.35 * vscale^1.35), ['\color[rgb]{' num2str(colors(1,:)) '}' labels{1} ' (n=' num2str(erp1n) ')' char(10) '\color[rgb]{' num2str(colors(2,:)) '}' labels{2} ' (n=' num2str(erp2n) ')' char(10) '\color[rgb]{' num2str(colors(3,:)) '}difference'] , 'Color', colors(1,:), 'HorizontalAlignment', align, 'VerticalAlignment', 'top');
end

% visualising p-values
pbars = [];
if permute > 0
    pbarwidth = 1/epochs1{1}.srate;
    for i = 1:length(pvals)
        if pvals(i) >= 0.05, continue;
        else color = maptorange(pvals(i), [.05 0], [1 .75], 'exp', 2); end
        
        patchx = [xvals(i)-pbarwidth/2, xvals(i)+pbarwidth/2, xvals(i)+pbarwidth/2, xvals(i)-pbarwidth/2];

        if smoothing
            patchy = [0 0 yticksize yticksize];
            patchc(1,1,:) = repmat(color,1,3);
            patchc(1,2,:) = repmat(color,1,3);
            patchc(1,3,:) = ones(1,3);
            patchc(1,4,:) = ones(1,3);
            pbars = [pbars, patch(patchx, patchy, patchc, 'EdgeColor', 'none'), patch(patchx, patchy*-1, patchc, 'EdgeColor', 'none')];
        else
            pbars = [pbars, rectangle('Position', [xvals(i)-pbarwidth/2, -yticksize/2, pbarwidth, yticksize], 'EdgeColor', 'none', 'FaceColor', repmat(color, 1, 3))];
        end
    end
end

% setting figure drawing order, removing original axes, scaling figure to
% fill window, setting font size
set(gca, 'Children', [curve1, curve2, curved, xaxis, yaxis, ylabelp, ylabeln, xscale, xlabel, chanlabel, legend, pbars, std1fill, std2fill, curve1std1, curve1std2, curve2std1, curve2std2]);
set(gca, 'Visible', 'off');
set(findall(gca,'type','text'), 'FontSize', fontsize);
if newfig, set(gca, 'Position', [0 .05 1 .90]); end

end