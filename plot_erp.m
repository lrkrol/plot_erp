% h = plot_erp(epochs, channel [, varargin])
%
%       Plots ERPs from any number of given epoched EEG sets, for a single
%       channel. Can optionally calculate and plot a difference wave,
%       standard errors, and statistics.
%
% In:
%       epochs - 1-by-n cell containing 1-by-m cells or structs of epoched 
%                EEG set(s) to average and plot. the length of epochs determines 
%                the amount of ERPs plotted, the length of epochs{1} the amount 
%                of sets merged together and averaged for the first ERP curve.
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
%       plotdiff - whether or not to plot the difference, i.e. the first
%                  ERP minus the second one. it will discard all other ERPs
%                  (0|1, default 0).
%       permute - number of permutations per sample for statistics to be
%                 calculated between the first two ERPs. if > 0,
%                 permutation tests will be calculated for every sample and
%                 p-values will be plotted as grey bars behind the plot. 
%                 this will discard all other ERPs. (default 0 disables)
%       labels - cell of legend entries for the ERPs. (default: no legend)
%       colors - array of colors for the ERPs. (default lines(length(epochs)))
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
%       smoothing - whether or not to antialias the curves and use gradients
%                   for the p-values. switch off if you wish to save the
%                   figure in a vector format. (0|1, default 1)
%       newfig - enable to plot figure in new window, switch off to plot
%                in existing (sub)figure. (0|1, default 1)
%
% Out:  
%       h - handle of the generated figure
%
% Usage example:
%       >> plot_erp({{EEGa1, EEGa2}, {EEGb1}, {EEGc1}}, 'Cz', 'labels', ...
%                   {'case a', 'case b', 'case c'}, 'colors', winter(3))
%       >> plot_erp({ALLEEG(1:8:145), ALLEEG(8:8:152)}, 'Fz', ...
%                   'plotdiff', 1, 'plotstd', 'fill')
% 
%                       Laurens R Krol, 2015
%                       Team PhyPA, Biological Psychology and Neuroergonomics,
%                       Berlin Institute of Technology
%                       CC BY-SA 4.0

% 2017-04-04 lrk
%   - Merged plotDiffERP into plotMultERP to form plot_erp,
%     i.e. added difference wave option and statistics
%   - Added standard error also to the difference wave
% 2017-03-13 lrk
%   - Added function to plot standard errors of the mean (plotstd)
%   - Script now accepts cells of structs as well as of cells
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

function h = plot_erp(epochs, channel, varargin)

% parsing input
p = inputParser;

addRequired(p, 'epochs', @iscell);
addRequired(p, 'channel', @ischar);

addParamValue(p, 'avgmode', 'auto', @(x) any(validatestring(x,{'auto', 'within', 'across'})));
addParamValue(p, 'plotstd', 'none', @(x) any(validatestring(x,{'fill', 'lines', 'none'})));
addParamValue(p, 'plotdiff', 0, @isnumeric);
addParamValue(p, 'permute', 0, @isnumeric);
addParamValue(p, 'labels', {}, @(x) (length(x) == length(epochs)));
addParamValue(p, 'colors', lines(length(epochs)), @isnumeric);
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
plotdiff = p.Results.plotdiff;
permute = p.Results.permute;
plotstd = p.Results.plotstd;
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

% transforming epochs to cell in case of struct
for i = 1:length(epochs)
    if isstruct(epochs{i})
        newepochs = cell(1, length(epochs{i}));
        for ii = 1:length(epochs{i})
            newepochs{ii} = epochs{i}(ii);
        end
        epochs{i} = newepochs;
    end
end

% setting average mode
if strcmp(avgmode, 'auto')
    if length(epochs{1}) > 1, avgmode = 'within';
    else avgmode = 'across'; end
end

% setting default labels
if isempty(labels)
    for n = 1:length(epochs)
        labels = [labels, {sprintf('condition %d', n)}];
    end
end

% getting ERP data
erps = [];
stderrs = [];
vars = [];
numepochs = [];
for n = 1:length(epochs)
    % looping through separate ERPs
    conditionerp = [];
    conditionnumepochs = 0;
    for m = 1:length(epochs{n})
        % looping through all datasets making up the individual ERPs
        % getting channel data
        channelidx = find(ismember({epochs{n}{m}.chanlocs.labels}, channel));
        if isempty(channelidx)
            warning('Channel %s not found in epochs{%d}{%d} (%s)', channel, n, m, epochs{n}{m}.setname);
        else 
            if strcmp(avgmode, 'within')
                % taking average of each dataset, or
                conditionerp = [conditionerp; mean(epochs{n}{m}.data(channelidx,:,:), 3)];
                conditionnumepochs = conditionnumepochs + 1;
            else
                % appending all epochs together
                conditionerp = [conditionerp; squeeze(epochs{n}{m}.data(channelidx,:,:))'];
                conditionnumepochs = conditionnumepochs + size(epochs{n}{m}.data, 3);
            end
        end
    end
    
    % saving non-averaged data of the first two ERPs for permutation tests
    if permute > 0 && n == 1
        erp1permute = conditionerp;
    elseif permute > 0 && n == 2
        erp2permute = conditionerp;
    end
    
    % averaging each ERP, saving sample size
    erps = [erps; mean(conditionerp, 1)];
    numepochs = [numepochs, conditionnumepochs];
    
    % getting standard errors
    if size(conditionerp, 1) == 1
        warning('plotstd: Cannot calculate standard error on single sample');
        stderrs = [stderrs; zeros(1, length(conditionerp))];
    else
        stderrs = [stderrs; std(conditionerp, 1) / sqrt(conditionnumepochs)];
        if plotdiff
            % standard error of the difference of the mean uses variance instead
            vars = [vars; var(conditionerp, 1)];
        end
    end
end

if any(all(isnan(erps), 2)), error('Could not generate all ERPs'); end

% calculating sample-by-sample statistics
pvals = [];
if permute > 0
    if size(erps, 1) > 2,
        warning('permute: Taking only first two ERPs, discarding the rest');
        erps = erps(1:2, :);
        stderrs = stderrs(1:2, :);
    end    
    w = waitbar(0, '', 'Name', 'Plot ERP');
    for i = 1:size(erp1permute, 2)
        w = waitbar(i / size(erp1permute, 2), w, sprintf('Statistical testing: sample %d of %d', i, size(erp1permute, 2)));
        pvals = [pvals, permutationTest(erp1permute(:,i), erp2permute(:,i), permute)];
    end
    delete(w);
end

% calculating difference curve, checking/changing relevant settings
if plotdiff
    if size(erps, 1) < 2
        warning('plotdiff: Cannot generate difference for single ERP');
    elseif size(erps, 1) > 2,
        warning('plotdiff: Taking only first two ERPs, discarding the rest');
        erps = erps(1:2, :);
        stderrs = stderrs(1:2, :);
        vars = vars(1:2, :);
    end
    
    if size(colors, 1) ~= 3
        warning('plotdiff: Taking default colour scheme because size(colors, 1) is not 3');
        colors = [0 51 153; 51 153 102; 255 102 51] ./ 255;
    end
    
    if numel(labels) ~= 2
        warning('plotdiff: Taking default labels because numel(labels) is not 2');
        labels = {'condition 1', 'condition 2', 'difference (1-2)'};
    else
        labels = [labels, {'difference'}];
    end
    
    % adding difference to ERP matrix
    erps = [erps; erps(1,:) - erps(2,:)];
    stderrs = [stderrs; sqrt(vars(1,:) / numepochs(1) + vars(2,:) / numepochs(2))];
end

% getting x axis limits, applying delay correction
xmin = epochs{1}{1}.xmin - delaycorrection;
xmax = epochs{1}{1}.xmax - delaycorrection;

% getting x axis indices
x = xmin:1/epochs{1}{1}.srate:xmax;

% drawing figure
if newfig, h = figure('units', 'pixels', 'Position', figpos, 'Color', 'w');
else h = NaN; end

if smoothing, lsmoothing = 'on';
else lsmoothing = 'off'; end

hold on;

% plotting ERP curves
for i = 1:size(erps, 1)
    curves(i) = plot(x, erps(i,:), 'Color', colors(i,:), 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
end

% plotting standard error curves
stdfills = [];
stdlines1 = [];
stdlines2 = [];
if strcmp(plotstd, 'fill')
    fillx = [x, fliplr(x)];
    for i = 1:size(stderrs, 1)
        filly = [erps(i,:) + stderrs(i,:), fliplr(erps(i,:) - stderrs(i,:))];
        stdfills(i) = fill(fillx, filly, colors(i,:));
        if smoothing, set(stdfills(i), 'facealpha', .25); end
        set(stdfills(i), 'edgecolor', 'none');
    end
elseif strcmp(plotstd, 'lines')
    for i = 1:size(stderrs, 1)
        alpha = .75;
        c = colors(i,:);
        c = [c(1) + alpha * (1 - c(1)), c(2) + alpha * (1 - c(2)), c(3) + alpha * (1 - c(3))];
        stdlines1(i) = plot(x, erps(i,:) + stderrs(i,:), 'Color', c, 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
        stdlines2(i) = plot(x, erps(i,:) - stderrs(i,:), 'Color', c, 'LineSmoothing', lsmoothing, 'LineWidth', linewidth);
    end
end

% getting data ranges
if strcmp(plotstd, 'none')
    ymax = max(max(erps));
    ymin = min(min(erps));
else
    ymax = max(max(erps + stderrs));
    ymin = min(min(erps - stderrs));
end
yrange = ymax - ymin;
xrange = xmax - xmin;

% setting plot limits
ylim([-max(abs([ymax ymin])), max(abs([ymax ymin]))])
xlim([xmin xmax]);

% increasing x limits slightly to give more drawing room
xmargin = xrange / 10;
xlim([xmin - xmargin, xmax + xmargin]);

% getting xticks
if xticksize == 0
    % from figure, or
    xticks = get(gca, 'XTick');
else
    % generating custom xticks, spaced around 0
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

% setting length of the tick lines
xticklinelength = yrange / 50;
yticklinelength = (xlim / ylim) * xticklinelength;

% drawing the y-axis
if yticksize == 0
    yticksize = yrange/4;
    if yticksize > 1, yticksize = round(yticksize);
    else yticksize = roundn(yticksize, -1); end
end
yaxisx = [-yticklinelength,  yticklinelength,  0        ,  0         ,  -yticklinelength,  yticklinelength];
yaxisy = [yticksize       ,  yticksize      ,  yticksize,  -yticksize,  -yticksize      ,  -yticksize     ];
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
        xaxisx = [xaxisx,  xticks(t),  xticks(t)      ,  xticks(t)       ];
        xaxisy = [xaxisy,  0        ,  xticklinelength,  -xticklinelength];
    else
        xaxisx = [xaxisx,  xticks(t),  xticks(t)      ,  xticks(t)       ,  xticks(t),  xticks(t+1)];
        xaxisy = [xaxisy,  0        ,  xticklinelength,  -xticklinelength,  0        ,  0          ];
    end
end
xaxis = line(xaxisx, xaxisy, 'Color', [0 0 0]);

% drawing x-axis scale and label
if xscalepos > 0
    xscaley = [-yticksize+xticklinelength, -yticksize-xticklinelength, -yticksize, -yticksize, -yticksize+xticklinelength, -yticksize-xticklinelength];
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
            legendx = xmin;
            align = 'left';
        case 2
            legendx = 0;
            align = 'right';
        case 3
            legendx = 0;
            align = 'center';
        case 4
            legendx = 0;
            align = 'left';
        case 5
            legendx = xmax;
            align = 'right';
    end
    
    % making distinction between regular ERPs (which have a sample size) and the difference wave (which does not)
    if plotdiff, numnondiffcurves = size(erps, 1) - 1;
    else numnondiffcurves = size(erps, 1); end
    
    legendtext = [];
    for i = 1:numnondiffcurves
        legendtext = [legendtext, '\color[rgb]{' num2str(colors(i,:)) '}' labels{i} ' (n=' num2str(numepochs(i)) ')'];
        if i < numnondiffcurves
            legendtext = [legendtext char(10)]; % line break
        end
    end
    if plotdiff
        legendtext = [legendtext, char(10) '\color[rgb]{' num2str(colors(size(erps, 1),:)) '}' labels{size(erps, 1)}];
    end
    legend = text(double(legendx), double(yticksize * -1.35 * vscale^1.35), legendtext, 'HorizontalAlignment', align, 'VerticalAlignment', 'top');
end

% visualising p-values
pbars = [];
if permute > 0
    pbarwidth = 1/epochs{1}{1}.srate;
    for i = 1:length(pvals)
        if pvals(i) >= 0.05, continue;
        else color = maptorange(pvals(i), [.05 0], [1 .75], 'exp', 2); end
        
        patchx = [x(i)-pbarwidth/2, x(i)+pbarwidth/2, x(i)+pbarwidth/2, x(i)-pbarwidth/2];

        if smoothing
            % drawing patch objects with gradients
            patchy = [0 0 yticksize yticksize];
            patchc(1,1,:) = repmat(color,1,3);
            patchc(1,2,:) = repmat(color,1,3);
            patchc(1,3,:) = ones(1,3);
            patchc(1,4,:) = ones(1,3);
            pbars = [pbars, patch(patchx, patchy, patchc, 'EdgeColor', 'none'), patch(patchx, patchy*-1, patchc, 'EdgeColor', 'none')];
        else
            % drawing solid rectangles
            pbars = [pbars, rectangle('Position', [x(i)-pbarwidth/2, -yticksize/2, pbarwidth, yticksize], 'EdgeColor', 'none', 'FaceColor', repmat(color, 1, 3))];
        end
    end
end

% setting figure drawing order, removing original axes, scaling figure to
% fill window, setting font size
set(gca, 'Children', [curves, xaxis, yaxis, ylabelp, ylabeln, xscale, xlabel, chanlabel, legend, stdfills, stdlines1, stdlines2, pbars]);
set(gca, 'Visible', 'off');
set(findall(gcf,'type','text'), 'FontSize', fontsize);
if newfig, set(gca, 'Position', [0 .05 1 .90]); end

end