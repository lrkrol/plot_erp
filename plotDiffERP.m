% h = plotDiffERP(epochs1, epochs2, channel [, name1, name2, delaycorrection, yticksize, vscale, xscalepos, legendpos, smoothing])
%
%       Plots two ERPs from epoched EEG data, and their difference wave,
%       for a single channel.
%
% In:
%       epochs1 - 1-by-n cell of epoched EEG set(s) to average and plot.
%                 multiple sets are first averaged within, then across sets.
%       epochs2 - 1-by-m cell of epoched EEG set(s) to average, plot, and
%                 subtract from epochs1; the difference will also be plotted.
%       channel - index of channel to plot
%
% Optional:
%       name1 - legend entry for the first curve (default 'condition 1')
%       name2 - legend entry for the second curve (default 'condition 2')
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
%       >> plotDiffERP({ALLEEG(1), ALLEEG(3)}, {ALLEEG(2), ALLEEG(4)}, 14)
% 
%                       Laurens R Krol, 2015
%                       Team PhyPA, Department of Biological Psychology and
%                       Neuroergonomics, Berlin Institute of Technology

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

function h = plotDiffERP(epochs1, epochs2, channel, name1, name2, delaycorrection, yticksize, vscale, xscalepos, legendpos, smoothing)

% setting defaults
if (~exist('delaycorrection', 'var')); delaycorrection = 0; end
if (~exist('name1', 'var')); name1 = 'condition 1'; end
if (~exist('name2', 'var')); name2 = 'condition 2'; end
if (~exist('yticksize', 'var')); yticksize = 0; end
if (~exist('vscale', 'var')); vscale = 1; end
if (~exist('xscalepos', 'var')); xscalepos = 7; end
if (~exist('legendpos', 'var')); legendpos = 4; end
if (~exist('smoothing', 'var')); smoothing = 'on'; end

% averaging over epochs
erp1 = [];
erp1n = 0;
for n = 1:length(epochs1)
    erp1 = [erp1; mean(epochs1{n}.data(channel,:,:), 3)];
    erp1n = erp1n + size(epochs1{n}.data, 3);
end
erp1 = mean(erp1, 1);

erp2 = [];
erp2n = 0;
for n = 1:length(epochs2)
    erp2 = [erp2; mean(epochs2{n}.data(channel,:,:), 3)];
    erp2n = erp2n + size(epochs2{n}.data, 3);
end
erp2 = mean(erp2, 1);

diff = erp1 - erp2;

% applying delay correction
xmin = epochs1{1}.xmin - delaycorrection;
xmax = epochs1{1}.xmax - delaycorrection;

% getting x values
x = xmin:1/epochs1{1}.srate:xmax;

% drawing figure
h = figure('Color', [1 1 1]);
colors = [4/255 129/255 255/255; 53/255 181/255 0/255; 255/255 4/255 141/255];
curve1 = plot(x, erp1, 'Color', colors(1,:), 'LineSmoothing', smoothing);
hold on;
curve2 = plot(x, erp2, 'Color', colors(2,:), 'LineSmoothing', smoothing);
curved = plot(x, diff, 'Color', colors(3,:), 'LineSmoothing', smoothing);

% getting data ranges
ymax = max([erp1, erp2, diff]);
ymin = min([erp1, erp2, diff]);
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
    xposticks = xticks(xticks < 0.0001);    % taking rounding errors into account
    xpos1 = xposticks(end-1);
    xpos2 = xposticks(end);
    if xscalepos == 2
        xaxislabely = xaxislabely * -1;
        textposy = textposy * -1;
    end
elseif xscalepos == 3 || xscalepos == 7
    xposticks = xticks(xticks > 0.0001);
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
chanlabel = text(0, double(yticksize * 1.4 * vscale^1.35), ['\bf' epochs1{1}.chanlocs(channel).labels], 'HorizontalAlignment', 'center');

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

    legend = text(double(x), double(yticksize * -1.35 * vscale^1.35), ['\color[rgb]{' num2str(colors(1,:)) '}' name1 ' (n=' num2str(erp1n) ')' char(10) '\color[rgb]{' num2str(colors(2,:)) '}' name2 ' (n=' num2str(erp2n) ')' char(10) '\color[rgb]{' num2str(colors(3,:)) '}difference'] , 'Color', colors(1,:), 'HorizontalAlignment', align, 'VerticalAlignment', 'top');
end

% setting figure drawing order, removing original axes, scaling figure to fill window
set(gca, 'Children', [curve1, curve2, curved, xaxis, yaxis, labelyp, labelyn, xaxislabel, labelx, chanlabel, legend]);
set(gca, 'Visible', 'off');
set(gca, 'Position', [0 .05 1 .90]);

end