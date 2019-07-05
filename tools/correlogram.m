function [bin, count, dt] = correlogram(spike_time, spike_time_ref, range, num_bins, remove_self)

if nargin<5
    remove_self=false;  % remove self for autocorrelogram
end

%% calc diff
[N, d] = size(spike_time);
[Nr, dr] = size(spike_time_ref);

% make spike_time a column vector
if N<d
    spike_time = spike_time';
end

% make spike_time a row vector
if Nr>dr
    spike_time_ref = spike_time_ref';
end

D=bsxfun(@minus, spike_time, spike_time_ref);

if remove_self
    % replace zero in diagonal with nan
    [N, ~] = size(D);
    for n=1:N
        D(n,n) = nan;
    end
end
% D(1:10, 1:10)

%% set bin & window size


% get values  smaller than a max window size
idxValues=D>=range(1) & D<=range(2);

values = D(idxValues);
%D(1:10, 1:10)

%% calc histogram
bin = linspace(range(1),range(2),num_bins);
dt = bin(2)-bin(1);

count = hist(values, bin);

return 


%% Example. autocorreogram for ch_31a (OFF, 25 Hz)
[bin, count, dt] = correlogram(ch_31a, ch_31a, [-1 1], 2001, true);
bar(bin,count, 'k')
xlabel('\Deltat (sec)')
ylabel('count')
box off
title (sprintf('autocoorelogram (\\Deltat=%.3f)',dt))

set(gcf, 'paperposition', [0 0 8 5])
set(gcf, 'papersize', [8 5])
saveas(gcf, sprintf('%scell_%dHz_auto-correlogram.pdf', CELL_TYPE,fps))
saveas(gcf, sprintf('%scell_%dHz_auto-correlogram.png', CELL_TYPE,fps))


%% Example. crosscorreogram for ch_22a and ch_31b
% channelNames = {'ch_12a'}    {'ch_22b'}    {'ch_31a'}    {'ch_32a'}
%ch_idx = 1 % 12a
%ch_idx = 2 % 22b
ch_idx = 4 % 32a
ch_idx_ref = 3 % 31a
[bin, count,dt] = correlogram(eval(channelNames{ch_idx}), eval(channelNames{ch_idx_ref}), [-0.2 0.2], 401);

clf
bar(bin,count, 'k')
xlabel('\Deltat (sec)')
ylabel('count')
box off
title (sprintf('cross-corelogram (\\Deltat=%.3f)',dt))

% % avg_time_diff=bin*count'/sum(count)
% % ylim=get(gca,'ylim')
% % hold on; plot(avg_time_diff*[1 1], ylim, 'k--', 'linewidth',5)

%
% calculate average time difference
range_to_calc_avg=[-0.1 0.1]
idxToCalc= bin>range_to_calc_avg(1) & bin<range_to_calc_avg(2);

bin_chosen = bin(idxToCalc);
count_chosen = count(idxToCalc);
avg_time_diff=bin_chosen*count_chosen'/sum(count_chosen)

hold on
bar(bin_chosen,count_chosen, 'b')
ylim=get(gca,'ylim');
plot(avg_time_diff*[1 1], ylim, 'b--', 'linewidth', 6)
title (sprintf('cross-corelogram %s - %s (\\Deltat=%.3f), avg=%.3f', channelNames{ch_idx}, channelNames{ch_idx_ref},dt,avg_time_diff), 'Interpreter', 'none')


set(gcf, 'paperposition', [0 0 8 5])
set(gcf, 'papersize', [8 5])
saveas(gcf, sprintf('%scell_%dHz_cross-correlogram_%s_%s.pdf', CELL_TYPE,fps,channelNames{ch_idx}, channelNames{ch_idx_ref}))
saveas(gcf, sprintf('%scell_%dHz_cross-correlogram_%s_%s.png', CELL_TYPE,fps,channelNames{ch_idx}, channelNames{ch_idx_ref}))


