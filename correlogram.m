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


%% Example. autocorreogram for ch_22a
[bin, count, dt] = correlogram(ch_22a, ch_22a, [-0.2 0.2], 401, true);
bar(bin,count, 'k')
xlabel('\Deltat (sec)')
ylabel('count')
box off
title (sprintf('autocoorelogram (\\Deltat=%f)',dt))

%% Example. crosscorreogram for ch_22a and ch_31b
[bin, count,dt] = correlogram(ch_22a, ch_31b, [-0.2 0.2], 401);
bar(bin,count, 'k')
xlabel('\Deltat (sec)')
ylabel('count')
box off
title (sprintf('autocoorelogram (\\Deltat=%f)',dt))

