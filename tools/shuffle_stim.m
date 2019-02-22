function  [X, spikes, num_total_spikes] = shuffle_stim(stim, spikeTrain, num_samples_per_window, random_shift_range)

%% set default params for shuffle
if nargin<4
    random_shift_range = num_samples_per_window * [1 10];
end

%%
shift_min = random_shift_range(1);
shift_max = random_shift_range(2); 

%% random shift
random_shift = round(diff(random_shift_range) * rand(1)  + random_shift_range(1));

%% calc ev of STC
[X, spikes, num_total_spikes] = collect_spike_triggered_stim(stim(random_shift+1:end-(shift_max-random_shift),:), spikeTrain(1:end-shift_max,:), num_samples_per_window);
