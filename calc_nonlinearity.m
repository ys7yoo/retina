function [generator, firing_rate] = calc_nonlinearity(stim, spikeTrain, sta, num_bin)

% bin stimulus
stim_bin = makeStimRows(stim, num_bin);
generator = stim_bin*sta(:);

% bin spikes
f_bin = makeStimRows(spikeTrain, num_bin);
% take avarage
firing_rate = mean(f_bin,2)/num_bin;

%plot(generator, firing_rate,'.k')

%xlabel('generator signal')
%ylabel('spikes / bin')

return 
%%

