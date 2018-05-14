function STA = calcSTA(stim, spikeTime, fps, W) 

% spikeTime - spike Time after stim onse
% W - 0.5*fps;    % window size 



% convert to index
spikeTimeIdx = round(spikeTime*fps);

[height, width, T] = size(stim);

if (max(spikeTimeIdx) > T)
    warning('Maximum spike time exceeds recoding time. Check!!!')
    
    % [QUICK FIX]
    % remove spikeTimeIdx out of stim range
    spikeTimeIdx = spikeTimeIdx(spikeTimeIdx<=T);
end
%%  calc STA

% Step 1. find stimulus pattern that triggers a spike
stimTriggeredSpike = [];

% remove spikeTimeIdx less than W
spikeTimeIdx = spikeTimeIdx(spikeTimeIdx>W);

numSpike = length(spikeTimeIdx);

stimTriggeredSpike = zeros (height,width,W,numSpike);
for i=1:length(spikeTimeIdx)
    idx= spikeTimeIdx(i);
        stimTriggeredSpike(:,:,:,i) = stim(:,:,idx-W+1:idx);
end

STA = mean(stimTriggeredSpike,4); % 3-dim array height x width x W



                                 