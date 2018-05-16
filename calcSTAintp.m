function [STA, gridT] = calcSTAintp(stim, stimTime, spikeTime, W, fps) 

% stimTime - stimlus time 
% spikeTime - spike Time 
% W -     % window size in sec


gridT = -W:1/fps:0;




% % % convert to index
% % spikeTimeIdx = round(spikeTime*fps);
% % 
% % [height, width, T] = size(stim);
% % 
% % if (max(spikeTimeIdx) > T)
% %     warning('Maximum spike time exceeds recoding time. Check!!!')
% %     
% %     % [QUICK FIX]
% %     % remove spikeTimeIdx out of stim range
% %     spikeTimeIdx = spikeTimeIdx(spikeTimeIdx<=T);
% % end
%%  calc STA

% Step 1. find stimulus pattern that triggers a spike
stimTriggeredSpike = [];

% remove spikeTimeIdx less than W
% spikeTimeIdx = spikeTimeIdx(spikeTimeIdx>W);

%numSpike = length(spikeTimeIdx);

%stimTriggeredSpike = zeros (height,width,W,numSpike);
cnt = 1;
for st=spikeTime
    
    stimInterp = interp1(stimTime, stim, st-gridT);  % 1st dim should be time!

    if sum(sum(isnan(stimInterp)))==0
        stimTriggeredSpike(:,:,cnt) = stimInterp;
        cnt = cnt + 1;
    end
end

STA = mean(stimTriggeredSpike,3); % 2-dim array (height x width) x W



                                 