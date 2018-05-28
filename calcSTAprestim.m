function [STA, gridT] = calcSTAprestim(stim, stimTime, spikeTime, W, fps) 

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
for st=1:length(spikeTime)
  if spikeTime(st) > stimTime(W*fps+1) && spikeTime(st) < stimTime(end)
    PreStims = [];
    PreStimInx = [];
    Stims = spikeTime(st)-stimTime;
    for i=1:length(Stims)
      if Stims(i)>0
        PreStims(1,i) = Stims(i);
      end
    end
        minPreStim = min(PreStims);
        PreStimInx = find(Stims==minPreStim);
       for k=1:length(gridT)
        stimTriggeredSpike(k,:,cnt) = stim(PreStimInx-length(gridT)+k,:); 
       end
        cnt = cnt + 1;   
  end
end

STA = mean(stimTriggeredSpike,3); % 2-dim array (height x width) x W



                                 