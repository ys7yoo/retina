
%% Let's analyze one cell
% ON cell (ch47b,58a,58b,66c,66d,76a,77b,85c,86c)
% OFF cell (ch16b,17a,17b,28a,37a,46a,47a,48b,55b,65b,66a,66b,67a,68b,77a,85a,86a)


%% load data
loadData

%% choose one channel from channelNames (MANUAL for now)
disp(channelNames)
channelName = 'ch_47b'
%channelName = 'ch_48a'
%channelName = 'ch_65b'
%channelName = 'ch_66a'
%channelName = 'ch_66b'
spikeTime=eval(channelName)

%% call here

% choose spike during the stim only!
W = 0.5; % Window in sec
[STA, gridT] = calcSTAintp(stim, A1a, spikeTime, W, fps);


% % W = 0.5*fps;    % window size 
% % %W = 0.2*fps;    % window size 
% % % choose spike during the stim only!
% % validIdx = spikeTime>A1a(1) & spikeTime<=A1a(end);
% % spikeTime=spikeTime(validIdx);
% % STA = calcSTA(stim, spikeTime-A1a(1), fps, W);



%% plot STA
clf
subplot(121)
imshow(STA')
xlabel('t')
ylabel('pixel')
title('STA')

STAmean = mean(STA,1);
subplot(322)
imshow(reshape(STAmean,height,width)',[])   % stim is rotate 90 degree!
xlabel('y') %xlabel('x')         % swap x and y
ylabel('x') %ylabel('y')
title(sprintf('mean STA across time'))
axis on xy


% mark the pixel with largest variance 
[XX,YY]= meshgrid(1:height);  % grid for ploting
[mm, maxIdx] = max(STAmean);

hold on;
%plot(XX(maxIdx),YY(maxIdx),'+r', 'markersize', 8, 'linewidth', 2)
plot(YY(maxIdx),XX(maxIdx),'+r', 'markersize', 8, 'linewidth', 2)

[mm, minIdx] = min(STAmean);
[XX,YY]= meshgrid(1:height);

hold on;
%plot(XX(minIdx),YY(minIdx),'xb', 'markersize', 8, 'linewidth',4)
plot(YY(minIdx),XX(minIdx),'xb', 'markersize', 8, 'linewidth',4)



subplot(324)
plot(gridT, STA(:,maxIdx),'r')
xlabel('t')
ylabel('STA')
title('STA for the pixel with the largest variance')
box off

subplot(326)
plot(gridT, STA(:,minIdx),'b')
xlabel('t')
ylabel('STA')
title('STA for the pixel with the largest variance')
box off
    
    saveas(gcf, sprintf('Num%d_STA_%s.pdf',NUM_EXP,channelName))
    
    
% % %% plot STA
% % clf
% % subplot(121)
% % imshow(STA')
% % xlabel('t')
% % ylabel('pixel')
% % title('STA')
% % 
% % STAmean = mean(STA,1);
% % subplot(322)
% % imshow(reshape(STAmean,height,width),[])
% % xlabel('x')
% % ylabel('y')
% % title(sprintf('mean across time %.2f',mean(STAmean)))
% % 
% % 
% % % check variance to see if there is any change 
% % STAvar = var(STA,[],1);
% % subplot(324)
% % imshow(reshape(STAvar,height,width),[])
% % xlabel('x')
% % ylabel('y')
% % title('variance across time')
% % 
% % % mark the pixel with largest variance 
% % [mm, maxIdx] = max(STAvar);
% % [XX,YY]= meshgrid(1:13);
% % 
% % hold on;
% % plot(XX(maxIdx),YY(maxIdx),'+r')
% % %axis xy
% % 
% % % plot STA for the pixel with the largest variance
% % 
% % subplot(326)
% % plot(gridT, STA(:,maxIdx)-mean(STA(:,maxIdx)))
% % hold on
% % plot(gridT([1,end]),[0 0], ':k') % plot mean as reference
% % xlabel('t')
% % ylabel('STA')
% % title('STA for the pixel with the largest variance')
% % %xaxis=(-W+1:0)/fps;
% % %set(gca,'xticklabel',(W-1*get(gca,'xtick'))/fps)
% % box off
% % 
% % set(gcf, 'paperposition', [0 0 9 8])
% % set(gcf, 'papersize', [9 8])
% % saveas(gcf, sprintf('Num%d/STA_%s.pdf',NUM_EXP,channelName))
