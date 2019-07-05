
%% Let's analyze all the cells


%% load data
loadData

%% run STA for all channels!


mkdir(sprintf('Num%d',NUM_EXP))
for i=1:length(channelNames)
    channelName = channelNames{i}
    spikeTime=eval(channelName);

    % choose spike during the stim only!
    W = 0.5; % Window in sec
    [STA, gridT] = calcSTAintp(stim, A1a, spikeTime, W, fps);



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

    set(gcf, 'paperposition', [0 0 9 8])
    set(gcf, 'papersize', [9 8])
    saveas(gcf, sprintf('Num%d/STA_%s.pdf',NUM_EXP,channelName))
end
