
%% Let's analyze all the cells


%% load data
loadData

%% run STA for all channels!


mkdir(sprintf('Num%d',NUM_EXP))
for i=1:length(channelNames)
    channelName = channelNames{i}
    spikeTime=eval(channelName);

    % choose spike during the stim only!
    W = 0.5 % Window in sec
    [STA, gridT] = calcSTAintp(stim, A1a, spikeTime, W, fps);



    %% plot STA
    clf
    subplot(121)
    imshow(STA')
    xlabel('t')
    ylabel('pixel')
    title('STA')

    % check variance to see if there is any change 
    STAvar = var(STA,[],1);
    subplot(222)
    imshow(reshape(STAvar,height,width),[])
    xlabel('x')
    ylabel('y')
    title('variance across time')

    % mark the pixel with largest variance 
    [mm, maxIdx] = max(STAvar);
    [XX,YY]= meshgrid(1:13);

    hold on;
    plot(XX(maxIdx),YY(maxIdx),'+r')
    %axis xy

    % plot STA for the pixel with the largest variance

    subplot(224)
    plot(gridT, STA(:,maxIdx)-mean(STA(:,maxIdx)))
    hold on
    plot(gridT([1,end]),[0 0], ':k') % plot mean as reference
    xlabel('t')
    ylabel('STA')
    title('STA for the pixel with the largest variance')
    %xaxis=(-W+1:0)/fps;
    %set(gca,'xticklabel',(W-1*get(gca,'xtick'))/fps)
    box off

    set(gcf, 'paperposition', [0 0 9 8])
    set(gcf, 'papersize', [9 8])
    saveas(gcf, sprintf('Num%d/STA_%s.pdf',NUM_EXP,channelName))
end
