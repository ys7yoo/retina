
%% Let's analyze all the cells


%% load data
loadData

%% run STA for all channels!


mkdir(sprintf('Num%d',NUM_EXP))
for i=1:length(channelNames)
    channelName = channelNames{i}
    spikeTime=eval(channelName);

    %% call STA here
    W = 0.5*fps;    % window size 
    %W = 0.2*fps;    % window size 
    STA = calcSTA(stim, spikeTime-A1a(1), fps, W);



    %% plot STA
    STAstack = reshape(STA,[height*width, W]);
    clf
    subplot(121)
    imshow(STAstack)
    xlabel('t')
    ylabel('pixel')
    title('STA')

    % check variance to see if there is any change 
    STAvar = var(STAstack,[],2);
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
    plot(STAstack(maxIdx,:)-mean(STAstack(maxIdx,:)))
    xlabel('t')
    ylabel('STA')
    title('STA for the pixel with the largest variance')
    xaxis=(-W+1:0)/fps;
    set(gca,'xticklabel',(W-1*get(gca,'xtick'))/fps)


    set(gcf, 'paperposition', [0 0 9 8])
    set(gcf, 'papersize', [9 8])
    saveas(gcf, sprintf('Num%d/STA_%s.pdf',NUM_EXP,channelName))
end
