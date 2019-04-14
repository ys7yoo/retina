function RFs = calc_and_plot_RF_from_STA(sta_all_channels, sta_num_samples, width, height, fps, channel_names)
%% calc and plot RF from STA (2019. 1. 27) => separated to a function (2019. 2.13)
%channel_index_to_analyze = 1:length(channel_names);

FLIP_XY=true;

addpath tools
clear RFs;
%for n=1:length(sta_all_channels)
cnt = 1;

for n = 1:length(sta_all_channels)
    %n = channel_index_to_analyze
    %% 
    disp(sprintf('processing %s...',channel_names{n}))
    
    
    %% pre-processing 
    clf
    hist(sta_all_channels{n}(:))
    sig = std(sta_all_channels{n}(:));
    ylim = get(gca,'ylim');
    hold on;plot(0.5+2.58*sig*[1 1], ylim, 'r--');plot(0.5-2.58*sig*[1 1], ylim, 'r--')
    box off
    %title(sprintf('histogram of STA (\\sigma=%.2f)', sig))
    %sta_all_channels{n}
    
    %% plot outside of m+- sig
    
    
    %%
    clf
    
    
    [pos_RF, neg_RF, strongest_RF] = calc_RF_from_STA_slice(sta_all_channels{n}, sta_num_samples, width, height, fps, FLIP_XY);
    strongest_RF.channel_name = channel_names{n};
    set(gcf, 'paperposition', [0 0 24 9])
    set(gcf, 'papersize', [24 9])

    saveas(gcf, sprintf('RF_from_STA_%s.png',channel_names{n}))
    saveas(gcf, sprintf('RF_from_STA_%s.pdf',channel_names{n}))
    
    % saves the largest RF from each cell
    RFs{n} = strongest_RF;
    %RFs{cnt} = strongest_RF;
    %cnt = cnt + 1;
    
end
