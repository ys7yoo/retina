
% compare STA slice and fitted ellipses
close all

%% plot STAs
gridT = (-sta_num_samples+1:0)/sampling_rate;

subplot(221)
channel_idx = 1;
plot(gridT,reshape(sta_all_channels{channel_idx},[sta_num_samples, width*height])); ylim([0 1]); box off
xlabel('pre spike (s)')
ylabel('STA')
title(sprintf('%s', channel_names{channel_idx}), 'Interpreter', 'none')

    
subplot(222)
channel_idx = 2;
plot(gridT,reshape(sta_all_channels{channel_idx},[sta_num_samples, width*height])); ylim([0 1]); box off
xlabel('pre spike (s)')
ylabel('STA')
title(sprintf('%s', channel_names{channel_idx}), 'Interpreter', 'none')




%for channel_idx=1:length(channel_names)
for channel_idx=1:2
    if isfield(RFs{channel_idx}, 'slice')
        slice_with_max_response = RFs{channel_idx}.slice;

        max_response_time = gridT(slice_with_max_response);

        STA_slice = get_slice_from_STA(sta_all_channels{channel_idx}, sta_num_samples, width, height, slice_with_max_response)

        subplot(2,2,channel_idx+2)
        imshow(STA_slice)
        title(sprintf('max response at t=%.2f', max_response_time), 'Interpreter', 'none')

        axis xy

        hold on
        plot_RF(RFs{channel_idx})
    end
    
%     subplot(133) % overlay RF ellipses
%     hold on
%     plot_RF(RFs{channel_idx})
%     axis equal
end


% save figure


set(gcf, 'paperposition', [0 0 24 20])
set(gcf, 'papersize', [24 20])

fig_filename= sprintf('RFs_%s_vs_%s.pdf',channel_names{1}(4:end),channel_names{2}(4:end));
saveas(gcf, fig_filename)
    
