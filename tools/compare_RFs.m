
% compare STA slice and fitted ellipses
clf

%for channel_idx=1:length(channel_names)
for channel_idx=1:2
    slice_with_max_response = RFs{channel_idx}.slice;
    STA_slice = get_slice_from_STA(sta_all_channels{channel_idx}, sta_num_samples, width, height, slice_with_max_response)

    subplot(1,3,channel_idx)
    imshow(STA_slice)
    title(sprintf('%s (slice %d)', channel_names{channel_idx}, slice_with_max_response), 'Interpreter', 'none')
    
    axis xy
    
    hold on
    plot_RF(RFs{channel_idx})
    
    subplot(133) % overlay RF ellipses
    hold on
    plot_RF(RFs{channel_idx})
    axis equal
end
