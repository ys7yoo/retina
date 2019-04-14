function STA_slice = get_slice_from_STA(STA, sta_num_samples, width, height, slice_idx)

STA_reshaped=reshape(STA, [sta_num_samples, width, height]);

STA_slice = permute(STA_reshaped(slice_idx,:,:), [2 3 1]);


return

%% example 


