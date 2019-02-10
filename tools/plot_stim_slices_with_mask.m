function plot_stim_slices_with_mask(stim, sta_num_samples, XX, YY, w, h)

num_pixels = length(XX)

%%
stim_max = max(stim(:));
stim_min = min(stim(:));

st = reshape(stim, sta_num_samples, []);

c = 5;
r = ceil(sta_num_samples/c);





for i=1:sta_num_samples
    subplot(r,c,i)
    
    slice = nan*zeros(h,w);
    % fill slice
    for j=1:num_pixels
        slice(YY(j),XX(j)) = st(i,j);
    end
    
    
    imagesc(slice, [stim_min stim_max])
    %surf(reshape(st(i,:,:),h,w))
    colormap gray
    axis ij
end

return 


%%
