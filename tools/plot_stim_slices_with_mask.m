function plot_stim_slices_with_mask(stim, sta_num_samples, XX, YY, w, h, FLIP_XY, MEA_xy, RF)

if nargin< 7 
    FLIP_XY = false;
end

if nargin<8
    MEA_xy = [];
end

if nargin<9
    RF=[];
end

num_pixels = length(XX);

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
    axis xy
    hold on
    
    %% mark MEA
    if ~isempty(MEA_xy)
        
        if FLIP_XY
            plot(MEA_xy(1), MEA_xy(2), 'go')
        else
            plot(MEA_xy(2), MEA_xy(1), 'go')
        end
    end
    
    %% plot RF
    if ~isempty(RF) && isfield(RF, 'slice')
        if RF.slice == i
%             plot_RF(RF, FLIP_XY)
            switch RF.type
                case 'ON'            
                    plot_ellipse(RF.mean, RF.cov, 'r--')
                case 'OFF'
                    plot_ellipse(RF.mean, RF.cov, 'b--')
            end
        end
    end
    
    %% zoon in masks
    %axis([min(XX(:)) max(XX(:)) min(YY(:)) max(YY(:))])
    
    %%
    if FLIP_XY
        view([90 -90])
    end
    
    box off
end

return 


%%
