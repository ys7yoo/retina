function plot_histogram_for_each_slice(STA, T, fps)




%% Step 1. calc noise sigma
sig = std(STA(:));

%% reshape to get slice

STA = reshape(STA, T, []);

%[T, num_pixels] = size(STA);
gridT = (-T+1:0)/fps;

%% find positive / negative peak slides
[max_val, max_slice_idx] =  max(max(STA, [], 2));
[min_val, min_slice_idx] =  min(min(STA, [], 2));


bins = linspace(min_val, max_val, 10)
%% Smooth STA for each slice

%STA = smooth_STA_slice(STA, 0.5);



%% For each slice
%c=ceil(sqrt(T));
%c=floor(sqrt(T));
c=5;
r=ceil(T/c);
cnt=1;

% max_val = max(STA(:));
% min_val = min(STA(:));

pos_RF = [];
neg_RF = [];
cnt_pos_RF = 0; pos_RFs = []; %clear pos_RFs
cnt_neg_RF = 0; neg_RFs = []; %clear neg_RFs
strongest_pos_RF = [];
strongest_neg_RF = [];
for t=1:T
    slice  = STA(t,:);
    subplot(r,c,t)
    hist(slice, bins)
    
    title (sprintf('%.0f ms', gridT(t)*1000))
    cnt = cnt + 1;
    
end

return

    if ~FLIP_XY
        imagesc(reshape(slice, height, width), [min_val max_val])
    else
        imagesc(reshape(slice, height, width)', [min_val max_val])
    end
    %axis tight equal