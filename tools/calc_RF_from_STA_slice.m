function [pos_RFs, neg_RFs, strongest_RF] = calc_RF_from_STA_slice(STA, T, width, height, fps, FLIP_XY)



[XX, YY] = meshgrid(1:width,1:height);
% XX=XX(:);
% YY=YY(:);
% 
% height = max(YY);
% width  = max(XX);

%% Step 1. calc noise sigma
sig = std(STA(:));

%% reshape to get slice

STA = reshape(STA, T, []);

%[T, num_pixels] = size(STA);
gridT = (-T+1:0)/fps;

%% find positive / negative peak slides
[max_val, max_slice_idx] =  max(max(STA, [], 2));
[min_val, min_slice_idx] =  min(min(STA, [], 2));

% decide cell type according to the strength of the signal
if (max_val-0.5) > (0.5 - min_val)
    cell_type = 1;  %'ON'
    % find peak center
    sta_slice = STA(max_slice_idx,:);
    [~, maxIdx] = max(sta_slice);
    peak_center = [max_slice_idx, floor((maxIdx-1)/height)+1, mod(maxIdx, height)];
else
    cell_type = 0;  %'OFF'
    % find peak center
    sta_slice = STA(min_slice_idx,:);
    [~, minIdx] = min(sta_slice);
    peak_center = [min_slice_idx, floor((minIdx-1)/height)+1, mod(minIdx, height)];
end
    





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
strongest_RF = [];
for t=1:T
    slice  = STA(t,:);
    subplot(r,c,cnt)
    if ~FLIP_XY
        imagesc(reshape(slice, height, width), [min_val max_val])
    else
        imagesc(reshape(slice, height, width)', [min_val max_val])
    end
    %axis tight equal
    axis xy
    %colorbar
    box off
    axis off
    
    %% Step 2. calc weighted centers 
    % noise)
    slice_smoothed = smooth_STA_slice(slice, 1.0, width, height);
    [pos_center, pos_cov, num_pos_pixels, neg_center, neg_cov, num_neg_pixels] =  calc_weighted_centers(slice_smoothed, width, height, 0.5, 2.58*sig);
    %[pos_center, pos_cov, num_pos_pixels, neg_center, neg_cov, num_neg_pixels] =  calc_weighted_centers(slice_smoothed, width, height, 0.5, 1.96*sig);  % 95% confidence interval
    hold on
    
    %% Step 3. plot ellipses (95% significant)
    if sum(isnan(pos_center))==0
        eig_values = eig(pos_cov);
        
        if sum(eig_values>0) == 2
            %plot(pos_center(1), pos_center(2), '+r', 'markersize', 5)
            plot_ellipse(pos_center, pos_cov, 'r-', FLIP_XY);
        
            
            cnt_pos_RF = cnt_pos_RF + 1;
        
            pos_RFs{cnt_pos_RF}.mean= pos_center;
            pos_RFs{cnt_pos_RF}.cov= pos_cov;
            pos_RFs{cnt_pos_RF}.eig= sum(eig_values); % size of RF along the principal axes
            pos_RFs{cnt_pos_RF}.slice = t;
            pos_RFs{cnt_pos_RF}.num_pixels = num_pos_pixels;
            
            if (cell_type == 1) && (t == max_slice_idx)
                strongest_RF = pos_RFs{cnt_pos_RF};
                strongest_RF.type = 'ON';
            end
        end
        
    end
    if sum(isnan(neg_center))==0
        eig_values = eig(neg_cov);
     
        if sum(eig_values>0) == 2
            %plot(neg_center(1), neg_center(2), '+b', 'markersize', 5)
            plot_ellipse(neg_center, neg_cov, 'b-', FLIP_XY);
        
            cnt_neg_RF = cnt_neg_RF + 1;

            neg_RFs{cnt_neg_RF}.mean= neg_center;
            neg_RFs{cnt_neg_RF}.cov= neg_cov;
            neg_RFs{cnt_neg_RF}.eig= sum(eig_values); % size of RF
            neg_RFs{cnt_neg_RF}.slice = t;
            neg_RFs{cnt_neg_RF}.num_pixels = num_neg_pixels;
            
            if (cell_type == 0) && (t == min_slice_idx)
                strongest_RF = neg_RFs{cnt_neg_RF};
                strongest_RF.type = 'OFF';
            end
            
        end
        
    end
    
    %title (sprintf('time bin %d',t))
    title (sprintf('%.0f ms', gridT(t)*1000))
    cnt = cnt + 1;
    set(gca, 'xlim', [1 height]);
    set(gca, 'ylim', [1 height]);
    %axis equal
    
    %% Store largest mean and cov of largest RF
    
    
    %% if fails
    if ~isfield(strongest_RF, 'type')
        if cell_type == 1
            strongest_RF.type = 'ON';
        else
            strongest_RF.type = 'OFF';
        end
    end
    
end

colormap gray



strongest_RF.peak_center = peak_center;




return 
%% 

set(gcf, 'paperposition', [0 0 24 18])
set(gcf, 'papersize', [24 18])

saveas(gcf, sprintf('STA_slice.png'))
saveas(gcf, sprintf('STA_slice.pdf'))
