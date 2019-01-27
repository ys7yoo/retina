function [pos_RFs, neg_RFs, strongest_RF] = calc_RF_from_STA_slice(STA, X, Y, fps)

[T, num_pixels] = size(STA);
gridT = (-T+1:0)/fps;

[YY, XX] = meshgrid(1:X,1:Y);
XX=XX(:);
YY=YY(:);

height = max(YY);
width  = max(XX);

%% Step 1. calc noise sigma
sig = std(STA(:));


%% find positive / negative peak slides
[max_val, max_slice_idx] =  max(max(STA, [], 2));
[min_val, min_slice_idx] =  min(min(STA, [], 2));

% decide cell type according to the strength of the signal
if (max_val-0.5) > (0.5 - min_val)
    cell_type = 1;  %'ON'
else
    cell_type = 0;  %'OFF'
end
    



%% Smooth STA for each slice

%STA = smooth_STA_slice(STA, 0.5);



%% For each slice
%c=ceil(sqrt(T));
%c=floor(sqrt(T));
c=5;
r=ceil(T/c);
cnt=1;

max_val = max(STA(:));
min_val = min(STA(:));

pos_RF = [];
neg_RF = [];
cnt_pos_RF = 0; pos_RFs = []; %clear pos_RFs
cnt_neg_RF = 0; neg_RFs = []; %clear neg_RFs
for t=1:T
    slice  = STA(t,:);
    subplot(r,c,cnt)
    imagesc(reshape(slice, X,Y), [min_val max_val])
    %axis tight equal
    
    %colorbar
    
    %% Step 2. calc weighted centers 
    % noise)
    slice_smoothed = smooth_STA_slice(slice, 1.0);
    [pos_center, pos_cov, num_pos_pixels, neg_center, neg_cov, num_neg_pixels] =  calc_weighted_centers(slice_smoothed, X, Y, 2.58*sig);
    hold on
    
    %% Step 3. plot ellipses (95% significant)
    if sum(isnan(pos_center))==0
        %plot(pos_center(1), pos_center(2), '+r', 'markersize', 5)
        eig_values=plot_ellipse(pos_center, pos_cov, 'r-');
        
        if sum(eig_values) > 0
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
        %plot(neg_center(1), neg_center(2), '+b', 'markersize', 5)
        eig_values=plot_ellipse(neg_center, neg_cov, 'b-');
     
        if sum(eig_values) > 0 
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
    
end

colormap gray







return 
%% 

set(gcf, 'paperposition', [0 0 24 18])
set(gcf, 'papersize', [24 18])

saveas(gcf, sprintf('STA_slice.png'))
saveas(gcf, sprintf('STA_slice.pdf'))
