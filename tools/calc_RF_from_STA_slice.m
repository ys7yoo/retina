function [pos_RF, neg_RF, largest_RF] = calc_RF_from_STA_slice(STA, X, Y)

[T, num_pixels] = size(STA);

[YY, XX] = meshgrid(1:X,1:Y);
XX=XX(:);
YY=YY(:);



%% Step 1. calc noise sigma
sig = std(STA(:));







%% For each slice
r=floor(sqrt(T));
c=ceil(T/r);
cnt=1;

max_val = max(STA(:));
min_val = min(STA(:));

cnt_pos_RF = 0; clear pos_RF
cnt_neg_RF = 0; clear neg_RF
for t=1:T
    slice  = STA(t,:);
    subplot(r,c,cnt)
    imagesc(reshape(slice, X,Y), [min_val max_val])
    cnt = cnt + 1;
    colorbar
    
    %% Step 2. calc weighted centers for out of 99% significant interval (out of
    % noise)
    [pos_center, pos_cov, neg_center, neg_cov] =  calc_weighted_centers(slice, X, Y, 2.58*sig);
    hold on
    
    %% Step 3. plot ellipses (95% significant)
    if sum(isnan(pos_center))==0
        plot(pos_center(1), pos_center(2), '+r', 'markersize', 5)
        eig_values=plot_ellipse(pos_center, pos_cov, 'r-');
        
        cnt_pos_RF = cnt_pos_RF + 1;
        
        pos_RF{cnt_pos_RF}.mean= pos_center;
        pos_RF{cnt_pos_RF}.cov= pos_cov;
        pos_RF{cnt_pos_RF}.eig= sum(eig_values); % size of RF
    end
    if sum(isnan(neg_center))==0
        plot(neg_center(1), neg_center(2), '+b', 'markersize', 5)
        eig_values=plot_ellipse(neg_center, neg_cov, 'b-');
     
        cnt_neg_RF = cnt_neg_RF + 1;
        
        neg_RF{cnt_neg_RF}.mean= neg_center;
        neg_RF{cnt_neg_RF}.cov= neg_cov;
        neg_RF{cnt_neg_RF}.eig= sum(eig_values); % size of RF
    end
    
    title (sprintf('time bin %d',t))
    
    
    %% Store largest mean and cov of largest RF
    
end

colormap gray



%% calc and return largest RF of this cell
if nargout>2
    
    num_pos_RF = length(pos_RF);
    num_neg_RF = length(neg_RF);
    if num_pos_RF==0 && num_neg_RF==0
        largest_RF.type = 'NONE';
        largest_RF.mean= NaN;
        largest_RF.cov= NaN;
        largest_RF.eig= NaN; %
    else
        % find the largest pos eig
        pos_eigs=[];
        if num_pos_RF > 0 
            for i=1:num_pos_RF       
                pos_eigs=[pos_eigs pos_RF{i}.eig];
            end
            [max_pos_eig, max_idx] = max(pos_eigs);
        else
            max_pos_eig = -1;
        end

        if max_pos_eig>0
            largest_pos_RF.mean = pos_RF{max_idx}.mean;
            largest_pos_RF.cov = pos_RF{max_idx}.cov;
            largest_pos_RF.eig = pos_RF{max_idx}.eig;
        else
            largest_pos_RF.mean = NaN;
            largest_pos_RF.cov = NaN;
            largest_pos_RF.eig = 0;
        end

        % find the largest neg eig    
        neg_eigs=[];
        if num_neg_RF > 0
            for i=1:num_neg_RF
                neg_eigs=[neg_eigs neg_RF{i}.eig];
            end
            [max_neg_eig, max_idx] = max(neg_eigs);
        else
            max_neg_eig = -1;
        end

        if max_neg_eig>0
            largest_neg_RF.mean = neg_RF{max_idx}.mean;
            largest_neg_RF.cov = neg_RF{max_idx}.cov;
            largest_neg_RF.eig = neg_RF{max_idx}.eig;
        else
            largest_neg_RF.mean = NaN;
            largest_neg_RF.cov = NaN;
            largest_neg_RF.eig = 0;
        end

        % return the larger of the two
        if max_pos_eig > max_neg_eig
            largest_RF = largest_pos_RF
            largest_RF.type = 'ON'
        else
            largest_RF = largest_neg_RF
            largest_RF.type = 'OFF'
        end
    end
    
end



return 
%% 

set(gcf, 'paperposition', [0 0 24 18])
set(gcf, 'papersize', [24 18])

saveas(gcf, sprintf('STA_slice.png'))
saveas(gcf, sprintf('STA_slice.pdf'))
