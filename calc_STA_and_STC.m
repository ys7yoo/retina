function [sta, stc_eig_val, stc_eig_vec, S] = calc_STA_and_STC(stim, spike_train, n, sta_to_project_out)

% input:
%       Stim = (time) x (space)
%       spike_train = (time) x (spikes)
%       n = number of samples for analysis
%       project_out_sta = bool whether to project out sta or not

if nargin<4
    sta_to_project_out = [];
end

%% store spike-tiggered stims in to X with spike numbers in spikes
[X, spikes, num_total_spikes] = collect_spike_triggered_stim(stim, spike_train, n);

%% calc STA
sta = spikes'*X/num_total_spikes;

%% calc STC (FINAL ALGORITHM)
if nargout>1
    %% preprocessing  
%     % MEAN SHOULD BE SUBTRACTED BEFORE CALLING!
%     % 1) subtract mean
%     X = bsxfun(@minus, X, sta);
        
    % 2) project out sta, if requested
    if ~isempty(sta_to_project_out)
        X = project_out_components(X, sta_to_project_out);
    else
        X = project_out_components(X, sta);       
    end
       
    %% calc STC
    switch nargout 
        
        case 2   % STA and STC eigen value only (COVARIANCE NOT NEEDED)
            % weight X with number of spikes
            Xs = bsxfun(@times, X, sqrt(spikes));
            
            [~, D, ~] = svd(Xs);
            stc_eig_val = diag(D).^2;

        case 3  % STA and STC eigen value & eiven vectors (COVARIANCE NOT NEEDED)
            % weight X with number of spikes
            Xs = bsxfun(@times, X, sqrt(spikes));

            [~, D, stc_eig_vec] = svd(Xs);
            
            stc_eig_val = diag(D).^2;
            
            stc_eig_vec = flip_column_sign(stc_eig_vec, sta);  % flip according to sta (for better visualization)
            
        case 4 % full algorithm with covariance
            S = X'*bsxfun(@times, X, spikes);  % DO NOT DIVIDE by num_total_spikes FOR BETTER NUMERICAL PRECISION. So, actually S is a scatter matrix.

            [stc_eig_vec, d, ~] = svd(S);

            stc_eig_val = diag(d);            
            
            stc_eig_vec = flip_column_sign(stc_eig_vec, sta);  % flip according to sta (for better visualization)
    end

end
        
% % 
% % switch nargout
% %     case 1  % STA only
% %         
% %     case 2  % STA and STC eigen value only (COVARIANCE NOT NEEDED)
% % 
% %         % subtract mean & scale down 
% %         %SS = bsxfun(@minus, SS, sta) / (num_bins-1);
% %         X = bsxfun(@minus, X, sta);       % DO NOT DIVIDE FOR BETTER NUMERICAL PRECISION
% % 
% %         [~, D, ~] = svd(X);
% %         stc_eig_val = diag(D).^2;
% % 
% %     case 3  % STA and STC eigen value & eiven vectors (COVARIANCE NOT NEEDED)
% % 
% %         % subtract mean & scale down 
% %         X = bsxfun(@minus, X, sta);       % DO NOT DIVIDE FOR BETTER NUMERICAL PRECISION
% % 
% %         [~, D, stc_eig_vec] = svd(X);
% %         stc_eig_val = diag(D).^2;
% % 
% %         
% %         
% %     case 4  % STA and full STC
% %         %% calc STC
% %         
% %         % subtract mean
% % %         X = bsxfun(@minus, X, sta);       % DO NOT DIVIDE FOR BETTER NUMERICAL PRECISION
% % %         
% % %         covariance = X'*X;
% %          
% %         % implementation 1 (when spikes are either 1 or 0, Matlab cov function is faster!
% %         %covariance = cov(X);
% % 
% %         % implementation 2
% %         covariance = X'*bsxfun(@times, X, spikes) - sta*sta'*num_total_spikes;  % DO NOT DIVIDE FOR BETTER NUMERICAL PRECISION
% %         %covariance = X'*bsxfun(@times, X, spikes)/(num_total_spikes-1) - sta*sta'*num_total_spikes/(num_total_spikes-1);
% % 
% %         % implementation 3
% %         %rowlen = size(SS,2);
% %         %stc = X'*(X.*repmat(spikes,1,rowlen))/(num_total_spikes-1) - sta*sta'*num_total_spikes/(num_total_spikes-1);
% %         
% %         
% %         [stc_eig_vec, d, ~] = svd(covariance);
% %     
% %         stc_eig_val = diag(d);
% %         %stc_eig_val = stc_eig_val(stc_eig_val>1e-5);
% %     
% % end

% dim_space = size(stim,2);     % dimension of the stimulus
% if dim_space ~= 1
%     sta = reshape(sta, n, dim_space);
%     
% % %     if nargout > 1
% % %         if dim ~= 1
% % %             % unpack STC into nxn blocks
% % %             for i=1:dim
% % %                 for j=i:dim
% % %                     STC{i,j} = stc(n*(i-1)+1:n*i, n*(j-1)+1:n*j);
% % %                     STC{j,i} = STC{i,j};
% % %                 end
% % %             end
% % %             stc= STC;
% % %         end
% % %     end
% end


return 






%% to debug

[sta1, stc1] = calc_STA_and_STC(Stim,sp,n);  

clf
subplot(121)
plot(sta1)
subplot(122)
imshow(stc1)
%%

[sta2, stc2] = calc_STA_and_STC([Stim Stim],sp,n);  

clf
subplot(121)
plot(sta2)
subplot(122)
imshow(stc2)

