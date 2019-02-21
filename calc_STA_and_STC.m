function [sta, stc_eig_val, stc_eig_vec, S] = calc_STA_and_STC(stim, spike_train, n, project_out_sta)

% input:
%       Stim = (time) x (space)
%       spike_train = (time) x (spikes)
%       n = number of samples for analysis
%       project_out_sta = bool whether to project out sta or not

if nargin<4
    project_out_sta = true;
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
    if project_out_sta
        X = project_out_components(X, sta);
    end
       
    %% calc STC
    switch nargout 
        
        case 2   % STA and STC eigen value only (COVARIANCE NOT NEEDED)
            stc_eig_val = calc_STC(X, spikes);

        case 3  % STA and STC eigen value & eiven vectors (COVARIANCE NOT NEEDED)
            [stc_eig_val, stc_eig_vec] = calc_STC(X, spikes);
            
            stc_eig_vec = flip_column_sign(stc_eig_vec, sta);  % flip according to sta (for better visualization)
            
        case 4 % full algorithm with covariance
            [stc_eig_val, stc_eig_vec, S] = calc_STC(X, spikes)
                        
            stc_eig_vec = flip_column_sign(stc_eig_vec, sta);  % flip according to sta (for better visualization)
            
    end

end
        

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

