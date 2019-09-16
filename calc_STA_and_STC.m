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

% save here for further analysis
%save ch33b.mat X spikes num_total_spikes

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
            [stc_eig_val, stc_eig_vec, S] = calc_STC(X, spikes);
                        
            stc_eig_vec = flip_column_sign(stc_eig_vec, sta);  % flip according to sta (for better visualization)
            
    end

end
        

return 


%% for detail analysis for each channel
% data saved by the following command: save ch33b.mat X spikes num_total_spikes

clear
load ch33b.mat

% calc STA
sta = spikes'*X/num_total_spikes;

% calc STFC
[stc_eig_val, stc_eig_vec, S] = calc_STC(X, spikes);
                        
stc_eig_vec = flip_column_sign(stc_eig_vec, sta);  % flip according to sta (for better visualization)

% plot eig val
clf
subplot(221)
plot(stc_eig_val, 'o--')
box off

subplot(222)
hist(stc_eig_val)
box off


% calc scores for first r dimension
r = 2
score = X*stc_eig_vec(:,1:r)

subplot(223)
sc1 = scatter(score(:,1), score(:,2), '.'); %, '.', 'alpha', 0.2)
sc1.MarkerEdgeAlpha=0.4;

hold on
cov12 = cov(score(:,1), score(:,2))
plot_ellipse([0 0], cov12)

xlabel('ev1')
ylabel('ev2')
axis equal

axis ([-2.5 2.5 -2 2])

subplot(224)
sc1 = scatter(score(:,end), score(:,end-1), '.'); %, '.', 'alpha', 0.2)
sc1.MarkerEdgeAlpha=0.4;

hold on
plot_ellipse([0 0], cov(score(:,end), score(:,end-1)))

xlabel(sprintf('ev%d',length(stc_eig_val)))
ylabel(sprintf('ev%d',length(stc_eig_val)-1))
axis equal
axis ([-2.5 2.5 -2 2])


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

