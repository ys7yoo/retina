%% Fit GLM with coupling filters!

%% set path to GLM tools

% This is where the GLMspiketools is installed
% You can get the package from https://github.com/ys7yoo/GLMspiketools
basedir = [getenv('HOME') '/src/GLMspiketools']

% Add a bunch sub-directories (with absoluate path names)
addpath([basedir '/glmtools_fitting/']);
addpath([basedir '/glmtools_misc/']);
addpath([basedir '/nlfuns/']);
addpath([basedir '/glmtools_spline/']);

addpath tools
%% load data! 

clear



% for loading 13x13 stim (20180828)
% load_data_13x13 

% for loading 26x26 stim (20180905)
load_data_26x26 


% % %% First, run STA independently for all the channels 
% % %W = 1; % Window in sec
% % %W = 0.8; % Window in sec
% % %W = 0.5; % Window in sec
% % W = 0.6; % Window in sec
% % 
% % 
% % clear STAall
% % clear STAs  % STAs will be stored here
% % 
% % for i = 1:length(channel_names)
% %     
% %     channel_name = channel_names{i}
% %     spike_time=eval(channel_name);
% %     
% %     [STA, gridT] = calcSTAprestim(stim, A1a, spike_time, W, fps);
% %     %save STA
% %     STAall{i} = STA;
% % 
% %     
% %     %figure(i)
% %     clf
% %     STA_max_var = find_STA_with_max_var(STA, gridT, width, height)
% %     
% %     STAs(:,i) = STA_max_var; % store STA with maximum variance (center)
% %     
% %     % save figure 
% %     set(gcf, 'paperposition', [0 0 12 10])
% %     set(gcf, 'papersize', [12 10])
% % 
% %     saveas(gcf, sprintf('%scell_%dHz_%s_STA.pdf', CELL_TYPE, fps, channel_name))
% %     
% % end
% % 
% % % plot STAs of multiple cells
% % clf
% % plot(gridT, STAs, 'linewidth',2)
% % ylabel('stim intensity')
% % xlabel('t')
% % title('STAs with the largest variance')
% % box off
% % 
% % legend(channel_names, 'location', 'NW','Interpreter', 'none'); legend boxoff
% % 
% % 
% % set(gcf, 'paperposition', [0 0 6 4])
% % set(gcf, 'papersize', [6 4])
% % 
% % saveas(gcf, sprintf('%scell_%dHz_STAs.pdf', CELL_TYPE, fps))



%% For further analysis, need to convert spike time to spike train

N = length(channel_names); % number of neurons 

binStim = size(stim,1);
assert (binStim==length(A1a))
spike_train = zeros(binStim,N);

for n = 1:N
    
    channel_name = channel_names{n};
    disp(channel_name)
    spike_time=eval(channel_name);
    
    
    % let's check for each stim time bin
    for i = 2:binStim
        t0 = A1a(i-1);
        t1 = A1a(i);

        % find stim time that occured during t0 and t1
        idx = find(spike_time>t0 & spike_time<=t1);

        spike_train(i,n) = length(idx);  % multiple spikes may occur in a bin
%         if ~isempty(idx)
%             %disp('found')
%             spike_train(i,n) = 1;
%         end

    end
end


clf
imagesc(spike_train', [0 1]); colormap gray

% % Simple code below does NOT work due to time jitter!!!
% % grab spike_time during stim & convert to bin idx
% idx = (spike_time >= A1a(1)) & (spike_time < A1a(end));
% binIdx = round(spike_time(idx)*fps);  
% %idx = ceil(spike_time*fps) % does it matter?
% spike_train(binIdx) = 1;
% spike_train = spike_train(:);
% 
% % cut spike train length to be the same as stim time length
% T = size(StimChosen,1);
% if length(spike_train) > T
%     spike_train = spike_train(1:T);
% end



%% Let's calc STA and STC
close all
%clear STA_max_var, max_idx

% sta params
%sta_num_samples = 16;
sta_num_samples = 10; % for 30Hz
gridT = (-sta_num_samples+1:0)/fps;


%channels_to_analyze = 1:length(channel_names);
%channels_to_analyze = 8;
%% 1. calc STA first (STC is too slow)

disp('running STA for all the channels...')
clear sta_all_channels
tic
for n=1:length(channel_names)
    
    disp(channel_names{n})
    
    % calc STA and STC
    [sta_all_channels{n}] = calc_STA_and_STC(stim, spike_train(:,n), sta_num_samples);
    
end
toc


%% lanch app to analyze individual channels and 
% decide which channels to analyze further


sta_app

%%
% ON cells: ch13c, ch37c, 58b, ch75d, ch77d
% OFF cells: ch16b, 26a, 27a, 37b, 47b

%channel_index_to_analyze = [8, 40, 63, 95] % ON 
%channel_index_to_analyze = [11, 22, 25, 39, 51] % OFF

% combined 
%channel_index_to_analyze = [8, 40, 63, 95, 11, 22, 25, 39, 51] % OFF

% better way to call
selected_channel_names = {'13c', '37c', '58b', '75d', '77d', '16b', '26a', '27a', '37b', '47b'}
addpath tools
channel_index_to_analyze = calc_channel_index(channel_names, selected_channel_names)

%% found 8 channels

channel_index_to_analyze = [25    28    91    92    97   101   103   106]

%% 
channel_index_to_analyze = 1:length(channel_names);
%% calc and plot RF from STA (2019. 1. 27)

FLIP_XY=true;

addpath tools
clear RFs;
%for n=1:length(sta_all_channels)
cnt = 1;
for n = channel_index_to_analyze
    %% 
    disp(sprintf('processing %s...',channel_names{n}))
    
    
    %% pre-processing 
    clf
    hist(sta_all_channels{n}(:))
    sig = std(sta_all_channels{n}(:));
    ylim = get(gca,'ylim');
    hold on;plot(0.5+2.58*sig*[1 1], ylim, 'r--');plot(0.5-2.58*sig*[1 1], ylim, 'r--')
    box off
    %title(sprintf('histogram of STA (\\sigma=%.2f)', sig))
    %sta_all_channels{n}
    
    %% plot outside of m+- sig
    
    
    %%
    clf
    
    
    [pos_RF, neg_RF, strongest_RF] = calc_RF_from_STA_slice(sta_all_channels{n}, sta_num_samples, width, height, fps, FLIP_XY);
    strongest_RF.channel_name = channel_names{n};
    set(gcf, 'paperposition', [0 0 24 9])
    set(gcf, 'papersize', [24 9])

    saveas(gcf, sprintf('RF_from_STA_slice_%s.png',channel_names{n}))
    saveas(gcf, sprintf('RF_from_STA_slice_%s.pdf',channel_names{n}))
    
    % saves the largest RF from each cell
    RFs{n} = strongest_RF;
    %RFs{cnt} = strongest_RF;
    %cnt = cnt + 1;
    
end



%% analyze RF centers
[ab, cd] = fit_RF_center_onto_MEA(RFs, channel_names);


%% plot mosaic

clf; hold on
for n=channel_index_to_analyze
    plot_RF(RFs{n}, FLIP_XY)  % drawing codes are moved to this function
end
xlabel('x')
ylabel('y')
title('receptive field mosaic')
axis xy
axis ([1 width+2  2 height+4])

%plot_MEA(offset_X, offset_Y)
plot_MEA_param(ab, cd);


set(gcf, 'paperposition', [0 0 24 20])
set(gcf, 'papersize', [24 20])

saveas(gcf, sprintf('mosaic.png'))
saveas(gcf, sprintf('mosaic.pdf'))



%% save RF mosaic info

save mosaic channel_names RFs width height ab cd




%% overlap all RFs with MEA as origins
clf; hold on
for n=channel_index_to_analyze
    %%
    %plot_RF_overlap(RFs{n}, FLIP_XY)  % drawing codes are moved to this function
    %
    if ~isfield(RFs{n}, 'type')
        continue
    end

    %%
    switch RFs{n}.type
        case 'ON'
            LINE_STYLE = 'r-';
            TEXT_COLOR = [1 0 0];
        case 'OFF'
            LINE_STYLE = 'b-';
            TEXT_COLOR = [0 0 1];        
    end

    xy = calc_MEA_location_from_channel_name(channel_names{n}, ab, cd)

    
    center = RFs{n}.mean-xy;

    %plot(xy(2), xy(1), 'o')
    %plot_ellipse(RFs{n}.mean, RFs{n}.cov, LINE_STYLE, FLIP_XY);
    %plot_ellipse(xy, RFs{n}.cov, LINE_STYLE, FLIP_XY);
    plot_ellipse(center, RFs{n}.cov, LINE_STYLE, FLIP_XY);
    
    % put names 
    text(center(2),center(1), channel_names{n}(4:end), 'HorizontalAlignment','center')
        
end
%plot_MEA_param(ab, cd)
title ('RF with MEA as origin')


% ROI +-6 seems to be OK

%roi_size = 6; % +- 6 pixels . TOO SLOW
roi_size = 3; % +- 3 pixels

plot([-1 1]*roi_size, [-1 -1]*roi_size, 'k--', 'linewidth',2)
plot([-1 1]*roi_size, [1 1]*roi_size, 'k--', 'linewidth',2)
plot([-1 -1]*roi_size, [-1 1]*roi_size, 'k--', 'linewidth',2)
plot([1 1]*roi_size, [-1 1]*roi_size, 'k--', 'linewidth',2)


set(gcf, 'paperposition', [0 0 24 20])
set(gcf, 'papersize', [24 20])

saveas(gcf, sprintf('RFs_overlap.png'))
saveas(gcf, sprintf('RFs_overlap.pdf'))


%return


%% calc STA & STC  # ONLY FOR SELECTED CHANNELS and inside of RF!
clear sta_ROI
clear stc_RF
num_significant_evs = [];
close all

shift_min = sta_num_samples*10;
shift_max = sta_num_samples*50;



for n = channel_index_to_analyze

    if ~isfield(RFs{n},'mean')
        continue;
    end
    
    
    %% generate RF mask 
    %     USE_MASK=0;
    USE_MASK=1;
    
    % [NEW] ROI mask is generated from the MEA position
    xy = calc_MEA_location_from_channel_name(channel_names{n}, ab, cd);   
    [mask, XX, YY] = generate_ROI_mask_from_MEA(xy, roi_size, width, height);
    
    % [OLD] mask was generated by RF
    %[mask, XX, YY] = generate_ellipse_mask(RFs{n}.mean, RFs{n}.cov, width, height);   
    %[mask, XX, YY] = generate_ellipse_mask(RFs{n}.mean, RFs{n}.cov*2^2, width, height);   % consider larger area
    
    
    % apply mask to select stim to be analyzed
    stim_chosen=stim(:,mask(:))-0.5;
    spike_train_chosen = spike_train(:,n);
    
    %% STC analysis
    % cleaned up code 
    [sta_ROI{n}, ev, u] = calc_STA_and_STC(stim_chosen(1:end-shift_max,:), spike_train_chosen(1:end-shift_max), sta_num_samples);
    
    ev = ev(ev>1e-5);
    num_non_zero_eig_val = length(ev);
    
    u = u(:,1:num_non_zero_eig_val);
    
    %% re-fit ellipse to sta_ROI
    [pos_RF, neg_RF, strongest_RF] = calc_RF_from_STA_slice(sta_ROI{n}+0.5, sta_num_samples, max(XX(mask>0))-min(XX(mask>0))+1, max(YY(mask>0))-min(YY(mask>0))+1, fps, FLIP_XY);

    if ~isempty(strongest_RF)
        ROI_xy = [min(XX(mask>0))-1 min(YY(mask>0))-1];
        strongest_RF.mean = strongest_RF.mean + ROI_xy;
    end
            

    %% identify significant eigen values using bootstraping with nested hypothesis
    disp(['Searching for significant eigenvalues of ' channel_names{n}])
    tic;
    num_repeat=50;
    idx_significant_ev = find_significant_eigen_values(ev, u, stim_chosen, spike_train_chosen, sta_num_samples, num_repeat, [shift_min shift_max], sta_ROI{n});   
    toc
    
        
    num_significant_evs(n) = length(idx_significant_ev);
    disp(sprintf('%d significant igenvalues found',num_significant_evs(n)))
    if ~isempty(idx_significant_ev)
        idx_significant_ev
    end
    
    %% plot results for this channel
    %close all
    
    figure(1)
    if num_significant_evs(n) == 0 
        r=2;
    else
        r=1+ceil(num_significant_evs(n)/2);
    end
    c=2;
    clf
    subplot(r,c,1)
    plot(reshape(sta_ROI{n}, sta_num_samples,[]))
    ylabel('STA')
    box off
    title(sprintf('STA of %s',channel_names{n}),'Interpreter', 'none')

    subplot(r,c,2)
    plot(ev)
    hold on
%     plot(ev_upper, 'r--')
%     plot(ev_lower, 'r--')
    set(gca,'yscale','log')
    axis tight
    
%     XLIM=get(gca,'xlim');
%     plot(XLIM, ev_range(1)*[1 1], 'r--')
%     plot(XLIM, ev_range(2)*[1 1], 'r--')

    plot(idx_significant_ev, ev(idx_significant_ev), '*k')
    ylabel('eigen values of STC')
    box off
    %axis tight
        
    
    drawnow
    
    if ~isempty(idx_significant_ev)
        disp(sprintf('Found significant eigen values at %d',idx_significant_ev))
        

        for i = 1:length(idx_significant_ev)
            ii = idx_significant_ev(i);
                    
            figure(1)
            subplot(r,c,2+i)
            us = reshape(u(:,ii),sta_num_samples,[]);
            plot(us)

            ylabel('STC')
            title (sprintf('STC filter %d',ii))
            box off

            
            
            figure(2)   % plot spatial pattern in a separate figure
            plot_stim_slices_with_mask(u(:,ii), sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height, FLIP_XY, xy)
            hold on;
            %plot_MEA_param(ab,cd)
            %plot(xy(2), xy(1), 'ro')
            
            set(gcf, 'paperposition', [0 0 24 9])
            set(gcf, 'papersize', [24 9])

            saveas(gcf, sprintf('STC_in_ROI_%s_eig%d.png',channel_names{n},ii))
            saveas(gcf, sprintf('STC_in_ROI_%s_eig%d.pdf',channel_names{n},ii))
            
            
            
        end
        
        figure(3) % plot STA for comparison
        plot_stim_slices_with_mask(sta_ROI{n}, sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height, FLIP_XY, xy, RFs{n})

        % plot re-fit ellipse
        if ~isempty(strongest_RF)
            subplot(2,5, strongest_RF.slice)
            plot_RF(strongest_RF, ~FLIP_XY)
        end


        set(gcf, 'paperposition', [0 0 24 9])
        set(gcf, 'papersize', [24 9])

        saveas(gcf, sprintf('STA_in_ROI_%s.png',channel_names{n}))
        saveas(gcf, sprintf('STA_in_ROI_%s.pdf',channel_names{n}))
            

        figure(1) % go back to figure 1
    else % no significant eigen values 
         % plot us in gray for checking
         
         % plot the largest eig vec in gray
        figure(1)
        subplot(r,c,3)
        us = reshape(u(:,1),sta_num_samples,[]);
        plot(us, 'color', 0.6*[1 1 1])
        
        title('eig. vec. for the largest eig. val.')        
        box off

        
         % plot the smallest eig vec in gray
        subplot(r,c,4)
        us = reshape(u(:,length(ev)),sta_num_samples,[]);
        plot(us, 'color', 0.6*[1 1 1])

        title('eig. vec. for the smallest eig. val.')        
        box off
        
        drawnow
        
        
    end


    % save ONLY WHEN significant eigen value is found 
    if ~isempty(idx_significant_ev)
        %%
        figure(1)
        
        set(gcf, 'paperposition', [0 0 24 10*r])
        set(gcf, 'papersize', [24 10*r])

        saveas(gcf, sprintf('STC_in_ROI_%s.png',channel_names{n}))
        saveas(gcf, sprintf('STC_in_ROI_%s.pdf',channel_names{n}))
    end

    
end




return 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLEANED UP UPTO THIS POINT (2019. 2. 11)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% 2. calc STA & STC  # ONLY FOR SELECTED CHANNELS!
for n = channel_index_to_analyze
    
    
    
    
    %% calc STA and STC
    [sta_all_channels{n}, stc{n}] = calc_STA_and_STC(stim, spike_train(:,n), sta_num_samples);
    
    
    %% calc range of STC eigen values (THIS IS TOO SLOW)
    disp(['Calculating the significance interval of eigen values for ' channel_names{n}])
    tic;
    ev_range = calc_STC_eigenvalue_range(stim, spike_train(:,n), sta_num_samples, 50, sta_num_samples*[10 50])
    toc;
    
    
    
    %% analyze STA (code from find_STA_with_max_var.m)
    sta_var = var(sta_all_channels{n},1);

    % find pixel with the maximum variance in STA
    [~, idx_max_STA] = max(sta_var);

    % choose STA with maximum variance
    sta_max_var = sta_all_channels{n}(:,idx_max_STA);
    
    
    % analyze nonlinearity!
    [generator, firing_rate] =calc_nonlinearity(stim, spike_train(:,n), sta_all_channels{n}, 16);
    %[generator, firing_rate] =calc_nonlinearity(stim(:,idx_max_STA), spike_train(:,n), sta_all_channels{n}(:,idx_max_STA), 16);
    
    
    
    %% plot STA
    clf
    r=2; c=4;     
    subplot(r,c,1)
    imshow(sta_all_channels{n}')
    xlabel('t')
    ylabel('pixel')
    title('STA')


    subplot(r,c,2)
    imshow(reshape(sta_var,height,width)',[])   % stim is rotated by 90 degree!
    xlabel('y') %xlabel('x')         % swap x and y
    ylabel('x') %ylabel('y')
    title(sprintf('var STA across time'))
    axis on xy

    % mark the pixel with the largest variance 
    [XX,YY]= meshgrid(1:height);  % grid for ploting

    hold on;
    %plot(XX(idxMax),YY(idxMax),'+r', 'markersize', 8, 'linewidth', 2)
    plot(YY(idx_max_STA),XX(idx_max_STA),'+r', 'markersize', 8, 'linewidth', 2)


    subplot(r,c,3)
    plot(gridT, sta_max_var,'r')
    xlabel('t')
    ylabel('STA')
    title('STA for the pixel with the largest variance')
    box off
    
    % analyze nonlinearity!
    subplot(r,c,4)
    plot(generator, firing_rate, '.k')
    xlabel('generator signal')
    ylabel('spikes / bin')
    box off
    
    
    
    
    %% analyze STC result
    %STC = stc{n}{idx_max_STA,idx_max_STA};   % choose STC of the target pixel
    %[U, D, V] = svd(STC);
    disp('anlayzin STC...')
    tic;
    [U, D, ~] = svd(stc{n});             % takes about 40 sec
    toc
    
    ev = diag(D);
    ev = ev(ev>1e-5);
    
    % match the sign of U to be consistent with sta_max_var
%     sgn = sign(U*sta_max_var);
%     UU=U*diag(sgn);
%     VV=diag(sgn)*V;
    
    %% plot
    subplot(r,c,5)
    %imagesc(STC)
    plot_off_diag(stc{n})
    %colormap gray
    box off
    title('STC for the pixel with the largest variance')
    
    subplot(r,c,6)
    plot(ev, 'ok'); hold on
    %plot(1, ev(1), 'b*')
    %plot(length(ev), ev(end), 'r*')   
    
    % plot range 
    if exist('ev_range','var')
        hold on
        XLIM=get(gca,'xlim');
        plot(XLIM, ev_range(1)*[1 1], 'r--')
        plot(XLIM, ev_range(2)*[1 1], 'r--')
    end
    
    
    box off
    title ('eigen values')
    
    subplot(r,c,7)
    %plot(eig(stc{max_idx(n),max_idx(n)}), 'o')
    %plot(UU(:,[1,end]))
    plot(U(:,[1,end]))
    box off
    title ('eigen vectors')
    
    box off
    
    
%     clf
%     plot_STA_and_STC(sta_all_channels{n}, stc{n}, gridT, width, height)
    
    
    set(gcf, 'paperposition', [0 0 20 24])
    set(gcf, 'papersize', [20 24])
    
%     set(gcf, 'paperposition', [0 0 11 5])
%     set(gcf, 'papersize', [11 5])

    saveas(gcf, sprintf('%s_STA_and_STC.pdf', channel_names{n}))
    saveas(gcf, sprintf('%s_STA_and_STC.png', channel_names{n}))
    
    %% analysis of STC for all pixels (2018.10.24)
    num_pixels = size(stc{n},1);
    is_significant_pixel = zeros(1,num_pixels);
    clear idx_significant_pixels
    clear Us
    evs = [];
    for i=1:num_pixels  % for each pixel 
        % check eiven values out of the range 
        
        [U, D, V ] = svd(stc{n}{i,i});
        ev = diag(D);
    
        % store all info
        evs = [evs ev];
        Us{i} = U; 
            
        
        idx_pos = find(ev > ev_range(2));        
        idx_neg = find(ev < ev_range(1));        
        idx_pos_significant_pixels{i} = idx_pos;
        idx_neg_significant_pixels{i} = idx_neg;
        
        if ~isempty(idx_pos) | ~isempty(idx_neg)
            is_significant_pixel(i) = true;
        end
    end
    
    num_significant_pixels = sum(is_significant_pixel);
    disp(sprintf('eigen values are significant in %d pixels in %s',num_significant_pixels,channel_names{n}))
    
    [ev_max,idx_max_STC]=max(evs(1,:));
    [ev_min,idx_min_STC]=min(evs(end,:));
    
    
    clf
    r=3;c=2;
    subplot(r,c,1)
    plot(evs(:,find(~is_significant_pixel)), ':', 'color', 0.5*[1 1 1]) 
    ylabel('eivenvalue')
    hold on
    plot(evs(:,find(is_significant_pixel)), '--.')        
    XLIM=get(gca,'xlim');
    plot(XLIM, ev_range(1)*[1 1], 'r--')
    plot(XLIM, ev_range(2)*[1 1], 'r--')
    box off
    title(sprintf('eigen values are significant in %d pixels in %s',num_significant_pixels,channel_names{n}),'Interpreter', 'none')
    
    subplot(r,c,2)
    hist(evs(:),32); hold on; 
    ylabel('eivenvalue')
    xlabel('count')
    title ('histogram of eivenvalues')
    YLIM=get(gca,'ylim');
    plot(ev_range(1)*[1 1], YLIM, 'r--')
    plot(ev_range(2)*[1 1], YLIM, 'r--')
    box off
    
    subplot(r,c,3)
    imagesc(reshape(evs(1,:),8,8)); colorbar; colormap gray; box off    
    %imagesc(reshape(evs(1,:),8,8), [ev_min ev_max]); colorbar; colormap gray; box off    
    hold on
    plot(floor((idx_max_STC-1)/8)+1, mod(idx_max_STC-1,8)+1, 'or', 'markersize', 12, 'linewidth', 3)   
    plot(floor((idx_max_STA-1)/8)+1, mod(idx_max_STA-1,8)+1, '+r', 'markersize', 8, 'linewidth', 2)
    
    
    title('largest eiven values for each pixel')
    
    subplot(r,c,4)
    imagesc(reshape(evs(end,:),8,8)); colorbar; colormap gray; box off
    %imagesc(reshape(evs(end,:),8,8), [ev_min ev_max]); colorbar; colormap gray; box off
    hold on
    plot(floor((idx_min_STC-1)/8)+1, mod(idx_min_STC-1,8)+1, 'ob', 'markersize', 12, 'linewidth', 3)   
    plot(floor((idx_max_STA-1)/8)+1, mod(idx_max_STA-1,8)+1, '+r', 'markersize', 8, 'linewidth', 2)
    
    title('smallest eiven values for each pixel')
    
    subplot(r,c,5)
    plot(Us{idx_max_STC}(:,idx_pos_significant_pixels{idx_max_STC}),'r')
    %plot(Us{idx_max_STC}(:,1),'r')
    title('significant eigen vector(s)') % of the pixel with the largest eiven value')
    set(gca,'ylim', [-1 1])
    box off
    
    subplot(r,c,6)
    plot(Us{idx_min_STC}(:,idx_neg_significant_pixels{idx_min_STC}),'b')
    %plot(Us{idx_min_STC}(:,end),'b')
    title('significant eigen vector(s)') % of the pixel with the smallest eiven value')
    set(gca,'ylim', [-1 1])
    box off
    

    set(gcf, 'paperposition', [0 0 20 24])
    set(gcf, 'papersize', [20 24])
%     set(gcf, 'paperposition', [0 0 8 10])
%     set(gcf, 'papersize', [8 10])
    
    saveas(gcf, sprintf('%s_%dHz_%s_STA_and_STC_all_pixels.pdf',CELL_TYPE, fps, channel_names{n}))
    saveas(gcf, sprintf('%s_%dHz_%s_STA_and_STC_all_pixels.png',CELL_TYPE, fps, channel_names{n}))
    

end






% saveas(gcf, sprintf('mosaic_%s.png',channel_names{n}))
% saveas(gcf, sprintf('mosaic_%s.pdf',channel_names{n}))



%% Now, let's fit GLM!



%% Set up basis for coupling filters

% Make basis for self-coupling term
% ihbasprs.ncols = 2; % number of basis vectors
% ihbasprs.hpeaks = [.1, .2]; % peak of 1st and last basis vector
% ihbasprs.ncols = 3; % number of basis vectors
% ihbasprs.hpeaks = [.1, .2, .3]; % peak of 1st and last basis vector

switch fps
    case 10   % for 10 Hz datasets        
        % param for self-coupling filters
        ihbasprs.ncols = 1; % number of basis vectors
        ihbasprs.hpeaks = 0.1; % peak of 1st and last basis vector
        
        %ihbasprs.b = .001;  % scaling (smaller -> more logarithmic scaling)
        ihbasprs.b = .02;  % scaling (smaller -> more logarithmic scaling)
        %ihbasprs.b = 0.1;  % scaling (smaller -> more logarithmic scaling)
        ihbasprs.absref = []; % absolute refractory period basis vector (optional)

        % param for cross-coupling filters
        ihbasprs2.ncols = 1;  % number of basis vectors
        ihbasprs2.hpeaks = [0.01, .1]; % put peak at 0.1s and "effective" 1st peak at 0
        %ihbasprs2.b = .001;  % smaller -> more logarithmic scaling
        %ihbasprs2.b = .01;  % smaller -> more logarithmic scaling
        ihbasprs2.b = .1;  % smaller -> more logarithmic scaling
        ihbasprs2.absref = []; % no abs-refracotry period for this one

    case 25
        % param for self-coupling filters
        ihbasprs.ncols = 1; % number of basis vectors
        ihbasprs.hpeaks = 0.04; % peak of 1st and last basis vector
        
        %ihbasprs.b = .001;  % scaling (smaller -> more logarithmic scaling)
        ihbasprs.b = .02;  % scaling (smaller -> more logarithmic scaling)
        %ihbasprs.b = 0.1;  % scaling (smaller -> more logarithmic scaling)
        ihbasprs.absref = []; % absolute refractory period basis vector (optional)

        % param for cross-coupling filters
        ihbasprs2.ncols = 1;  % number of basis vectors
        ihbasprs2.hpeaks = [0.01, 0.04*2]; % put peak at 0.1s and "effective" 1st peak at 0
        %ihbasprs2.b = .001;  % smaller -> more logarithmic scaling
        %ihbasprs2.b = .01;  % smaller -> more logarithmic scaling
        ihbasprs2.b = .1;  % smaller -> more logarithmic scaling
        ihbasprs2.absref = []; % no abs-refracotry period for this one
        
end

% Make basis for self-coupling filter
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dtSpike);
nht = length(iht); % number of bins

% Make basis for cross-coupling filter
[iht2,ihbas2,ihbasis2] = makeBasis_PostSpike(ihbasprs2,dtSpike);
nht2 = length(iht2);

% pad to put them into the same time bins
if nht2>nht
    % padd ih1 with zeros
    iht = iht2; zz = zeros(nht2-nht,ihbasprs.ncols);
    ihbas = [ihbas;zz]; ihbasis = [ihbasis;zz]; nht=nht2;
elseif nht2<nht
    % padd ih1 with zeros
    iht2 = iht; zz = zeros(nht-nht2,ihbasprs2.ncols);
    ihbas2 = [ihbas2;zz]; ihbasis2 = [ihbasis2;zz]; nht2=nht;
end    


% Set an example filter shapes

%wself = [-10; 3.2; -1]; % weights for self-coupling term
wself = [-3; 0.5];
wself = wself(1:ihbasprs.ncols);
ihself = ihbasis*wself; % self-coupling filter
wcpl = 0.5; % weights for cross-coupling term
ihcpl = ihbasis2*wcpl; % cross-coupling filter


% plot coupling filter basis
%clf; 

c = get(0, 'DefaultAxesColorOrder');

XLIM = [0 0.7];
clf
subplot(311); 
plot(iht, ihbasis, 'Color', c(1,:), 'linewidth', 2); hold on
plot(ihbasprs.hpeaks, 1, '*k', 'linewidth', 1) %, 'MarkerSize', 7);
xlabel('time after spike (s)'); title('post-spike basis'); box off
set(gca,'xlim',XLIM);
subplot(312); 
plot(iht2, ihbasis2, 'Color', c(2,:), 'linewidth', 2); hold on
plot(ihbasprs2.hpeaks, 1, '*k', 'linewidth', 1) %, 'MarkerSize', 7);
xlabel('time after spike (s)'); title('coupling basis'); box off
set(gca,'xlim',XLIM);
subplot(313); 
plot(iht, exp(ihself), iht, exp(ihcpl), iht, iht*0+1, 'k:', 'linewidth', 2);
legend('self-coupling', 'cross-coupling');
xlabel('time lag (s)');
ylabel('gain (sp/s)');
box off
set(gca,'xlim',XLIM);

set(gcf, 'paperposition', [0 0 6 9])
set(gcf, 'papersize', [6 9])
saveas(gcf, sprintf('coupling-filters-basis_%dHz.pdf',fps))
saveas(gcf, sprintf('coupling-filters-basis_%dHz.png',fps))


%% Fit GLM 
%  Initialize params for fitting --------------
% Set params for fitting, including bases 


%ggInit = makeFittingStruct_GLM(dtStim,dtSp);  % Initialize params for fitting struct 

nkt  = size(STA,1);  % number of bins for stimulus filter
nkbasis = 2
ggInit = makeFittingStruct_GLM(dtStim,dtSpike,nkt,nkbasis,STA)


% Initialize fields for coupling filters (using h bases computed above)
ggInit.ihbas = ihbas; % h self-coupling basis
ggInit.ihbas2 = ihbas2; % h coupling-filter basis
nhbasis = size(ihbas,2); % number of basis vectors in h basis
nhbasis2 = size(ihbas2,2); % number of basis vectors in h basis
ggInit.ihw = zeros(nhbasis,1); % init params for self-coupling filter
ggInit.ihw2 = zeros(nhbasis2,N-1); % init params for cross-coupling filter
ggInit.ih = [ggInit.ihbas*ggInit.ihw ggInit.ihbas2*ggInit.ihw2];
ggInit.iht = iht;
ggInit.dc = 0; % Initialize dc term to zero
ggInit.sps = spike_train(:,1); % spikes from 1 cell
ggInit.sps2 = spike_train(:,2:N); % spikes from all 
ggInit.couplednums = 2:N; % cell numbers of cells coupled to this one 


ggfit(1:N) = ggInit; % initialize fitting struct
kfit = zeros(nkt,64,N);
dcfit = zeros(N,1);
for jj = 1:N
    couplednums = setdiff(1:N,jj);  % cell numbers of cells coupled to this one (jj)

    % Set spike responses for this cell (jj) and coupled cells
    ggInit.sps = spike_train(:,jj);
    ggInit.sps2 = spike_train(:,couplednums);
    ggInit.couplednums = couplednums; % numbers of cells coupled to this one (for clarity)

    % Do ML fitting
    fprintf('==== Fitting filters to neuron %d ==== \n',  jj)
    opts = {'display', 'iter', 'maxiter', 1000};
    ggfit(jj) = MLfit_GLM(ggInit,stim,opts); % do ML fitting
    kfit(:,:,jj) = ggfit(jj).k;
    dcfit(jj) = ggfit(jj).dc;    
end


%% plot linear filter of the last cell (for check)
clf;colormap gray
for n=1:N
    %%
    subplot(2,N,n);imagesc(STAall{n}');  ylabel('pixel'); xlabel('time bin');
    title(sprintf('STA (%s)', channel_names{n}),'Interpreter', 'none')
    subplot(2,N,n+N);imagesc(ggfit(n).k'); ylabel('pixel'); xlabel('time bin')
    title(sprintf('k by GLM (%s)', channel_names{n}),'Interpreter', 'none')
end
set(gcf, 'paperposition', [0 0 8 8])
set(gcf, 'papersize', [8 8])
saveas(gcf, sprintf('%scell_%dHz_STA_vs_GLM.pdf', CELL_TYPE,fps))


%% plot coupling filters
clf reset; 
colors = get(gca,'colororder');
ncolrs = min(size(colors,1),N); % number of traces to plot for each cell 

% find max of y range
ymax = 0;
for n=1:N
    ymax = max([exp(ggfit(n).ih(:)); ymax]); 
end
ymax = min(ymax,10);



switch N
    case 7 
        num_row=3;num_col=3;
    otherwise
        num_row=1;num_col=N;
end
    
for jj = 1:N
    subplot(num_row,num_col,jj); 
    plot(ggfit(jj).iht, exp(ggfit(jj).ih(:,1)), '--k', 'linewidth', 2); hold on
    plot(ggfit(jj).iht, exp(ggfit(jj).ih(:,2:ncolrs)), '-', 'linewidth', 2); % estim
    title(sprintf(' cell %s',channel_names{jj}(4:end) )); axis tight; set(gca,'ylim',[0,ymax]);
    box off
    
    if jj==1 || jj==4 || jj ==7
        ylabel('gain (sp/s)'); 
    end
    xlabel('time after spike (s)');
end


set(gcf, 'paperposition', [0 0 3*num_col 3*num_row])
set(gcf, 'papersize', [3*num_col 3*num_row])
saveas(gcf, sprintf('%scell_%dHz_coupling_filters.pdf', CELL_TYPE,fps))

