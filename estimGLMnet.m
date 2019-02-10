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
    
    channel_name = channel_names{n}
    spike_time=eval(channel_name);
    
    
    % let's check for each stim time bin
    for i = 2:binStim
        t0 = A1a(i-1);
        t1 = A1a(i);

        % find stim time that occured during t0 and t1
        idx = find(spike_time>t0 & spike_time<=t1);

        if ~isempty(idx)
            %disp('found')
            spike_train(i,n) = 1;
        end

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
tic
for n=1:length(channel_names)
    
    disp(channel_names{n})
    
    % calc STA and STC
    [sta_all_channels{n}] = calc_STA_and_STC(stim, spike_train(:,n), sta_num_samples);
    
end
toc


%% lanch app to analyze individual channels and 
% decide which channels to analyze further

addpath tools
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



%% calc and plot RF from STA (2019. 1. 27)
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
    
    
    [pos_RF, neg_RF, strongest_RF] = calc_RF_from_STA_slice(sta_all_channels{n}, height, width, fps);
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



%% plot mosaic

clf; hold on
for n=channel_index_to_analyze
    switch RFs{n}.type
        case 'ON'
            plot_ellipse(RFs{n}.mean, RFs{n}.cov, 'r-');
            tt=text(RFs{n}.mean(1), RFs{n}.mean(2), RFs{n}.channel_name(4:end), 'HorizontalAlignment','center');
            tt.Color = [1 0 0];
        case 'OFF'
            plot_ellipse(RFs{n}.mean, RFs{n}.cov, 'b-');
            tt=text(RFs{n}.mean(1), RFs{n}.mean(2), RFs{n}.channel_name(4:end), 'HorizontalAlignment','center');
            tt.Color = [0 0 1];
    end
end
xlabel('x')
ylabel('y')
title('Receptive field mosaic')
axis image
axis ([1 width 1 height])

set(gcf, 'paperposition', [0 0 24 20])
set(gcf, 'papersize', [24 20])

saveas(gcf, sprintf('mosaic.png'))
saveas(gcf, sprintf('mosaic.pdf'))








%% calc STA & STC  # ONLY FOR SELECTED CHANNELS and inside of RF!
clear stc_RF
for n = channel_index_to_analyze
    
    close all
    
    %% generate RF mask 
    [mask, XX, YY] = generate_ellipse_mask(RFs{n}.mean, RFs{n}.cov, width, height);   
    
    %% calc STA and STC
    [sta_RF{n}, stc_RF{n}] = calc_STA_and_STC(stim(:,mask(:)), spike_train(:,n), sta_num_samples);
    

    tic;
    [u, d, ~] = svd(stc_RF{n});
    toc
    ev = diag(d);
    ev = ev(ev>1e-5);
    

    %% calc range of STC eigen values
    disp(['Calculating the significance interval of eigen values for ' channel_names{n}])
    tic;
    ev_range = calc_STC_eigenvalue_range(stim(:,mask(:)), spike_train(:,n), sta_num_samples, 50, sta_num_samples*[10 50])
    toc;
    
    idx_large_eig = find (ev > ev_range(2));
    idx_small_eig = find (ev < ev_range(1));
    
    
    %% to plot
    sta_RF_var = var(sta_RF{n},[],1);
    [~, sorted_index] = sort(sta_RF_var, 'descend');
    
    %
    clf
    subplot(231)
    plot(sta_RF{n})
    ylabel('STA')
    box off

    subplot(232)
    bar(sta_RF_var)
    xlabel('pixel index')
    ylabel('var(STA)')
    box off

    subplot(234)
    plot(ev)
    hold on
    XLIM=get(gca,'xlim');
    plot(XLIM, ev_range(1)*[1 1], 'r--')
    plot(XLIM, ev_range(2)*[1 1], 'r--')
    
    plot(idx_large_eig, ev(idx_large_eig), '*r')
    plot(idx_small_eig, ev(idx_small_eig), '*b')
    ylabel('eigen values')
    box off
    %axis tight

    
    if ~isempty(idx_large_eig)
        subplot(235)
        %idx_to_plot = 1;
        us = reshape(u(:,idx_large_eig),sta_num_samples,[]);
        plot(us)
        us_var = var(us, [], 1);
        
        %RF_area_in_pixel = size(us,2);
       % us_var = var(us, [], 1);
        %h = plot(us(:,sorted_index) - repmat(0.3*(1:RF_area_in_pixel),sta_num_samples,1), 'r');
        %h(max_idx).LineStyle = '--';
        %plot(reshape(U(:,1),sta_num_samples,[]))
        %axis off
        
        title ('STC filter with the large eig. val.')
        box off
        
        for ii = idx_large_eig
            figure
            plot_stim_slices_with_mask(u(:,ii), sta_num_samples, XX(mask(:)>0), YY(mask(:)>0))
            
            set(gcf, 'paperposition', [0 0 24 9])
            set(gcf, 'papersize', [24 9])

            saveas(gcf, sprintf('STC_inside_RF_%s_eig%d.png',channel_names{n},ii))
            saveas(gcf, sprintf('STC_inside_RF_%s_eig%d.pdf',channel_names{n},ii))
        end
        
        figure(1)
        
    end

    if ~isempty(idx_small_eig)
        subplot(236)
        %idx_to_plot = length(ev);
        us = reshape(u(:,idx_small_eig),sta_num_samples,[]);
        plot(us)
        us_var = var(us, [], 1);
        
        %RF_area_in_pixel = size(us,2);
        %h = plot(us(:,sorted_index) - repmat(0.3*(1:RF_area_in_pixel),sta_num_samples,1), 'b');
        %plot(reshape(us,sta_num_samples,[]))      
        %axis off
        title ('STC filter with the small eig. val.')
        box off
        
        
        for ii = idx_small_eig
            figure
            plot_stim_slices_with_mask(u(:,ii), sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height)
            
            set(gcf, 'paperposition', [0 0 24 9])
            set(gcf, 'papersize', [24 9])

            saveas(gcf, sprintf('STC_inside_RF_%s_eig%d.png',channel_names{n},ii))
            saveas(gcf, sprintf('STC_inside_RF_%s_eig%d.pdf',channel_names{n},ii))
        end
        
        figure(1)
    end
    
%     subplot(236)
% %     plot3(repmat((1:sta_num_samples)',1,25), repmat(xx',sta_num_samples,1), repmat(yy',sta_num_samples,1) + reshape(U(:,1),sta_num_samples,[])) 
%     bar(us_var)
% %     xlabel('t')
% %     ylabel('x')
% %     zlabel('y')

    set(gcf, 'paperposition', [0 0 24 15])
    set(gcf, 'papersize', [24 15])

    saveas(gcf, sprintf('STC_inside_RF_%s.png',channel_names{n}))
    saveas(gcf, sprintf('STC_inside_RF_%s.pdf',channel_names{n}))
    
%     %% plot filters separately 
%     figure(2)
%     plot_stim_slices(u(:,1), sta_num_samples)
%     
%     figure(3)
%     plot_stim_slices(u(:,length(ev)), sta_num_samples)
    
end







% % %% compare variances
% % figure()
% % plot(sta_RF_var, us_var, 'o')
% % xlabel('var(STA)')
% % ylabel('var(STC)')
% % box off



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

