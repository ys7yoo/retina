% script to analyze coupling effects
% data are stored under coupling_data/

%% set path to GLM tools

% This is where the GLMspiketools is installed
% You can get the package from https://github.com/ys7yoo/GLMspiketools
basedir = [getenv('HOME') '/src/GLMspiketools']

% Add a bunch sub-directories (with absoluate path names)
addpath([basedir '/glmtools_fitting/']);
addpath([basedir '/glmtools_misc/']);
addpath([basedir '/nlfuns/']);
addpath([basedir '/glmtools_spline/']);

%% load data (stim and spikes)
clear
addpath tools


% exp_param.exp_date = exp_date_str;
% data_folder_name = fullfile(base_folder_name, exp_date_str);


% base_folder_name = 'coupling_data/20180524'
% base_folder_name = 'coupling_data/20180618'
base_folder_name = 'coupling_data/20180828'
% base_folder_name = 'coupling_data/20180905' % response is too weak
% base_folder_name = 'coupling_data/20181004' % response is too weak

[stim, spike_train, channel_names, exp_param] = load_data(base_folder_name);

exp_param.num_electrodes_per_dim = 8;
exp_param.inter_electrode_space = 200 % manually set (100 um)
%exp_param.electrod_diameter = 30;


%sampling_rate, width, height
sampling_rate = exp_param.sampling_rate;
width = exp_param.num_pixels;
height = exp_param.num_pixels;

%% Run STA
%% Let's calc STA and STC
close all
%clear STA_max_var, max_idx

% sta params
%sta_num_samples = 16;
sta_num_samples = 10; % for 30Hz
% gridT = (-sta_num_samples+1:0)/sampling_rate;

%sta_num_samples = ceil(sampling_rate*0.66);  % consider 660 ms prior to spike
%sta_num_samples = floor(sampling_rate*1);  % consider 660 ms prior to spike

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

fps = sampling_rate;
sta_app

%% calc and plot RF from STA
RFs = calc_and_plot_RF_from_STA(sta_all_channels, sta_num_samples, width, height, sampling_rate, channel_names)



% %% analyze RF centers
% [ab, cd] = fit_RF_center_onto_MEA(RFs, channel_names);


%% plot mosaic
FLIP_XY=true;
channel_index_to_analyze = 1:length(channel_names)

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
%calc_MEA_loca
if exist('ab')
    plot_MEA_param(ab, cd);
end

set(gcf, 'paperposition', [0 0 24 20])
set(gcf, 'papersize', [24 20])

saveas(gcf, sprintf('mosaic.png'))
saveas(gcf, sprintf('mosaic.pdf'))


%% overlap all RFs with MEA as origins
if exist('ab')
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

        xy = calc_MEA_location_from_channel_name(channel_names{n}, ab, cd);


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
end


%% Before doing further analysis, let's compare RFs based on STA
compare_RFs


% % %% calc STA & STC  # ONLY FOR SELECTED CHANNELS and inside of RF!
% % clear sta_ROI
% % clear stc_RF
% % num_significant_evs = [];
% % close all
% % 
% % shift_min = sta_num_samples*10;
% % shift_max = sta_num_samples*50;
% % 
% % 
% % 
% % for n = channel_index_to_analyze
% % 
% %     if ~isfield(RFs{n},'mean')
% %         continue;
% %     end
% %     
% %     
% %     %% generate RF mask 
% %     USE_MASK=0;
% %     %USE_MASK=1;
% % 
% %     if USE_MASK==0
% %         mask = ones(size(stim,2),1);
% %         [XX, YY] = meshgrid(1:width, 1:height);
% %     else
% %     
% %         % [NEW] ROI mask is generated from the MEA position
% %         xy = calc_MEA_location_from_channel_name(channel_names{n}, ab, cd);   
% %         [mask, XX, YY] = generate_ROI_mask_from_MEA(xy, roi_size, width, height);
% % 
% %         % [OLD] mask was generated by RF
% %         %[mask, XX, YY] = generate_ellipse_mask(RFs{n}.mean, RFs{n}.cov, width, height);   
% %         %[mask, XX, YY] = generate_ellipse_mask(RFs{n}.mean, RFs{n}.cov*2^2, width, height);   % consider larger area
% %     end
% %     
% %     
% %     % apply mask to select stim to be analyzed
% %     stim_chosen=stim(:,mask(:))-0.5;
% %     spike_train_chosen = spike_train(:,n);
% %     
% %     %% STC analysis
% %     % cleaned up code 
% %     [sta_ROI{n}, ev, u] = calc_STA_and_STC(stim_chosen(1:end-shift_max,:), spike_train_chosen(1:end-shift_max), sta_num_samples);
% %     
% %     ev = ev(ev>1e-5);
% %     num_non_zero_eig_val = length(ev);
% %     
% %     u = u(:,1:num_non_zero_eig_val);
% %     
% %     %% re-fit ellipse to sta_ROI
% %     [pos_RF, neg_RF, strongest_RF] = calc_RF_from_STA_slice(sta_ROI{n}+0.5, sta_num_samples, max(XX(mask>0))-min(XX(mask>0))+1, max(YY(mask>0))-min(YY(mask>0))+1, fps, FLIP_XY);
% % 
% %     if ~isempty(strongest_RF)
% %         ROI_xy = [min(XX(mask>0))-1 min(YY(mask>0))-1];
% %         strongest_RF.mean = strongest_RF.mean + ROI_xy;
% %     end
% %             
% % 
% %     %% identify significant eigen values using bootstraping with nested hypothesis
% %     disp(['Searching for significant eigenvalues of ' channel_names{n}])
% %     tic;
% %     num_repeat=50;
% %     idx_significant_ev = find_significant_eigen_values(ev, u, stim_chosen, spike_train_chosen, sta_num_samples, num_repeat, [shift_min shift_max], sta_ROI{n});   
% %     toc
% %     
% %         
% %     num_significant_evs(n) = length(idx_significant_ev);
% %     disp(sprintf('%d significant igenvalues found',num_significant_evs(n)))
% %     if ~isempty(idx_significant_ev)
% %         idx_significant_ev
% %     end
% %     
% %     %% plot results for this channel
% %     %close all
% %     
% %     figure(1)
% %     if num_significant_evs(n) == 0 
% %         r=2;
% %     else
% %         r=1+ceil(num_significant_evs(n)/2);
% %     end
% %     c=2;
% %     clf
% %     subplot(r,c,1)
% %     plot(reshape(sta_ROI{n}, sta_num_samples,[]))
% %     ylabel('STA')
% %     box off
% %     title(sprintf('STA of %s',channel_names{n}),'Interpreter', 'none')
% % 
% %     subplot(r,c,2)
% %     plot(ev)
% %     hold on
% % %     plot(ev_upper, 'r--')
% % %     plot(ev_lower, 'r--')
% %     set(gca,'yscale','log')
% %     axis tight
% %     
% % %     XLIM=get(gca,'xlim');
% % %     plot(XLIM, ev_range(1)*[1 1], 'r--')
% % %     plot(XLIM, ev_range(2)*[1 1], 'r--')
% % 
% %     plot(idx_significant_ev, ev(idx_significant_ev), '*k')
% %     ylabel('eigen values of STC')
% %     box off
% %     %axis tight
% %         
% %     
% %     drawnow
% %     
% %     if ~isempty(idx_significant_ev)
% %         disp(sprintf('Found significant eigen values at %d',idx_significant_ev))
% %         
% % 
% %         for i = 1:length(idx_significant_ev)
% %             ii = idx_significant_ev(i);
% %                     
% %             figure(1)
% %             subplot(r,c,2+i)
% %             us = reshape(u(:,ii),sta_num_samples,[]);
% %             plot(us)
% % 
% %             ylabel('STC')
% %             title (sprintf('STC filter %d',ii))
% %             box off
% % 
% %             
% %             
% %             figure(2)   % plot spatial pattern in a separate figure
% %             plot_stim_slices_with_mask(u(:,ii), sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height, FLIP_XY, xy)
% %             hold on;
% %             %plot_MEA_param(ab,cd)
% %             %plot(xy(2), xy(1), 'ro')
% %             
% %             set(gcf, 'paperposition', [0 0 24 9])
% %             set(gcf, 'papersize', [24 9])
% % 
% %             saveas(gcf, sprintf('STC_in_ROI_%s_eig%d.png',channel_names{n},ii))
% %             saveas(gcf, sprintf('STC_in_ROI_%s_eig%d.pdf',channel_names{n},ii))
% %             
% %             
% %             
% %         end
% %         
% %         figure(3) % plot STA for comparison
% %         plot_stim_slices_with_mask(sta_ROI{n}, sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height, FLIP_XY, xy, RFs{n})
% % 
% %         % plot re-fit ellipse
% %         if ~isempty(strongest_RF)
% %             subplot(2,5, strongest_RF.slice)
% %             plot_RF(strongest_RF, ~FLIP_XY)
% %         end
% % 
% % 
% %         set(gcf, 'paperposition', [0 0 24 9])
% %         set(gcf, 'papersize', [24 9])
% % 
% %         saveas(gcf, sprintf('STA_in_ROI_%s.png',channel_names{n}))
% %         saveas(gcf, sprintf('STA_in_ROI_%s.pdf',channel_names{n}))
% %             
% % 
% %         figure(1) % go back to figure 1
% %     else % no significant eigen values 
% %          % plot us in gray for checking
% %          
% %          % plot the largest eig vec in gray
% %         figure(1)
% %         subplot(r,c,3)
% %         us = reshape(u(:,1),sta_num_samples,[]);
% %         plot(us, 'color', 0.6*[1 1 1])
% %         
% %         title('eig. vec. for the largest eig. val.')        
% %         box off
% % 
% %         
% %          % plot the smallest eig vec in gray
% %         subplot(r,c,4)
% %         us = reshape(u(:,length(ev)),sta_num_samples,[]);
% %         plot(us, 'color', 0.6*[1 1 1])
% % 
% %         title('eig. vec. for the smallest eig. val.')        
% %         box off
% %         
% %         drawnow
% %         
% %         
% %     end
% % 
% % 
% %     % save ONLY WHEN significant eigen value is found 
% %     if ~isempty(idx_significant_ev)
% %         %%
% %         figure(1)
% %         
% %         set(gcf, 'paperposition', [0 0 24 10*r])
% %         set(gcf, 'papersize', [24 10*r])
% % 
% %         saveas(gcf, sprintf('STC_in_ROI_%s.png',channel_names{n}))
% %         saveas(gcf, sprintf('STC_in_ROI_%s.pdf',channel_names{n}))
% %     end
% % 
% %     
% % end
% % 
% % 
% % %% plot numbers of significant eivenvalues
% % close all
% % bar(num_significant_evs); xlabel('channel index'); ylabel('number of significant eigen values'); box off
% % 
% % title(sprintf('Sigificant eigen values found in %d channels.', sum(num_significant_evs>0)))
% % 
% % set(gca,'ytick',[0 1 2 3])
% % 
% % set(gcf, 'paperposition', [0 0 12 7])
% % set(gcf, 'papersize', [12 7])
% % 
% % saveas(gcf, 'num_significatn_eigen_values.pdf')
% % saveas(gcf, 'num_significatn_eigen_values.png')
% % 
% % 
% % %% 
% % channels_with_significant_eigen_values = channel_names(num_significant_evs>0)
% % 
% % save channels_with_significant_eigen_values channels_with_significant_eigen_values






%% Set up basis for coupling filters

% Make basis for self-coupling term
% ihbasprs.ncols = 2; % number of basis vectors
% ihbasprs.hpeaks = [.1, .2]; % peak of 1st and last basis vector
% ihbasprs.ncols = 3; % number of basis vectors
% ihbasprs.hpeaks = [.1, .2, .3]; % peak of 1st and last basis vector

dtStim = 1/sampling_rate;
dtSpike = 1/sampling_rate;

switch sampling_rate
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
        
     case 30
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
        
    otherwise
        error('define params')
        
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

set(gcf, 'paperposition', [0 0 12 18])
set(gcf, 'papersize', [12 18])
saveas(gcf, sprintf('coupling-filters-basis_%dHz.pdf',fps))
saveas(gcf, sprintf('coupling-filters-basis_%dHz.png',fps))


%% Fit GLM 
%  Initialize params for fitting --------------
% Set params for fitting, including bases 

N = length(channel_names);

%ggInit = makeFittingStruct_GLM(dtStim,dtSp);  % Initialize params for fitting struct 

%nkt  = size(STA,1);  % number of bins for stimulus filter
%nkt = 10;
nkt = sta_num_samples;
nkbasis = 2
ggInit = makeFittingStruct_GLM(dtStim,dtSpike,nkt,nkbasis,reshape(sta_all_channels{1}, [sta_num_samples, width*height]))


% Initialize fields for coupling filters (using h bases computed above)
ggInit.ihbas = ihbas; % h self-coupling basis
ggInit.ihbas2 = ihbas2; % h coupling-filter basis
nhbasis = size(ihbas,2); % number of basis vectors in h basis
nhbasis2 = size(ihbas2,2); % number of basis vectors in h basis
ggInit.ihw = zeros(nhbasis,1); % init params for self-coupling filter
ggInit.ihw2 = zeros(nhbasis2,N-1); % init params for cross-coupling filter
ggInit.ih = [ggInit.ihbas*ggInit.ihw ggInit.ihbas2*ggInit.ihw2];
ggInit.iht = iht;
%ggInit.dc = 0; % Initialize dc term to zero
ggInit.dc = 0.5; % Initialize dc term to zero
ggInit.sps = spike_train(:,1); % spikes from 1 cell
ggInit.sps2 = spike_train(:,2:N); % spikes from all 
ggInit.couplednums = 2:N; % cell numbers of cells coupled to this one 


ggfit(1:N) = ggInit; % initialize fitting struct
kfit = zeros(nkt,width*height,N);
dcfit = zeros(N,1);
for jj = 1:N
    couplednums = setdiff(1:N,jj);  % cell numbers of cells coupled to this one (jj)

    ggInit.k = reshape(sta_all_channels{jj}, nkt, [])-0.5;
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


%% plot linear filter estimgaed by GLM (in ggfit(n).k)
close all

for n=1:N
    figure(n)
    
    calc_RF_from_STA_slice(ggfit(n).k, sta_num_samples, width, height, sampling_rate, FLIP_XY)
    
    set(gcf, 'paperposition', [0 0 24 9])
    set(gcf, 'papersize', [24 9])

    saveas(gcf, sprintf('RF_from_GLM_%s.png',channel_names{n}))
    saveas(gcf, sprintf('RF_from_GLM_%s.pdf',channel_names{n}))
    
end
    

% clf;colormap gray
% for n=1:N
%     %%
%     %subplot(2,N,n);imagesc(STAall{n}');  ylabel('pixel'); xlabel('time bin');
% %     title(sprintf('STA (%s)', channel_names{n}),'Interpreter', 'none')
%     subplot(2,N,n+N);imagesc(ggfit(n).k'); ylabel('pixel'); xlabel('time bin')
%     title(sprintf('k by GLM (%s)', channel_names{n}),'Interpreter', 'none')
% end
% set(gcf, 'paperposition', [0 0 8 8])
% set(gcf, 'papersize', [8 8])
% saveas(gcf, sprintf('%s_%dHz_STA_vs_GLM.pdf', exp_date,fps))


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


set(gcf, 'paperposition', [0 0 6*num_col 6*num_row])
set(gcf, 'papersize', [6*num_col 6*num_row])
saveas(gcf, sprintf('coupling_filters_%dHz_.pdf', sampling_rate))
