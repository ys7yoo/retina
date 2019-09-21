%% Demo for calculating STA and STC

% set path to helper functions

addpath tools

%% load data! 

clear

% unified loaindg routine
%base_folder_name = 'data/20180618'
%base_folder_name = 'data/20180621'
base_folder_name = 'data/20180626'
%base_folder_name = 'data/20180905' # 26x26
%base_folder_name = 'data/20180828_ReceptiveField_TemporalFilter'

[stim, spike_train, channel_names, exp_param] = load_data(base_folder_name);

exp_param.num_electrodes_per_dim = 8;
exp_param.inter_electrode_space = 200 % manually set (100 um)
%exp_param.electrod_diameter = 30;


%sampling_rate, width, height
sampling_rate = exp_param.sampling_rate;
width = exp_param.num_pixels_per_dim;
height = exp_param.num_pixels_per_dim;


%% Set up parameters for STA and STC
close all

% sta params
%sta_num_samples = 16;
sta_num_samples = 10; % for 10Hz

gridT = (-sta_num_samples+1:0)/sampling_rate;       % time steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. STA analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calc STA for all the channels
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

%% analyze all ch
channel_index_to_analyze = 1:length(channel_names)

%% choose channels to analyze
% selected_channel_names = {'45a'}
% %selected_channel_names = {'41a', '54e', '72a'}
% %selected_channel_names = {'13c', '37c', '58b', '75d', '77d', '16b', '26a', '27a', '37b', '47b'}
% 
% channel_index_to_analyze = calc_channel_index(channel_names, selected_channel_names)


%% fit elliptical RFs to STA results (2019. 1. 27)

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
    [pos_RF, neg_RF, strongest_RF] = calc_RF_from_STA_slice(sta_all_channels{n}, sta_num_samples, width, height, sampling_rate, FLIP_XY);
    strongest_RF.channel_name = channel_names{n};
    set(gcf, 'paperposition', [0 0 24 9])
    set(gcf, 'papersize', [24 9])

    subplot(2,5, strongest_RF.peak_center(1));
    plot(strongest_RF.peak_center(3), strongest_RF.peak_center(2), '+c')
    
    saveas(gcf, sprintf('RF_from_STA_slice_%s.png',channel_names{n}))
    saveas(gcf, sprintf('RF_from_STA_slice_%s.pdf',channel_names{n}))
    
    % saves the largest RF from each cell
    RFs{n} = strongest_RF;
    %RFs{cnt} = strongest_RF;
    %cnt = cnt + 1;
    
end



%% analyze RF centers for multiple channels 

if length(channel_index_to_analyze) > 1
    [ab, cd] = fit_RF_center_onto_MEA_general(RFs, channel_names);
    %[ab, cd] = fit_RF_center_onto_MEA(RFs, channel_names);  % for fixed scale

    if isinf(ab(1))
        disp('Cannot fit with a single channel. Put ROI manually!')
    end



    set(gcf, 'paperposition', [0 0 24 20])
    set(gcf, 'papersize', [24 20])

    saveas(gcf, sprintf('fit_RF_onto_MEA.png'))
    saveas(gcf, sprintf('fit_RF_onto_MEA.pdf'))


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

        xy = calc_MEA_location_from_channel_name(channel_names{n}, ab, cd);

        if isfield(RFs{n}, 'mean')
            center = RFs{n}.mean-xy;
            
            %plot(xy(2), xy(1), 'o')
            %plot_ellipse(RFs{n}.mean, RFs{n}.cov, LINE_STYLE, FLIP_XY);
            %plot_ellipse(xy, RFs{n}.cov, LINE_STYLE, FLIP_XY);
            plot_ellipse(center, RFs{n}.cov, LINE_STYLE, FLIP_XY);

        else
            center = RFs{n}.peak_center(2:3)-xy;
            
            if FLIP_XY
                plot(center(2), center(1), '+g');
            else
                plot(center(1), center(2), '+g');
            end
            
        end

        % put names 
        text(center(2),center(1), channel_names{n}(4:end), 'HorizontalAlignment','center')

    end
    axis xy
    %plot_MEA_param(ab, cd)

    % ROI +-6 seems to be OK

    %roi_size = 6; % +- 6 pixels . TOO SLOW
    %roi_size = 3; % +- 3 pixels
    %roi_size = mean([ab(1); -cd(1)]);       

    % measured spacing between electrodes
    inter_electrode_space_x=-cd(1)
    inter_electrode_space_y=ab(1)

    plot([-inter_electrode_space_y inter_electrode_space_y], [-inter_electrode_space_x -inter_electrode_space_x], 'k--', 'linewidth',2)
    plot([-inter_electrode_space_y inter_electrode_space_y], [inter_electrode_space_x inter_electrode_space_x], 'k--', 'linewidth',2)
    plot([-inter_electrode_space_y -inter_electrode_space_y], [-inter_electrode_space_x inter_electrode_space_x], 'k--', 'linewidth',2)
    plot([inter_electrode_space_y inter_electrode_space_y], [-inter_electrode_space_x inter_electrode_space_x], 'k--', 'linewidth',2)

    title (sprintf('RF with MEA as origin (inter-electrode space in pixels=(%.1f,%.1f)',inter_electrode_space_y,inter_electrode_space_x))

    set(gcf, 'paperposition', [0 0 24 20])
    set(gcf, 'papersize', [24 20])

    xlabel('x')
    ylabel('y')

    saveas(gcf, sprintf('RFs_overlap.png'))
    saveas(gcf, sprintf('RFs_overlap.pdf'))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. calc STC  # ONLY FOR SELECTED CHANNELS and inside of RF!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear sta_ROI
clear stc_RF
num_significant_evs = [];
clear stc_ev
clear stc_u
clear stc_significant_ev_idx
close all

shift_min = sta_num_samples*10;
shift_max = sta_num_samples*50;


%% set up ROI mask
USE_ROI_MASK=1;
    
% ROI mask based on the center of RF and radius
roi_radius = 1;  % for 3x3 mask
mask_name = '3x3';
    
% roi_radius = 2; % for 5x5 mask
% mask_name = '5x5';
    

for n = channel_index_to_analyze

%     if ~isfield(RFs{n},'mean')
%         continue;
%     end

    %% choose center 
    center = RFs{n}.peak_center(2:3);
    % center = RFs{n}.mean;
    
    %% set up ROI mask
    if USE_ROI_MASK
        [mask, XX, YY] = generate_ROI_mask_neighbors(width, height, center, roi_radius);
    end

    % Alternative ways to set ROI mask
    % [Option 1] mask was generated by RF
    %[mask, XX, YY] = generate_ellipse_mask(RFs{n}.mean, RFs{n}.cov, width, height);   
    %[mask, XX, YY] = generate_ellipse_mask(RFs{n}.mean, RFs{n}.cov*2^2, width, height);   % consider larger area
    
    % [Option 2] ROI mask is generated from the MEA position
%     xy = calc_MEA_location_from_channel_name(channel_names{n}, ab, cd);   
%     [mask, XX, YY] = generate_ROI_mask_from_MEA(xy, roi_size, width, height);
    
    if USE_ROI_MASK
        % apply mask to select stim to be analyzed
        stim_chosen=stim(:,mask(:))-0.5;
        spike_train_chosen = spike_train(:,n);
    end
    
    %% STC analysis
    % cleaned up code 
    [sta_ROI{n}, ev, u, ~, score] = calc_STA_and_STC(stim_chosen(1:end-shift_max,:), spike_train_chosen(1:end-shift_max), sta_num_samples, true, channel_names{n});
    
    ev = ev(ev>1e-5);
    num_non_zero_eig_val = length(ev);
    
    u = u(:,1:num_non_zero_eig_val);
    
%     %% re-fit ellipse to sta_ROI
%     [pos_RF, neg_RF, strongest_RF] = calc_RF_from_STA_slice(sta_ROI{n}+0.5, sta_num_samples, max(XX(mask>0))-min(XX(mask>0))+1, max(YY(mask>0))-min(YY(mask>0))+1, sampling_rate, FLIP_XY);
% 
%     if ~isempty(strongest_RF)
%         ROI_xy = [min(XX(mask>0))-1 min(YY(mask>0))-1];
%         strongest_RF.mean = strongest_RF.mean + ROI_xy;
%     end
            

    %% identify significant eigen values using bootstraping with nested hypothesis
    disp(['Searching for significant eigenvalues of ' channel_names{n}])
    tic;
    num_repeat=50;
    idx_significant_ev = find_significant_eigen_values(ev, u, stim_chosen, spike_train_chosen, sta_num_samples, num_repeat, [shift_min shift_max], sta_ROI{n});   
    toc
    
    idx_significant_ev = sort(idx_significant_ev);
    
    num_significant_evs(n) = length(idx_significant_ev);
    disp(sprintf('%d significant igenvalues found',num_significant_evs(n)))
    if ~isempty(idx_significant_ev)
        idx_significant_ev
    end
    
    %% SAVE STC RESULTS FOR ALL THE CHANNELS
    stc_ev{n} = ev;
    stc_u{n} = u;
    stc_significant_ev_idx{n} = idx_significant_ev;
    stc_score{n} = score;
    
    %% plot results for this channel
    %close all
    
    figure(1)
    if num_significant_evs(n) == 0 
        r=2;
    else
        %r=1+ceil(num_significant_evs(n)/2);
        r=1+num_significant_evs(n);
    end
    c=2;
    clf
    subplot(r,c,1)
    plot(gridT, reshape(sta_ROI{n}, sta_num_samples,[]))
    ylabel('STA')
    xlabel('time to spike (s)')
    box off
    title(sprintf('%s',channel_names{n}),'Interpreter', 'none')

    subplot(r,c,2)
    plot(ev, 'ko:', 'markersize', 5)
    hold on
%     plot(ev_upper, 'r--')
%     plot(ev_lower, 'r--')
%     set(gca,'yscale','log')
    axis tight
    
%     XLIM=get(gca,'xlim');
%     plot(XLIM, ev_range(1)*[1 1], 'r--')
%     plot(XLIM, ev_range(2)*[1 1], 'r--')

    plot(idx_significant_ev, ev(idx_significant_ev), '*r')
    ylabel('eigen values of STC')
    box off
    %axis tight
    
    drawnow
    
    if ~isempty(idx_significant_ev)
        disp(sprintf('Found significant eigen values at %d',idx_significant_ev))

        for i = 1:length(idx_significant_ev)
            ii = idx_significant_ev(i);
                    
            %%
            figure(1)
            subplot(r,c,i*2+1)
            us = reshape(u(:,ii),sta_num_samples,[]);
            plot(gridT,us)

            %ylabel('STC')
            ylabel(sprintf('STC filter %d',ii))
            xlabel('time to spike (s)')
            box off
            
            % plot histogram of score!
            subplot(r,c,i*2+2)
            hist(score(:,ii))

            xlabel('score')
            %xlabel(sprintf('score onto the STC filter %d',ii))
            ylabel('count')
            title(sprintf('histogram of scores (var=%.3f, eig val=%.3f)', var(score(:,ii)),ev(ii)))
            box off
            
            %%
            figure(2)   % plot spatial pattern in a separate figure
            clf
%             if exist('xy', 'var')
%                 plot_stim_slices_with_mask(u(:,ii), sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height, FLIP_XY, xy)
%             else
%                 plot_stim_slices_with_mask(u(:,ii), sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height, FLIP_XY)                
%             end
%             hold on;
            %plot_MEA_param(ab,cd)
            %plot(xy(2), xy(1), 'ro')
            
            % fit ellipse fot STC filters
            [STC_pos_RF, STC_neg_RF, STC_strongest_RF] = calc_RF_from_STC_slice(u(:,ii), sta_num_samples, max(XX(mask>0))-min(XX(mask>0))+1, max(YY(mask>0))-min(YY(mask>0))+1, sampling_rate, FLIP_XY);
            
            if ~isempty(STC_strongest_RF)
                disp(STC_strongest_RF)
                %plot_RF(strongest_RF, ~FLIP_XY)
            end            
        
            
            set(gcf, 'paperposition', [0 0 24 9])
            set(gcf, 'papersize', [24 9])

            saveas(gcf, sprintf('STC_in_ROI_%s_%s_eig%d.png',mask_name, channel_names{n},ii))
            saveas(gcf, sprintf('STC_in_ROI_%s_%s_eig%d.pdf',mask_name, channel_names{n},ii))          
            
        end
        
        %%
        figure(3) % plot STA for comparison
        if exist('xy', 'var')
            plot_stim_slices_with_mask(sta_ROI{n}, sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height, FLIP_XY, xy, RFs{n})
        else
            plot_stim_slices_with_mask(sta_ROI{n}, sta_num_samples, XX(mask(:)>0), YY(mask(:)>0), width, height, FLIP_XY, [], RFs{n})
        end
        hold on

%         % plot re-fit ellipse
%         if ~isempty(strongest_RF)
%             subplot(2,5, strongest_RF.slice)
%             plot_RF(strongest_RF, ~FLIP_XY)
%         end


        set(gcf, 'paperposition', [0 0 24 9])
        set(gcf, 'papersize', [24 9])

        saveas(gcf, sprintf('STA_in_ROI_%s_%s.png', mask_name, channel_names{n}))
        saveas(gcf, sprintf('STA_in_ROI_%s_%s.pdf', mask_name, channel_names{n}))
            

        figure(1) % go back to figure 1
    else % no significant eigen values 
         % plot us in gray for checking
         
         % plot the smallest eig vec in gray
        figure(1)
        subplot(r,c,3)
        us = reshape(u(:,length(ev)),sta_num_samples,[]);
%        us = reshape(u(:,1),sta_num_samples,[]);
        
        
        plot(gridT, us, 'color', 0.6*[1 1 1])
        
%         title('eig. vec. for the largest eig. val.')        
        title('eig. vec. for the smallest eig. val.')        
        box off

        
         % plot the 2nd smallest eig vec in gray
        subplot(r,c,4)
%         us = reshape(u(:,length(ev)),sta_num_samples,[]);
        us = reshape(u(:,length(ev)-1),sta_num_samples,[]);
        plot(gridT, us, 'color', 0.6*[1 1 1])

        title('eig. vec. for the 2nd smallest eig. val.')        
        
        box off
        
        drawnow
    end

    % save ONLY WHEN significant eigen value is found 
    if ~isempty(idx_significant_ev)
        %%
        figure(1)
        
        set(gcf, 'paperposition', [0 0 24 10*r])
        set(gcf, 'papersize', [24 10*r])

        saveas(gcf, sprintf('STC_in_ROI_%s_%s.png', mask_name, channel_names{n}))
        saveas(gcf, sprintf('STC_in_ROI_%s_%s.pdf', mask_name, channel_names{n}))
    end
  
end


%% Save STC Results
save STC.mat stc_u stc_ev stc_significant_ev_idx

%% plot evs in one plot
clf; 
subplot(121)
hold on
for n=1:length(stc_ev)
    plot(stc_ev{n})
end
legend (get_channel_names(channel_names, num_significant_evs>0))

subplot(121)



%% plot numbers of significant eivenvalues
figure
bar(num_significant_evs);   box off

title(sprintf('Sigificant eigen values found in %d cells.', sum(num_significant_evs>0)))

set(gca,'ytick',[0 1 2 3])
ylabel('number of significant eigen values');

set(gca,'xtick', channel_index_to_analyze)
% xlabel('channel index');
set(gca, 'xticklabel', get_channel_names(channel_names, channel_index_to_analyze))
xlabel('channel names')

set(gcf, 'paperposition', [0 0 12 7])
set(gcf, 'papersize', [12 7])

saveas(gcf, 'num_significatn_eigen_values.pdf')
saveas(gcf, 'num_significatn_eigen_values.png')


%% 
channels_with_significant_eigen_values = channel_names(num_significant_evs>0)

save channels_with_significant_eigen_values channels_with_significant_eigen_values


%% plot STA  and STC results on the mosaic

clf; 

for n=channel_index_to_analyze
    channel_names{n}
    stc_significant_idx = stc_significant_ev_idx{n}
    
    if ~isfield(RFs{n}, 'type')
        warning(sprintf('skip channel %s', channel_names{n}))
        continue
    end
        
    % plot RF from STA
    if  strcmp(RFs{n}.type,'ON') % ON cell
        subplot(221); hold on 
        plot_RF(RFs{n}, FLIP_XY)
        STC_SUBPLOT=223;
    else %strcmp(RFs{n}.type,'OFF') % OFF cell        
        subplot(222); hold on
        plot_RF(RFs{n}, FLIP_XY)  % drawing codes are moved to this function
        STC_SUBPLOT=224;
    end
    
    % plot STC eig signs
    subplot(STC_SUBPLOT); hold on

    if isfield(RFs{n}, 'mean')
        center = RFs{n}.mean;
    else
        center = RFs{n}.peak_center(2:3);
    end
        
    if ~isempty(stc_significant_idx)
        for idx = stc_significant_idx
            if idx < (size(u,1)/2)   % larger
                tt=text(center(2), center(1), channel_names{n}(4:end), 'HorizontalAlignment','left');
                tt.Color = [1 0 0];                
                disp('larger ev')
            else   % smaller
                tt=text(center(2), center(1), channel_names{n}(4:end), 'HorizontalAlignment','right');
                tt.Color = [0 0 1];
                disp('smaller ev')
            end
        end
    else
        tt=text(center(2), center(1), channel_names{n}(4:end), 'HorizontalAlignment','center');
        tt.Color = [0 0 0];
        disp('no significant ev')
    end
end

% STA info
subplot(221)
xlabel('x')
ylabel('y')
title('receptive field from STA (ON cell)')

plot_MEA_param(ab, cd);
axis xy
%axis ([1 width  1 height])

subplot(222)
xlabel('x')
ylabel('y')
title('receptive field from STA (OFF cell)')

plot_MEA_param(ab, cd);
axis xy
%axis ([1 width  1 height])


% STC info
subplot(223) 
axis xy
plot_MEA_param(ab, cd);
%axis ([1 width  1 height])


title('STC with larger(red)/smaller(blue) eig. val. (ON cell)')

subplot(224) 
axis xy
plot_MEA_param(ab, cd);
%axis ([1 width  1 height])

title('STC with larger(red)/smaller(blue) eig. val. (OFF cell)')


set(gcf, 'paperposition', [0 0 10 8])
set(gcf, 'papersize', [10 8])

saveas(gcf, sprintf('mosaic_STA_STC.png'))
saveas(gcf, sprintf('mosaic_STA_STC.pdf'))
    
return 


