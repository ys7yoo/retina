function  [evs, num_spikes] = bootstrap_STC_eigen_value(stim, spikeTrain, num_samples_per_window, num_random_shift, random_shift_range, sta_to_project_out)

    %% set default params for shuffle
    if nargin<4
        num_random_shift = 1;
    end
    
    if nargin<5
        random_shift_range = num_samples_per_window * [1 10];
    end
    
    if nargin<6
        sta_to_project_out = [];
    end
    
    % to store eigen values for each repeat
    dim = size(stim,2) * num_samples_per_window;
    evs = zeros(dim, num_random_shift);           
       
    for n=1:num_random_shift
        %% random shift & collect spike-triggered stim       
        [X, spikes, num_total_spikes] = shuffle_stim(stim, spikeTrain, num_samples_per_window, random_shift_range);
        

        %% project if requested 
        if ~isempty(sta_to_project_out)
            X = project_out_components(X, sta_to_project_out);
        end
        
        ev = calc_STC(X, spikes);
    
        
        %[~, ev] = calc_STA_and_STC(stim(random_shift+1:end-(shift_max-random_shift),:), spikeTrain(1:end-shift_max,:), num_samples_per_window, sta_to_project_out);
                
        % store evs
        evs(1:length(ev),n) = ev;    

        % store evs
        if nargout>1
            num_spikes(:,n) = num_total_spikes;
        end
                        
    end


%     % calc confidence interval using non-zero eigen values
%     evs(evs<1e-6) = nan;
% 
% 
%     mm = nanmean(evs,2);
%     ss = nanstd(evs,[],2);
%     
%     Confidence interval calculation is upto the user!
%     % 95% confidence interval for each eig
%     ev_upper = mm+1.96*ss;
%     ev_lower = mm-1.96*ss;
%     
%     % 99% confidence interval for each eig
%     ev_upper = mm+2.576*ss;
%     ev_lower = mm-2.576*ss;
    

%     % 99% significance interval
%     ev_range = [mm(end)-2.576*ss(end) mm(1)+2.576*ss(1)]; 
    
% %     % Option 1. calc range by min and max 
% %     ev_range=[mean(ev_min) mean(ev_max)];
%     
%     % Option 2. calc range by mean and var
%     
%     ev_mean_global = mean(ev_mean);
%     ev_var_global = mean(ev_var);
%     ev_std = sqrt(ev_var_global);
% 
% %     % 95% significance interval
% %     ev_range = ev_mean_global + ev_std*1.96*[-1 1];
%     
%     % 99% significance interval
%     ev_range = ev_mean_global + ev_std*2.58*[-1 1];
    
    
    
    
return



%% HOW TO USE
ev_range = calc_STC_eigenvalue_range(stim, spikeTrain(:,n), STA_num_samples, 50, [1 10]);


%% Further analysis (with mask)
[ev_range, evs] = calc_STC_eigenvalue_range(stim(:,mask(:)), spike_train(:,n), sta_num_samples, 50, sta_num_samples*[10 50]);


%% test with different pojections 

% save sample_X X spikes
% load sample_X X spikes
ev0 = calc_STC(X, spikes);

ev1 = calc_STC(project_out_components(X, mean(X)), spikes); 
% the last eig goes to zero!
% the slope becomes slightly steeper



LINE_STYLES = {'-.','-'};
MARKERS = {'o', '+'};
COLORS = {'r', 'b'};
clf

subplot(221)
h=plot([ev0 ev1]);
set(h(1), 'LineStyle', LINE_STYLES{1})
%set(h(1), 'Marker', MARKERS{1})
set(h(1), 'Color', COLORS{1})
set(h(2), 'LineStyle', LINE_STYLES{2})
%set(h(2), 'Marker', MARKERS{2})
set(h(2), 'Color', COLORS{2})
box off

ylabel('eigen value')
xlabel('index')

subplot(222)
h=qqplot(ev0, ev1)
xlabel('percentile of ev_0')
ylabel('percentile of ev_1')


% bin=0:10:(max([ev0 ev1])+1);
% cnt = histc([ev0 ev1], bin);
% bar(bin, cnt)
% %plot(bin, cnt)
% box off


subplot(223)
h=plot([ev0 ev1]);
set(h(1), 'LineStyle', LINE_STYLES{1})
set(h(1), 'Marker', MARKERS{1})
set(h(1), 'Color', COLORS{1})
set(h(2), 'LineStyle', LINE_STYLES{2})
set(h(2), 'Marker', MARKERS{2})
set(h(2), 'Color', COLORS{2})
box off

ylabel('eigen value')
xlabel('index')

set(gca, 'xlim', [0 10])


subplot(224)
h=plot([ev0 ev1]);
set(h(1), 'LineStyle', LINE_STYLES{1})
set(h(1), 'Marker', MARKERS{1})
set(h(1), 'Color', COLORS{1})
set(h(2), 'LineStyle', LINE_STYLES{2})
set(h(2), 'Marker', MARKERS{2})
set(h(2), 'Color', COLORS{2})
box off

ylabel('eigen value')
xlabel('index')

set(gca, 'xlim', length(ev0) + [-10 0])


set(gcf, 'paperposition', [0 0 24 20])
set(gcf, 'papersize', [24 20])

saveas(gcf, sprintf('projection_and_eigen_values.png'))
saveas(gcf, sprintf('projection_and_eigen_values.pdf'))
            

%%
mean(evs(:))

clf
hist(evs(:),100)
box off

evs_mean = mean(evs(:))
hold on
YLIM = get(gca,'ylim');
plot(evs_mean*[1 1], YLIM, 'r--');

%% new plot
close all

CONFIDENCE=2.576; %% 99% 
%CONFIDENCE=1.96;  % 95% 

subplot(121)
% plot(ev, 'linewidth', 2); hold on;
plot(evs, 'color', 0.5*[1 1 1]); hold on 
plot(mm+CONFIDENCE*ss, 'r--')
plot(mm-CONFIDENCE*ss, 'r--')
set(gca,'yscale', 'log')

subplot(122)
% plot(ev, 'linewidth', 2); hold on;
plot(evs, 'color', 0.5*[1 1 1]); hold on 
plot(mm+CONFIDENCE*ss, 'r--')
plot(mm-CONFIDENCE*ss, 'r--')
set(gca,'yscale', 'log')

set(gca,'xlim', length(ev)+ [-5 0])

% subplot(223)
% hist(evs(1,:))
% 
% subplot(224)
% hist(evs(end,:))



        
%%
subplot(221)
%hist(mean(evs,2))
hist(evs(:),100)

hold on
YLIM = get(gca,'ylim');
plot(evs_mean*[1 1], YLIM, 'r--');
plot(ev_range(1)*[1 1], YLIM, 'r--');
plot(ev_range(2)*[1 1], YLIM, 'r--');
plot(mean(min(evs))*[1 1], YLIM, 'g--');
plot(mean(max(evs))*[1 1], YLIM, 'g--');


subplot(222)
hist(std(evs,[],2))
title('histogram of std')

subplot(223)
hist(min(evs))
mean(min(evs))
title('histogram of min')

subplot(224)
hist(max(evs))
mean(max(evs))
title('histogram of max')
%% Further analysis
[ev_range, evs, stc, sta] = calc_STC_eigenvalue_range(stim, spikeTrain(:,n), STA_num_samples, 50, [1 10]);


close all

for i=1:10
    subplot(5,2,i)
    plot(evs{i}{61}, 'o')
    
    hold on
    XLIM=get(gca,'xlim');
    plot(XLIM, ev_range(1)*[1 1], 'r--')
    plot(XLIM, ev_range(2)*[1 1], 'r--')
end


%% let's compare with no shift
[~, evs, stc, sta] = calc_STC_eigenvalue_range(stim, spikeTrain(:,n), STA_num_samples, 1, [1 1]);


close all
r=1;c=2;

idx_pixel=61;

subplot(r,c,1)
imagesc(stc{idx_pixel,idx_pixel})
%colormap gray
box off
title('STC for the pixel with the largest variance')

subplot(r,c,2)
plot(evs{1}{idx_pixel}, 'ok'); hold on
% plot(1, ev(1), 'b*')
%plot(length(ev), ev(end), 'r*')   
hold on
XLIM=get(gca,'xlim');
plot(XLIM, ev_range(1)*[1 1], 'r--')
plot(XLIM, ev_range(2)*[1 1], 'r--')


box off
title ('eigen values')


box off
    



%% example plot (histogram of eiven values)

close all
hist(reshape(cell2mat(evs{1}),[],1)); box off;

% set(gcf, 'paperposition', [0 0 6 3.5])
% set(gcf, 'papersize', [6 3.5])

set(gcf, 'paperposition', [0 0 8 6])
set(gcf, 'papersize', [8 6])


saveas(gcf, 'eigenvalue_histogram.png')
saveas(gcf, 'eigenvalue_histogram.pdf')
    