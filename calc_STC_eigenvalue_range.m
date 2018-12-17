function  [ev_range, evs, stc, sta] = calc_STC_eigenvalue_range(stim, spikeTrain, num_samples_per_window, num_random_shift, random_shift_range)

    %% set default params for shuffle
    if nargin<4
        num_random_shift = 1;
    end
    
    if nargin<5
        random_shift_range = num_samples_per_window * [1 10];
    end
    
    
    ev_min=[];
    ev_max=[];
    
    for n=1:num_random_shift
        %% random shift
        random_shift = round(diff(random_shift_range) * rand(1)  + random_shift_range(1));
    
        
        %% calc STA and STC
        [sta, stc] = calc_STA_and_STC(stim(random_shift:end,:), spikeTrain(1:end-random_shift+1,:), num_samples_per_window);



        %% analyze STC result
        if iscell(stc)
            for i=1:size(stc,1)
                [~, D, ~] = svd(stc{i,i});
                ev{i} = diag(D);
                
                % store max and min 
                ev_max = [ev_max max(ev{i})];
                ev_min = [ev_min min(ev{i})];
            end
        else
            [~, D, ~] = svd(stc);
            ev = diag(D);
            
            % store max and min 
            ev_max = [ev_max max(ev{i})];
            ev_min = [ev_min min(ev{i})];
                
        end
        
        % store evs
        evs{n} = ev;
                        
    end
    
    %% for debug
    if 1==0
        %%
        close all
        subplot(121)
        hist(ev_min)
        subplot(122)
        hist(ev_max)
    end

    
    % calc range!
    ev_range=[mean(ev_min) mean(ev_max)];
    
return



%% HOW TO USE
ev_range = calc_STC_eigenvalue_range(stim, spikeTrain(:,n), STA_num_samples, 50, [1 10]);


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
    