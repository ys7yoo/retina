function [idx_found] = find_significant_eigen_values(ev, u, stim, spike_train, num_samples_per_window, num_random_shift, random_shift_range, sta_to_project_out, idx_candidate)

% CONFIDENCE=3.090; % 99.9%  one-side
CONFIDENCE=2.326;  % 99%  one-side

%CONFIDENCE=2.576; %% 99%  two-sided
%CONFIDENCE=1.96;  % 95% 
%http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_HypothesisTest-Means-Proportions/BS704_HypothesisTest-Means-Proportions3.html

MAX_NUM_COMPONENTS=10;


% project STA component


ev = ev(ev>1e-5);
num_eigen_values = length(ev);
u = u(:,1:num_eigen_values);

if nargin < 8
    sta_to_project_out = [];
end
    
if nargin < 9
    idx_candidate = [1 num_eigen_values];
end
% idx_candidate  % log for debug

% [idx_candidate_large idx_candidate_small] % for debug

idx_found = [];
component_to_project_out = sta_to_project_out;

flag_continue = true;

while flag_continue && idx_candidate(1)<=MAX_NUM_COMPONENTS && idx_candidate(end) >= (num_eigen_values-MAX_NUM_COMPONENTS)

    %% bootstrapping by shuffling
    [evs, num_spikes] = bootstrap_STC_eigen_value(stim, spike_train, num_samples_per_window, num_random_shift, random_shift_range, component_to_project_out);
    
    %% decide significance interval
    mm = nanmean(evs,2);
    ss = nanstd(evs,[],2);

    %% identify most significant eigen value between the largest or smallest
    % eigven values
    z_score = (ev - mm(1:length(ev))) ./ ss(1:length(ev));

    if abs(z_score(1)) > abs(z_score(end)) % choose stronger one
        if z_score(1) > CONFIDENCE % above chance level?
            disp(sprintf('%dth eigenvalue is significant',idx_candidate(1)))                
            idx_found = [idx_found idx_candidate(1)];

            % for next 
            component_to_project_out = [component_to_project_out; u(:,1)'];

            idx_candidate = idx_candidate+[1 0];
            ev = ev(2:end);
            u = u(:,2:end);
        else
            % not found 
            flag_continue = false;
        end
    else
        if z_score(end) < -CONFIDENCE % above chance level?
            disp(sprintf('%dth eigenvalue is significant',idx_candidate(2)))
            idx_found = [idx_found idx_candidate(2)];

            % for next 
            component_to_project_out = [component_to_project_out; u(:,end)'];

            idx_candidate = idx_candidate+[0 -1];
            ev = ev(1:end-1);
            u = u(:,1:end-1);
        else
            % not found 
            flag_continue = false;
        end
    end

end

return





%% DEBUG

clf
subplot(211)
plot(ev, 'b-+', 'linewidth', 2); hold on; 
%plot(evs, '-', 'color', 0.5*[1 1 1])

CONFIDENCE=2.326;  % 99%  one-side

ev_upper = mm + CONFIDENCE*ss;
ev_lower = mm - CONFIDENCE*ss;

plot(ev_upper, 'r--')
plot(ev_lower, 'r--')

ylabel('eigen value')
xlabel('index')
box off

subplot(223)
plot(ev, 'b-+', 'linewidth', 2); hold on; 
plot(evs, '-', 'color', 0.5*[1 1 1])

CONFIDENCE=2.576; %% 99% 
%CONFIDENCE=1.96;  % 95% 
% 99%
ev_upper = mm + CONFIDENCE*ss;
ev_lower = mm - CONFIDENCE*ss;

plot(ev_upper, 'r--')
plot(ev_lower, 'r--')

ylabel('eigen value')
xlabel('index')
box off

set(gca,'xlim', [0 5])

subplot(224)
plot(ev, 'b-+', 'linewidth', 2); hold on; 
plot(evs, '-', 'color', 0.5*[1 1 1])

CONFIDENCE=2.576; %% 99% 
%CONFIDENCE=1.96;  % 95% 
% 99%
ev_upper = mm + CONFIDENCE*ss;
ev_lower = mm - CONFIDENCE*ss;

plot(ev_upper, 'r--')
plot(ev_lower, 'r--')

ylabel('eigen value')
xlabel('index')
box off

set(gca,'xlim', size(evs,1)+[-5 0])

%%
COMPONENT=3;
set(gcf, 'paperposition', [0 0 24 20])
set(gcf, 'papersize', [24 20])

saveas(gcf, sprintf('testing_eigen_values_%d.png',COMPONENT))
saveas(gcf, sprintf('testing_eigen_values_%d.pdf',COMPONENT))
            