function [idx_found] = find_significant_eigen_values(ev, u, stim, spike_train, num_samples_per_window, num_random_shift, random_shift_range, sta_to_project_out, idx_candidate)

CONFIDENCE=2.576; %% 99% 
%CONFIDENCE=1.96;  % 95% 

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
idx_candidate

% [idx_candidate_large idx_candidate_small] % for debug

idx_found = [];
component_to_project_out = sta_to_project_out;

flag_continue = true;

while flag_continue && idx_candidate(1)<=MAX_NUM_COMPONENTS && idx_candidate(end) >= (num_eigen_values-MAX_NUM_COMPONENTS)

    [mm, ss, evs, num_spikes] = calc_STC_eigenvalue_range(stim, spike_train, num_samples_per_window, num_random_shift, random_shift_range, component_to_project_out);

    % identify most significant eigen value between the largest or smallest
    % eigven values
    abs_z_score = abs((ev - mm(1:length(ev))) ./ ss(1:length(ev)));

    if abs_z_score(1) > abs_z_score(end)
        if abs_z_score(1) > CONFIDENCE % above chance level?
            disp(sprintf('%dth eigenvalue is significant',idx_candidate(1)))                
            idx_found = [idx_found idx_candidate(1)];

            % for next 
            component_to_project_out = [component_to_project_out; u(:,1)'];

            idx_candidate = idx_candidate+[1 0]
            ev = ev(2:end);
            u = u(:,2:end);
        else
            % not found 
            flag_continue = false;
        end
    else
        if abs_z_score(end) > CONFIDENCE % above chance level?
            disp(sprintf('%dth eigenvalue is significant',idx_candidate(2)))
            idx_found = [idx_found idx_candidate(2)];

            % for next 
            component_to_project_out = [component_to_project_out; u(:,end)'];

            idx_candidate = idx_candidate+[0 -1]
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

%CONFIDENCE=2.576; %% 99% 
%CONFIDENCE=1.96;  % 95% 
% 99%
ev_upper = mm + 2.576*ss;
ev_lower = mm - 2.576*ss;


% 95%
% ev_upper = mm + 1.96*ss;
% ev_lower = mm - 1.96*ss;

plot(ev,'ko-', 'linewidth', 2)
hold on
plot(ev_upper, 'r--')
plot(ev_lower, 'r--')
plot(evs, 'color', 0.5*[1 1 1])
set(gca,'yscale','log')

set(gca,'xlim', length(ev)+[-4 0])

%     XLIM=get(gca,'xlim');
%     plot(XLIM, ev_range(1)*[1 1], 'r--')
%     plot(XLIM, ev_range(2)*[1 1], 'r--')

% plot(idx_significant_ev, ev(idx_significant_ev), '*k')
ylabel('eigen values of STC')
box off