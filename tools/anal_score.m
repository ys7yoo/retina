function anal_score(channel_name)
%% for detail analysis for each channel
% data saved by the following command: save ch33b.mat X spikes num_total_spikes

%load ch_42b.mat
% channel_name = 'ch_42b'
load(channel_name)

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
xlabel('index')
ylabel('eigen value')

subplot(222)
hist(stc_eig_val)
%hist(sqrt(stc_eig_val))
box off
xlabel('eigen value')
ylabel('count')
title('histogram of eigen values')

% calc scores for each dim
r = 2
score = X*stc_eig_vec;

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
cov_small = cov(score(:,end), score(:,end-1))
plot_ellipse([0 0], cov_small)

xlabel(sprintf('ev%d',length(stc_eig_val)))
ylabel(sprintf('ev%d',length(stc_eig_val)-1))
axis equal
axis ([-2.5 2.5 -2 2])


% saveas(gcf, sprintf('%s_score.pdf', channel_name))
saveas(gcf, sprintf('%s_score.png', channel_name))


return


%% call for all the channels
for n = 1:length(channel_names)
    anal_score(channel_names{n})
end