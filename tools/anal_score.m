function anal_score(channel_name, hist_num_bins)
%% for detail analysis for each channel
% data saved by the following command: save ch33b.mat X spikes num_total_spikes

%load ch_42b.mat
% channel_name = 'ch_42b'
load(channel_name)

% calc STA
sta = spikes'*X/num_total_spikes;

% calc STFC
[stc_eig_val, stc_eig_vec] = calc_STC(X, spikes);
                        
stc_eig_vec = flip_column_sign(stc_eig_vec, sta);  % flip according to sta (for better visualization)

%  select only non-zero eigen values
r = length(find(stc_eig_val>1e-15));
stc_eig_val = stc_eig_val(1:r);
stc_eig_vec = stc_eig_vec(:,1:r);

% plot eig val
clf
subplot(221)
plot(stc_eig_val, 'o--', 'markersize', 5)
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
score = X*stc_eig_vec;


% plot in 2D plots
subplot(223)
sc1 = scatter(score(:,1), score(:,2), '.'); %, '.', 'alpha', 0.2)
sc1.MarkerEdgeAlpha=0.4;

hold on
cov12 = score(:,1:2)'*score(:,1:2)/(size(score,1)-1);
plot_ellipse([0 0], cov12)

xlabel('ev1')
ylabel('ev2')
axis equal

axis ([-2.5 2.5 -2 2])

subplot(224)
sc1 = scatter(score(:,r), score(:,r-1), '.'); %, '.', 'alpha', 0.2)
sc1.MarkerEdgeAlpha=0.4;

hold on
cov_small = score(:,r:-1:r-1)'*score(:,r:-1:r-1)/(size(score,1)-1);
plot_ellipse([0 0], cov_small)

xlabel(sprintf('ev%d',r))
ylabel(sprintf('ev%d',r-1))
axis equal
axis ([-2.5 2.5 -2 2])


set(gcf, 'paperposition', [0 0 8 6])
set(gcf, 'papersize', [8 6])

saveas(gcf, sprintf('%s_score.pdf', channel_name))
saveas(gcf, sprintf('%s_score.png', channel_name))

%% 
%% plot histogram of scores in separate figures

% before plotting, determine figure params: figure size, common limits for x-axis
FIGURE_W = 6;
FIGURE_H = 4.5;

if nargin<2
    hist_num_bins = 30;      % default value
end

X_MAX = ceil(max(max(abs(score(:,[1 2 r-1 r]))))*2)/2;
hist_bins = linspace(-X_MAX, X_MAX, hist_num_bins);
XLIM = X_MAX*[-1 1];


figure
idx=1
plot_score_histogram(score, idx, hist_bins, XLIM);
set(gcf, 'paperposition', [0 0 FIGURE_W FIGURE_H])
set(gcf, 'papersize', [FIGURE_W FIGURE_H])
saveas(gcf, sprintf('%s_score_hist_%d.pdf', channel_name, idx))
saveas(gcf, sprintf('%s_score_hist_%d.png', channel_name, idx))

figure
idx = 2
plot_score_histogram(score, idx, hist_bins, XLIM);
set(gcf, 'paperposition', [0 0 FIGURE_W FIGURE_H])
set(gcf, 'papersize', [FIGURE_W FIGURE_H])
saveas(gcf, sprintf('%s_score_hist_%d.pdf', channel_name, idx))
saveas(gcf, sprintf('%s_score_hist_%d.png', channel_name, idx))

figure
idx = r-1
plot_score_histogram(score, idx, hist_bins, XLIM);
set(gcf, 'paperposition', [0 0 FIGURE_W FIGURE_H])
set(gcf, 'papersize', [FIGURE_W FIGURE_H])
saveas(gcf, sprintf('%s_score_hist_%d.pdf', channel_name, idx))
saveas(gcf, sprintf('%s_score_hist_%d.png', channel_name, idx))

figure
idx = r
plot_score_histogram(score, r, hist_bins, XLIM);
set(gcf, 'paperposition', [0 0 FIGURE_W FIGURE_H])
set(gcf, 'papersize', [FIGURE_W FIGURE_H])
saveas(gcf, sprintf('%s_score_hist_%d.pdf', channel_name, idx))
saveas(gcf, sprintf('%s_score_hist_%d.png', channel_name, idx))


return


function XLIM = plot_score_histogram(score, idx, bins, XLIM)

if nargin<3
    hist(score(:,idx))
else
    % use the provided bins
    hist(score(:,idx), bins)
end
title (sprintf('histogram of score %d',idx))
xlabel (sprintf('score %d',idx))
ylabel('count')
box off


% set XLIM, if needed
if nargin>3
    set(gca, 'xlim', XLIM);
% else
%     % set symmetric limits
%     XLIM = get(gca, 'xlim');
%     XLIM_MAX = max(XLIM);
%     set(gca, 'xlim', XLIM_MAX*[-1 1]);
end

% return XLIM
XLIM = get(gca, 'xlim');


return








%% call for all the channels
for n = 1:length(channel_names)
    anal_score(channel_names{n})
end