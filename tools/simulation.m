clear

N1 = 1000
N2 = 1000
dim = 2

AXIS = [-3.5 3.5 -3.5 3.5]


%% generate data
X = randn(N1,dim)

clf
subplot(221)
plot(X(:,1), X(:,2), '+b')
box off


%% calc covariance matrix
Cov = X'*X / N1

hold on; plot_ellipse([0 0], Cov, 'b')

axis(AXIS)

%% introduce filter
%filt = [1 0.5; 0 1]


% suppress along 
% R = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)]
% filt = R'*diag([1 0.5])*R

filt = diag([1 0.5])



X_filtered = X*filt'


Cov_filtered = X_filtered'*X_filtered / N1

subplot(222)
plot(X_filtered(:,1), X_filtered(:,2), 'or')
box off

hold on; plot_ellipse([0 0], Cov_filtered, 'r')


axis(AXIS)

%% analysis
eig0 = svd(Cov)

eig1 = svd(Cov_filtered)

eig0 ./ eig1


%% plot mixed
%subplot(223)
figure(2)
clf; hold on
%plot(X_mix(:,1), X_mix(:,2), 'og')
plot(X(:,1), X(:,2), '+b')
plot(X_filtered(:,1), X_filtered(:,2), 'or')
box off

plot_ellipse([0 0], Cov, 'b--')
plot_ellipse([0 0], Cov_filtered, 'r--')

axis(AXIS)
xlabel('X_1')
ylabel('X_2')


set(gcf, 'paperposition', [0 0 12 10])
set(gcf, 'papersize', [12 10])
saveas(gcf, 'simulation.pdf')
saveas(gcf, 'simulation.png')

% %% Mix the dataset 
% 
% lamb = 0.5;
% X_mix = [X; X_filtered];
% Y_mix = [-ones(N1,1); ones(N2,1)];
% %X_mix = (1-lamb)*X + lamb*X_filtered;
% 
% Cov_mix = X_mix'*X_mix / size(X_mix,1)
% 
% eig_2 = svd(Cov_mix)
% 
% 
% subplot(223)
% plot(X_mix(:,1), X_mix(:,2), 'og')
% box off
% 
% hold on; plot_ellipse([0 0], Cov_mix, 'g')
% 
% axis(AXIS)






