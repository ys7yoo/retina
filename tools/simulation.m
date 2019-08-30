clear

N1 = 1000
N2 = 1000
dim = 2

%AXIS = [-3.5 3.5 -3.5 3.5]
AXIS = [-4 4 -4 4];


%% generate data
X = randn(N1,dim);

figure(1)
clf
subplot(221)
plot(X(:,1), X(:,2), '+b')
box off

xlabel('X_1')
ylabel('X_2')
title('samples for no spike')

%% calc covariance matrix
Cov = X'*X / N1

hold on; plot_ellipse([0 0], Cov, 'b')

axis(AXIS)

%% introduce filter
%filt = [1 0.5; 0 1]


% suppress along 
R = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)]
filt = R'*diag([1 0.5])*R

% filt = diag([1 0.5])
mu_filtered = [1.5, 1]


X_filtered = X*filt';


Cov_filtered = X_filtered'*X_filtered / N1

subplot(222)
plot(mu_filtered(1)+X_filtered(:,1), mu_filtered(2)+X_filtered(:,2), 'or')
box off

xlabel('X_1')
ylabel('X_2')

hold on; plot_ellipse(mu_filtered, Cov_filtered, 'r')

title('samples for spike')

axis(AXIS)

%% eigen analysis
eig0 = svd(Cov)
[U, D, V] = svd(Cov);


eig1 = svd(Cov_filtered)
[U1, D1, V1] = svd(Cov_filtered);

plot_principal_axes(Cov_filtered, mu_filtered)


eig0 ./ eig1


%% anlayze projected components

% project onto eiven vectors (U1)
projected=X_filtered*U1;

vars = var(projected);

subplot(2,2,3)
[count, bin] = hist(projected(:,1));
h1 = bar(bin, count,'facecolor', 'r');
xlabel('PC_1')
ylabel('count of projected values')
box off
title(sprintf('histogram of projected values for PC_1 (var=%.2f)',vars(1)))



subplot(2,2,4)
[count, bin] = hist(projected(:,2));
h2 = bar(bin, count,'facecolor', 'r');
xlabel('PC_2')
ylabel('count of projected values')
box off
title(sprintf('histogram of projected values for PC_2 (var=%.2f)',vars(2)))



set(gcf, 'paperposition', [0 0 20 18])
set(gcf, 'papersize', [20 18])
saveas(gcf, 'simulation1.pdf')
saveas(gcf, 'simulation1.png')


%% plot mixed samples in one image
%subplot(223)
figure(2)
clf; 
subplot(121)
hold on
%plot(X_mix(:,1), X_mix(:,2), 'og')
plot(X(:,1), X(:,2), '+b')
%plot(X_filtered(:,1), X_filtered(:,2), 'or')
plot(mu_filtered(1)+X_filtered(:,1), mu_filtered(2)+X_filtered(:,2), 'or')
box off


plot([0 mu_filtered(1)], [0 mu_filtered(2)], 'k-', 'LineWidth', 3)
text(mu_filtered(1)/2-0.1, mu_filtered(2)/2+0.2, 'STA', 'FontSize',20,  'HorizontalAlignment','right')

plot_ellipse([0 0], Cov, 'b--')
plot_ellipse(mu_filtered, Cov_filtered, 'r--')

axis(AXIS)
axis equal
xlabel('PC_1')
ylabel('PC_2')



subplot(122)
hist(projected)
title('histogram of projected values')
legend (sprintf('PC_1 (var=%.2f)', vars(1)), sprintf('PC_2 (var=%.2f)', vars(2)))
box off


set(gcf, 'paperposition', [0 0 32 10])
set(gcf, 'papersize', [32 10])
saveas(gcf, 'simulation2.pdf')
saveas(gcf, 'simulation2.png')

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






