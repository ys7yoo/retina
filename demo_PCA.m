num_samples = 5000;
mu = [4 6];
%mu = [0 0];
C = [1 0.5; 0.5 2];
R = chol(C);
Z = repmat(mu,num_samples,1) + randn(num_samples,2)*R

close all
plot(Z(:,1),Z(:,2),'.')
axis equal

%% 
CC=cov(Z)
%[U, D, V] = eig(CC);
[U, D, V] = svd(CC);

plot_directions(U)
D

%% 
coeff = pca(Z)

hold on
plot([0; coeff(1,1)], [0, coeff(2,1)], 'k--')
plot([0; coeff(1,2)], [0, coeff(2,2)], 'k--')
axis equal


%% 
plot_ellipse(mu, CC, 'r--')

%% 
mask = generate_ellipse_mask(mu, CC, 10, 10)
%figure(2)
imshow(mask>0)
%contour(mask)
axis xy


% hold on
% plot_ellipse(mu, CC, 'r--')


% hold on 
% plot(Z(:,1),Z(:,2),'.')
% axis equal
%% 