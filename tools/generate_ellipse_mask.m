function mask =  generate_ellipse_mask(mu, C, h, w)

% generate 95% confidence interval ellipse mask

[XX, YY] = meshgrid(1:w, 1:h);

%R = chol(C);

tt = [XX(:)-mu(1) YY(:)-mu(2)]*inv(chol(C));
phi = (tt(:,1)).^2 + (tt(:,2)).^2;

phi = reshape(phi, h, w);

%[U, D, V] = eig(C);
%phi = (XX-mu(1)).^2 / D(1,1).^2 + (YY-mu(2)).^2 / D(2,2).^2;


% Get the 95% confidence interval error ellipse
%chisquare_val = 2.4477;
mask = 2.4477^2 - phi;
%plot_ellipse(, , 'r-');


return 

%% hot to run this code
%mu = [0 0];
%C = [1 0.5; 0.5 2];
mu = [10 10];
C = [5 3; 3 10];

close all
mask = generate_ellipse_mask(mu, C, 26, 26);
%imshow(mask)
%contour(mask)
imshow(mask>0);
%surf(mask)
hold on;
plot_ellipse(mu, C, 'r--')
axis xy


