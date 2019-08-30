function plot_principal_axes(Cov, mu)

if nargin<2
    mu = zeros(2,1);  % default: no offset
end

if nargin<3
    LINE_TYPE='k';
end

[U, D, V] = svd(Cov);

% plot principal axes 
plot([mu(1) mu(1)+U(1,1)], [mu(2) mu(2)+U(2,1)], LINE_TYPE, 'linewidth', 2); hold on
plot([mu(1) mu(1)+U(1,2)], [mu(2) mu(2)+U(2,2)], LINE_TYPE, 'linewidth', 2)
axis equal