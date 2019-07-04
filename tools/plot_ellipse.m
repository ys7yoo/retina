function plot_ellipse(avg, covariance, LINE_TYPE, FLIP_XY)

if nargin<3
    LINE_TYPE = 'r--';
end
if nargin<4
    FLIP_XY = false;
end
    

%% simpler version based on linear algebra

%% new function for plot 
theta = linspace(0,2*pi);

circ = [cos(theta') sin(theta')];

if max(covariance(:)) > 0
    % disp(covariance)
    ellipse = avg + circ*chol(covariance+diag([eps, eps]))*diag([2.4477, 2.4477]);


    if ~FLIP_XY
        %plot(circ(:,1), circ(:,2), 'k--'); hold on
        plot(ellipse(:,1), ellipse(:,2), LINE_TYPE); 
    else
        plot(ellipse(:,2), ellipse(:,1), LINE_TYPE); 
    end
else % single pixel
    if ~FLIP_XY
        plot(avg(:,1), avg(:,2), ['o' LINE_TYPE])
    else
        plot(avg(:,2), avg(:,1), ['o' LINE_TYPE])
    end
        

end

return

%% modified code from http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

[eigenvec, eigenval ] = eig(covariance);

eigenval = diag(eigenval);


if sum(eigenval) == 0
    return
end

% Get the index of the largest eigenvector
[~, largest_eigenvec_ind_c] = max(eigenval); 
%[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = eigenval(2);
%     smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = eigenval(1);
%     smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end



% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + avg(1), r_ellipse(:,2) + avg(2), LINE_TYPE)
hold on;



return 

%% test for enlargement

close all
C = [1 0; 0 1];
eig(C)
plot_ellipse([0 0], C)

hold on 
plot_ellipse([0 0], C*2^2)

%%
[U, D, V] = svd(C)




%%  test


m = [7.9909    5.2512]
 
C = [0.0090,    0.0023
     0.0023,    0.4392]


% m = [5, 5]

% C = [ 1, 0.5
%       0.5, 2]
%% generate data
X=mvnrnd(m, C, 100)

clf
plot(X(:,1), X(:,2), 'o')

axis ([1 10 1 10])

hold on
plot_ellipse(m, C, 'r-')
%%
%     7.9040    5.3098
% 
%     0.0868   -0.0074
%    -0.0074    0.5825
   