function [m, C] = calc_mean_and_cov(X, prob)

% each row is a sample
[~, dim] = size(X);

m = prob'*X;

diff = bsxfun(@minus, X, m);


C = zeros(dim,dim);
for i = 1:dim
    for j = i:dim
        v = diff(:,i).*diff(:,j);
        C(i,j) = prob'*v;
        C(j,i) = C(i,j);
    end
end

    
return 

%% 
calc_mean_and_cov([1; 2; 3; 4], [0.25; 0.25; 0.25; 0.25])


%% 
X = [1 2; 3 5]
prob = [0.5; 0.5]
[m, C] =calc_mean_and_cov(X, prob)



    