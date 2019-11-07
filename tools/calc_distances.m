function [scores, score_dist, ortho_dist] =  calc_distances(X_centered, U, L)
% inputs
% X_centerd - rows of data samples with mean removed 
% U - columns of eivenvectors
% L - eigenvalues

%% calc score
scores = X_centered*U;

%% calc Mahalanobis distance along PCs
% make L a row vector
if size(L,1) ~= 1
    L = L';
end
score_dist = sqrt(sum(bsxfun(@times, scores.^2, 1./L),2));

%% calc orthogonal distance
X_projected = scores * U';  % the Perpendicular Foot
X_diff = X_centered - X_projected;
ortho_dist = sqrt(sum(X_diff.^2,2));

return



