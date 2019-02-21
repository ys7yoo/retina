function X = project_out_components(X, row_vectors)

%% more general implementation

if nargin<2
    row_vectors = mean(X);
end

num_vectors = size(row_vectors,1);

for i=1:num_vectors
    X = project_out_sta(X, row_vectors(i,:));
end


return

%% initial implementation
% normalize vectors
L2normSqr = sqrt(sum(row_vectors.^2,2));
row_vectors = diag(1./L2normSqr)*row_vectors;

for i=1:num_vectors

    coef = X*row_vectors(i,:)';

    X = X - coef*row_vectors(i,:); 
end

return 