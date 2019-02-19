function X_residue = project_out_sta(X, sta)

% X - rows of stim
% sta - a row vector to be projected out 

if nargin<2
    sta = mean(X);
end



L2normSqr = sta*sta';
coef = X*sta' / L2normSqr;

X_residue = X - coef*sta; 

return 


%% test code

X_residue = project_out_sta(X,sta);
