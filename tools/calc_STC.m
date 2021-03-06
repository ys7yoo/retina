function [stc_eig_val, stc_eig_vec, S] = calc_STC(X, spikes)


num_total_spikes = sum(spikes);
%% calc STC
switch nargout 
    case 1   % STA and STC eigen value only (COVARIANCE NOT NEEDED)
        Xs = bsxfun(@times, X, sqrt(spikes));

        %[~, D, ~] = svd(Xs);
        %stc_eig_val = diag(D).^2;
        s = svd(Xs);

        stc_eig_val = s.^2;
        

    case 2  % STA and STC eigen value & eiven vectors (COVARIANCE NOT NEEDED)

        % subtract mean & scale down 
        Xs = bsxfun(@times, X, sqrt(spikes));

        [~, D, stc_eig_vec] = svd(Xs);

        stc_eig_val = diag(D).^2;

        %stc_eig_vec = flip_column_sign(stc_eig_vec, sta-0.5);  % flip according to sta (for better visualization)

    case 3 % full algorithm with covariance
        S = X'*bsxfun(@times, X, spikes);  % DO NOT DIVIDE by num_total_spikes FOR BETTER NUMERICAL PRECISION. So, actually S is a scatter matrix.

        [stc_eig_vec, d, ~] = svd(S);

        stc_eig_val = diag(d);
        
        % devide S here
        S = S / (num_total_spikes-1);

        %stc_eig_vec = flip_column_sign(stc_eig_vec, sta-0.5);  % flip according to sta (for better visualization)
end

%% normalize eig values 
stc_eig_val = stc_eig_val / (num_total_spikes-1);

