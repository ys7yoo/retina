function [sta, stc_eig_val, stc_eig_vec, S] = calc_STA_and_STC(Stim, spike_train, n, sta_to_project_out)

% input:
%       Stim = (time) x (space)
%       spike_train = (time) x (spikes)
%       n = number of samples for analysis
%       sta_to_project_out

if nargin<4
    sta_to_project_out = [];
end

dim = size(Stim,2);     % dimension of the stimulus
spike_train(1:n-1,:) = 0;  % Ignore spikes before time n

% Compute spike-triggered STA and STC
spike_idx = find(spike_train>0);

spikes = spike_train(spike_idx);
% num_bins = length(spikes);
num_total_spikes = sum(spikes);


% construct design matrix
X = makeStimRows(Stim, n, spike_idx);


%% calc STA
% new implementation (when spikes are either 1 or 0, Matlab mean function is slightly faster!
%sta = mean(X);

sta = spikes'*X/num_total_spikes;

%% calc STC (FINAL ALGORITHM)

if nargout>1
    
    
    
    % subtract mean first
    X = bsxfun(@minus, X, sta);
    
    % project out sta, if requested
    if ~isempty(sta_to_project_out)
        X = project_out_sta(X, sta_to_project_out);
    end
    
    switch nargout 
        
        case 2   % STA and STC eigen value only (COVARIANCE NOT NEEDED)
            Xs = bsxfun(@times, X, sqrt(spikes));
            
            [~, D, ~] = svd(Xs);
            stc_eig_val = diag(D).^2;

        case 3  % STA and STC eigen value & eiven vectors (COVARIANCE NOT NEEDED)

            % subtract mean & scale down 
            Xs = bsxfun(@times, X, sqrt(spikes));

            [~, D, stc_eig_vec] = svd(Xs);
            stc_eig_val = diag(D).^2;
            
            
        case 4 % full algorithm with covariance
            S = X'*bsxfun(@times, X, spikes);  % DO NOT DIVIDE by num_total_spikes FOR BETTER NUMERICAL PRECISION. So, actually S is a scatter matrix.

            [stc_eig_vec, d, ~] = svd(S);

            stc_eig_val = diag(d);            
    end

end
        
% % 
% % switch nargout
% %     case 1  % STA only
% %         
% %     case 2  % STA and STC eigen value only (COVARIANCE NOT NEEDED)
% % 
% %         % subtract mean & scale down 
% %         %SS = bsxfun(@minus, SS, sta) / (num_bins-1);
% %         X = bsxfun(@minus, X, sta);       % DO NOT DIVIDE FOR BETTER NUMERICAL PRECISION
% % 
% %         [~, D, ~] = svd(X);
% %         stc_eig_val = diag(D).^2;
% % 
% %     case 3  % STA and STC eigen value & eiven vectors (COVARIANCE NOT NEEDED)
% % 
% %         % subtract mean & scale down 
% %         X = bsxfun(@minus, X, sta);       % DO NOT DIVIDE FOR BETTER NUMERICAL PRECISION
% % 
% %         [~, D, stc_eig_vec] = svd(X);
% %         stc_eig_val = diag(D).^2;
% % 
% %         
% %         
% %     case 4  % STA and full STC
% %         %% calc STC
% %         
% %         % subtract mean
% % %         X = bsxfun(@minus, X, sta);       % DO NOT DIVIDE FOR BETTER NUMERICAL PRECISION
% % %         
% % %         covariance = X'*X;
% %          
% %         % implementation 1 (when spikes are either 1 or 0, Matlab cov function is faster!
% %         %covariance = cov(X);
% % 
% %         % implementation 2
% %         covariance = X'*bsxfun(@times, X, spikes) - sta*sta'*num_total_spikes;  % DO NOT DIVIDE FOR BETTER NUMERICAL PRECISION
% %         %covariance = X'*bsxfun(@times, X, spikes)/(num_total_spikes-1) - sta*sta'*num_total_spikes/(num_total_spikes-1);
% % 
% %         % implementation 3
% %         %rowlen = size(SS,2);
% %         %stc = X'*(X.*repmat(spikes,1,rowlen))/(num_total_spikes-1) - sta*sta'*num_total_spikes/(num_total_spikes-1);
% %         
% %         
% %         [stc_eig_vec, d, ~] = svd(covariance);
% %     
% %         stc_eig_val = diag(d);
% %         %stc_eig_val = stc_eig_val(stc_eig_val>1e-5);
% %     
% % end


if dim ~= 1
    sta = reshape(sta, n, dim);
    
% %     if nargout > 1
% %         if dim ~= 1
% %             % unpack STC into nxn blocks
% %             for i=1:dim
% %                 for j=i:dim
% %                     STC{i,j} = stc(n*(i-1)+1:n*i, n*(j-1)+1:n*j);
% %                     STC{j,i} = STC{i,j};
% %                 end
% %             end
% %             stc= STC;
% %         end
% %     end
end


return 


function S= makeStimRows(Stim, n, flag)
%  S = makeStimRows(Stim, n, flag);
%
%  Converts raw stimulus to a matrix where each row is loaded with the full
%  space-time stimulus at a particular moment in time.  The resulting
%  matrix has length equal to the number of time points and width equal to
%  the (number of spatial dimensions) x (kernel size n).
%
%  Inputs: 
%   Stim = stimulus, first dimension is time, other dimensions are spatial
%          dimensions
%   n = size of temporal kernel; number of time samples to include in each
%       row of stimulus matrix.
%   flag (optional)
%        = 0, default behavior: padded w/ zeros at beginning so length of
%          output matrix matches size of Stim
%        = 1, no padding with zeros: length of S is length(Stim)-n+1.
%        = vector of indices, (e.g. indices of spikes times).  Return
%        a matrix containting only the spiking stimuli
%
%  Output: S = matrix where each row is the size of the linear kernel
%
%  Last updated:  6/30/2005, J. Pillow

% parse inputs
if nargin < 3
    flag = 0;
%     if nargin < 2
%         global n
%         if isempty(n)
%             error('ERROR -- makeStimRows:  n is undefined');
%         end
%     end
end

sz = size(Stim);
n2 = prod(sz(2:end));  % total dimensionality in spatial dimensions

% If necessary, convert Stim to a 2D matrix
if (n2 > sz(2))       % reshape to matrix if necessary
    Stim = reshape(Stim, sz(1), n2);
end

if flag == 0  % Compute with zero-padding. ----------------------------------
    S = zeros(sz(1), n2*n);
    for j=1:n2
        S(:,(n*(j-1)+1):(n*j)) = ...
            fliplr(toeplitz(Stim(:, j), [Stim(1,j) zeros(1,n-1)]));
    end
    
    
elseif (length(flag) == 1)  % compute only for stimuli at times >= n ----------
    S = zeros(sz(1)-n+1,n2*n);
    for j=1:n2
        S(:,(n*(j-1)+1):(n*j)) = ...
            fliplr(toeplitz(Stim(n:sz(1), j), Stim(n:-1:1,j)));
    end
    
    
else % compute for spiking stimuli --------------------------------------------
    if (min(flag) < 1) || (max(flag) > sz(1))
        error('makeStimRows:  3rd arg should be spike indices (vals are too high or too low): %d %d', ...
            min(flag), max(flag));
    end
    S = zeros(length(flag), n2*n);
    % Do for spinds < n
    nsp1 = length(flag(flag<n));
    for j = 1:nsp1
        S(j,:) = reshape([zeros(n-flag(j), n2); Stim(1:flag(j),:)], 1, n2*n);
    end
    for j = nsp1 +1:length(flag)
        S(j,:) = reshape(Stim(flag(j)-n+1:flag(j),:), 1, n2*n);
    end
end

return

%% to debug

[sta1, stc1] = calc_STA_and_STC(Stim,sp,n);  

clf
subplot(121)
plot(sta1)
subplot(122)
imshow(stc1)
%%

[sta2, stc2] = calc_STA_and_STC([Stim Stim],sp,n);  

clf
subplot(121)
plot(sta2)
subplot(122)
imshow(stc2)

