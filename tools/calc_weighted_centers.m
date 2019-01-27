function [center_pos, cov_pos, num_pos_pixels, center_neg, cov_neg, num_neg_pixels] = calc_weighted_centers(slice, X, Y, threshold)

ref_value = 0.5;

positive_values=max(slice - ref_value-threshold,0);
negative_values=max(ref_value-slice-threshold,0);

num_pos_pixels =  sum(positive_values>0);
num_neg_pixels =  sum(negative_values>0);

[XX, YY] = meshgrid(1:X,1:Y);
% XX=XX(:);
YY=YY(:);


[center_pos, cov_pos] = calc_mean_and_cov([XX(:), YY(:)], positive_values(:)/sum(positive_values(:)));
[center_neg, cov_neg] = calc_mean_and_cov([XX(:), YY(:)], negative_values(:)/sum(negative_values(:)));



% % 
% % 
% % pos_mean_x = positive_values*XX(:)./sum(positive_values);
% % pos_mean_y = positive_values*YY(:)./sum(positive_values);
% % 
% % neg_mean_x = negative_values*XX(:)./sum(negative_values);
% % neg_mean_y = negative_values*YY(:)./sum(negative_values);
% % 
% % 
% % 
% % % calc variance
% % pos_mean_x_sqr = positive_values*XX(:).^2./sum(positive_values);
% % pos_mean_y_sqr = positive_values*YY(:).^2./sum(positive_values);
% % 
% % pos_sig_x = sqrt(pos_mean_x_sqr - pos_mean_x.^2);
% % pos_sig_y = sqrt(pos_mean_y_sqr - pos_mean_y.^2);
% % 
% % neg_mean_x_sqr = negative_values*XX(:).^2./sum(negative_values);
% % neg_mean_y_sqr = negative_values*YY(:).^2./sum(negative_values);
% % 
% % neg_sig_x = sqrt(neg_mean_x_sqr - neg_mean_x.^2);
% % neg_sig_y = sqrt(neg_mean_y_sqr - neg_mean_y.^2);


return 

