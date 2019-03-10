function [mask, XX, YY] =  generate_ROI_mask_from_MEA(xy, roi_size, width, height)

% generate 95% confidence interval ellipse mask

[XX, YY] = meshgrid(1:width, 1:height);
XX = XX(:);
YY = YY(:);
%R = chol(C);


mask = (XX-xy(1) <= roi_size) & (XX-xy(1) >= -roi_size) & (YY-xy(2) <= roi_size) & (YY-xy(2) >= -roi_size);

mask = reshape(mask, height, width);

return 