function [mask, XX, YY] = generate_ROI_mask_neighbors(W, H, xy, radius)

[XX, YY] = meshgrid(1:W, 1:H);
%mask = zeros(H,W);

xy = round(xy);

% circular mask
%mask = ((XX-xy(1)).^2 + (YY-xy(2)).^2) <= radius*radius;

mask = abs(XX-xy(1))<=radius & abs(YY-xy(2))<=radius;






return 


%%


mask = generate_ROI_mask_neighbors(8, 8, RFs{n}.mean, 1.5)
imshow(mask)
hold on
plot(RFs{n}.mean(1), RFs{n}.mean(2), 'ob')

axis xy

set(gcf, 'paperposition', [0 0 9 8])
set(gcf, 'papersize', [9 8])

