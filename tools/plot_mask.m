function plot_mask(STA, height, width)

if nargin<2
    height = sqrt(num_pixels);
    width = sqrt(num_pixels);
end


%% plot thresholded mask for STA
T = size(STA,1);
sig = std(STA(:));



%% smooth STA for each slice 
STA_smoothed = smooth_STA_slice(STA, height, width);


%% 
r=floor(sqrt(T));
c=ceil(T/r);

figure(1)
clf; cnt = 1;
for t=1:T
    subplot(r,c,cnt)
    imshow(reshape(STA_smoothed(t,:,:),height,width)>0.5+2.58*sig)
    cnt = cnt + 1;
end

set(gcf, 'paperposition', [0 0 24 18])
set(gcf, 'papersize', [24 18])

saveas(gcf, sprintf('mask_pos.png'))
saveas(gcf, sprintf('mask_pos.pdf'))


figure(2)
clf; cnt = 1;
for t=1:T
    subplot(r,c,cnt)
    imshow(reshape(STA_smoothed(t,:,:),height,width)<0.5-2.58*sig)
    cnt = cnt + 1;
end

set(gcf, 'paperposition', [0 0 24 18])
set(gcf, 'papersize', [24 18])

saveas(gcf, sprintf('mask_neg.png'))
saveas(gcf, sprintf('mask_neg.pdf'))