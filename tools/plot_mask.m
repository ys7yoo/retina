function plot_mask(STA, sig_smooth, height, width)

if nargin<3
    num_pixels = size(STA,2);
    height = sqrt(num_pixels);
    width = sqrt(num_pixels);
end






%% smooth STA for each slice 
if sig_smooth > 0
    STA_smoothed = smooth_STA_slice(STA, sig_smooth, height, width);
else
    STA_smoothed = STA;
end


%% plot thresholded mask for STA
T = size(STA,1);
sig = std(STA(:));

r=floor(sqrt(T));
c=ceil(T/r);


figure(1)
clf; cnt = 1;
for t=1:T
    subplot(r,c,cnt)
    
    imagesc(reshape(STA_smoothed(t,:,:),height,width), [0.2 0.8])
    colormap gray
    
    cnt = cnt + 1;
end

set(gcf, 'paperposition', [0 0 24 18])
set(gcf, 'papersize', [24 18])

saveas(gcf, sprintf('sta_slice.png'))
saveas(gcf, sprintf('sta_slice.pdf'))



figure(2)
clf; cnt = 1;
for t=1:T
    subplot(r,c,cnt)
    
    mask_on = reshape(STA_smoothed(t,:,:),height,width)>0.5+2.58*sig;
    mask_off = reshape(STA_smoothed(t,:,:),height,width)<0.5-2.58*sig;
    imshow(mask_on | mask_off)
    cnt = cnt + 1;
end

set(gcf, 'paperposition', [0 0 24 18])
set(gcf, 'papersize', [24 18])

saveas(gcf, sprintf('mask.png'))
saveas(gcf, sprintf('mask.pdf'))

% 
% figure(2)
% clf; cnt = 1;
% for t=1:T
%     subplot(r,c,cnt)
%     imshow(mask_off)
%     cnt = cnt + 1;
% end
% 
% set(gcf, 'paperposition', [0 0 24 18])
% set(gcf, 'papersize', [24 18])
% 
% saveas(gcf, sprintf('mask_neg.png'))
% saveas(gcf, sprintf('mask_neg.pdf'))