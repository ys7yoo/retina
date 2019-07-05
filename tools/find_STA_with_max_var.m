function [STA_max_var, idxMax] = find_STA_with_max_var(STA, gridT, width, height)

% calc vaciance of STA across time
STAvar = var(STA,1);

% find pixel with the maximum variance in STA
[~, idxMax] = max(STAvar);

% choose STA with maximum variance
STA_max_var = STA(:,idxMax);



%% plot STA
subplot(121)
imshow(STA')
xlabel('t')
ylabel('pixel')
title('STA')


subplot(222)
imshow(reshape(STAvar,height,width)',[])   % stim is rotated by 90 degree!
xlabel('y') %xlabel('x')         % swap x and y
ylabel('x') %ylabel('y')
title(sprintf('var STA across time'))
axis on xy


% mark the pixel with the largest variance 
[XX,YY]= meshgrid(1:height);  % grid for ploting

hold on;
%plot(XX(idxMax),YY(idxMax),'+r', 'markersize', 8, 'linewidth', 2)
plot(YY(idxMax),XX(idxMax),'+r', 'markersize', 8, 'linewidth', 2)


% % % mark the pixel with the smallest variance 
% % [mm, idxMin] = min(STAvar);
% % [XX,YY]= meshgrid(1:height);
% % 
% % hold on;
% % %plot(XX(idxMin),YY(idxMin),'xb', 'markersize', 8, 'linewidth',4)
% % plot(YY(idxMin),XX(idxMin),'xb', 'markersize', 8, 'linewidth',4)


subplot(224)
plot(gridT, STA_max_var,'r')
xlabel('t')
ylabel('STA')
title('STA for the pixel with the largest variance')
box off

% % % choose STA with minimum variance
% % STA_min_var = STA(:,idxMin);
% % 
% % subplot(326)
% % plot(gridT, STA_max_var,'b')
% % xlabel('t')
% % ylabel('STA')
% % title('STA for the pixel with the largest variance')
% % box off
    
%saveas(gcf, sprintf('Num%d_STA_%s.pdf',NUM_EXP,channelName))
    