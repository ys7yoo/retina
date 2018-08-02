function plotSTA(STA, gridT, width, height)
%% plot STA
subplot(121)
imshow(STA')
xlabel('t')
ylabel('pixel')
title('STA')

STAvar = var(STA,1);
subplot(322)
imshow(reshape(STAvar,height,width)',[])   % stim is rotate 90 degree!
xlabel('y') %xlabel('x')         % swap x and y
ylabel('x') %ylabel('y')
title(sprintf('var STA across time'))
axis on xy

% mark the pixel with largest variance 
[XX,YY]= meshgrid(1:height);  % grid for ploting
[mm, maxIdx] = max(STAvar);

hold on;
%plot(XX(maxIdx),YY(maxIdx),'+r', 'markersize', 8, 'linewidth', 2)
plot(YY(maxIdx),XX(maxIdx),'+r', 'markersize', 8, 'linewidth', 2)

[mm, minIdx] = min(STAvar);
[XX,YY]= meshgrid(1:height);

hold on;
%plot(XX(minIdx),YY(minIdx),'xb', 'markersize', 8, 'linewidth',4)
plot(YY(minIdx),XX(minIdx),'xb', 'markersize', 8, 'linewidth',4)



subplot(324)
plot(gridT, STA(:,maxIdx),'r')
xlabel('t')
ylabel('STA')
title('STA for the pixel with the largest variance')
box off

subplot(326)
plot(gridT, STA(:,minIdx),'b')
xlabel('t')
ylabel('STA')
title('STA for the pixel with the largest variance')
box off
    
%saveas(gcf, sprintf('Num%d_STA_%s.pdf',NUM_EXP,channelName))
    