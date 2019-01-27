

chn_idx = 8
T = size(sta_all_channels{1},1);
gridT = (-T+1:0)/fps;


clf

subplot(211)
plot(gridT*1000, sta_all_channels{chn_idx})
set(gca, 'ylim', [0 1])
box off
xlabel('time to spike onset (ms)')


sig = std(sta_all_channels{chn_idx}(:));

subplot(212)

hist(sta_all_channels{chn_idx}(:), 200); box off
%set(gca,'yscale','log')
hold on
ylim = get(gca, 'ylim');

plot(0.5-2.58*sig*[1 1], ylim, '--r')
plot(0.5+2.58*sig*[1 1], ylim, '--r')


title(sprintf('histogram of STA (\\sigma=%.3f)',sig))


set(gcf, 'paperposition', [0 0 24 18])
set(gcf, 'papersize', [24 18])

saveas(gcf, sprintf('STA_histogram.png'))
saveas(gcf, sprintf('STA_histogram.pdf'))

