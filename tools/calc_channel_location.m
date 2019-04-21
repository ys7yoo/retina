function [x, y] = calc_channel_location(XX, YY, channel_no)
%%
% parse channel index

idx_col = floor(channel_no/10);
idx_row = mod(channel_no,10);

% mark channel 27
x = XX(idx_row,idx_col);
y = YY(idx_row,idx_col);
%text(, +0.5, sprintfc('%d',27), 'HorizontalAlignment','center')


return 

%% 
[x, y] = calc_channel_location(XX, YY, 27)