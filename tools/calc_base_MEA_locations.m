function [XX, YY] = calc_base_MEA_locations(exp_param)

%%
space_in_pixel = exp_param.inter_electrode_space / exp_param.pixel_size;

if isfield(exp_param, 'num_electrodes_per_dim')
    gridX = (1:exp_param.num_electrodes_per_dim);
else
    gridX = (1:8);   % default 8x8 MEA
end
    
gridX = (gridX-0.5)*space_in_pixel;
gridY = fliplr(gridX);          % Y decreases from the largest! (top)

% gridX = gridX + offset_x;
% gridX = gridY + offset_y;


[XX, YY] = meshgrid(gridX, gridY);




return 


%% test
[XX, YY] = calc_base_MEA_locations(exp_param)


% plot
%clf

plot(XX(:), YY(:), 'ok')

hold on
% 
channel_no = 11:10:81;

for co=channel_no
    [x,y] = calc_channel_location(XX, YY, co);
    text(x, y+0.5, sprintfc('%d',co), 'HorizontalAlignment','center')
end


%text(XX(1,:), YY(1,:)+0.5, sprintfc('%d',channel_no), 'HorizontalAlignment','center')

axis xy
box off


% mark channel 27
text(XX(7,2), YY(7,2)+0.5, sprintfc('%d',27), 'HorizontalAlignment','center')

